bool closedVolume = p_rgh.needReference();
dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

{
    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
      + SmassFluid[i]                 // kopling dua-arah adsorpsi (Bab 5.1)
    );
    // poreSourceImplicit: linearisasi sumber gas pori thd tekanan kanal; nol saat konvergen (p = p*).
    if (GporeFluid.set(i))
    {
        p_rghDDtEqn += fvm::Sp(GporeFluid[i], p_rgh) - GporeFluid[i]*pStarPoreFluid[i];
    }

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        p_rghEqn.setReference
        (
            pRefCell,
            compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
        );

        p_rghEqn.solve
        (
            p_rgh.select
            (
                (
                    finalIter
                 && corr == nCorr-1
                 && nonOrth == nNonOrthCorr
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA + p_rghEqn.flux();

            // massSourceDiag: sumber massa yang benar-benar dipakai persamaan ini (sink > 0)
            SmassUsedTot[i] = fvc::domainIntegrate(SmassFluid[i]).value();
            if (GporeFluid.set(i))
            {
                SmassUsedTot[i] += fvc::domainIntegrate(GporeFluid[i]*(p_rgh - pStarPoreFluid[i])).value();
            }

            // Relaksasi p_rgh hanya pada iterasi luar yang hasilnya TIDAK dipakai penyelesaian species/energi
            // terakhir (REGISTER PB0-V2-2b..2d, 2026-10-04). prevIter p_rgh disimpan sekali per langkah
            // (storeOldFluidFields.H), jadi tiap iterasi ter-relaksasi hanya mencapai faktor x perubahan
            // langkah, tak peduli jumlah iterasi luar. YEqn/EEqn iterasi final memakai rho (termo, dari p
            // ter-relaksasi) dan phi (tak ter-relaksasi) dari iterasi SEBELUM final: keduanya tak memenuhi
            // kontinuitas. Pada evakuasi itu merusak neraca species (galat -28 %, Y 0,15 -> 0,05) dan
            // membuat suhu bergantung acuan entalpi (0,32 K). Iterasi final sendiri memang sudah tanpa
            // relaksasi (OpenFOAM), tetapi itu terlambat untuk Y dan he.
            // Bila interfaceTolerance aktif, final dipicu kontrol kopling (bukan hitungan), sehingga
            // "dua iterasi terakhir" tak diketahui di muka: relaksasi p_rgh dimatikan seluruhnya.
            if (!finalIter && oCorr < nOuterCorr - 2 && interfaceTol <= 0)
            {
                p_rgh.relax();
            }

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

}

pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    if (!compressible)
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
    else
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /compressibility;
        thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
        rho = thermo.rho();
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
}
else
{
    thermo.correctRho(psi*p - psip0,  rhoMin, rhoMax);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrors.H"

rho = thermo.rho();

// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
