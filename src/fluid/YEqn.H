/*---------------------------------------------------------------------------*\
  YEqn.H

  Transport adsorbat sebagai FRAKSI MASSA Y_i (Bab 5.1, bentuk konservatif):

      d(rho Y_i)/dt + div(phi Y_i) - div(rho D_i grad Y_i) = R_i

  D_i = difusivitas molekular (DmolList) + bagian turbulen (mut/Sct).
  Konsentrasi C_i = rho Y_i terikat EOS -> merespons swing tekanan (TVSA).
  Penyerapan adsorbat ke solid terjadi via BC interface (speciesCoupledMixed),
  bukan source volumetrik. Kehilangan massa bulk masuk ke kontinuitas via
  SmassFluid (kopling dua-arah, fluidAdsorptionSource.H).
\*---------------------------------------------------------------------------*/

tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

if (Y.size())
{
    mvConvection = tmp<fv::convectionScheme<scalar>>
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
}

{
    reaction.correct();
    Qdot = reaction.Qdot();
    volScalarField Yt
    (
        IOobject("Yt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Yt", dimless, 0)
    );

    forAll(Y, j)
    {
        if (j != inertIndex && composition.active(j))
        {
            volScalarField& Yi = Y[j];

            // speciesCoupling implicit: adsorbat diselesaikan BERSAMA C padatan
            // di solveSolidSpecies.H (speciesImplicitCouplerFluid.H).
            if (speciesCouplerPtr)
            {
                bool implisit = false;
                forAll(Yi.boundaryField(), patchi)
                {
                    if (Yi.boundaryField()[patchi].useImplicit())
                    {
                        implisit = true;
                    }
                }
                if (implisit)
                {
                    Yt += Yi;
                    continue;
                }
            }

            // Difusivitas efektif: molekular + turbulen
            const dimensionedScalar Dmol
            (
                "Dmol", dimArea/dimTime, DmolList[i][j]
            );
            volScalarField DEff
            (
                "DEff." + Yi.name(),
                rho*Dmol + turbulence.mut()/SctFluid
            );

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(DEff, Yi)
             ==
                reaction.R(Yi)
              + fvOptions(rho, Yi)
            );

            if (SadvFluid[i].size() == Y.size() && SadvFluid[i].set(j))
            {
                YiEqn -= SadvFluid[i][j];      // + sumber adveksi pori di ruas kanan (tahap V2)
            }

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].clamp_min(0);
    }
}

// Kopling dua-arah: sink massa adsorpsi interface -> kontinuitas/pEqn
if (twoWayCoupling)
{
    #include "fluidAdsorptionSource.H"
}
