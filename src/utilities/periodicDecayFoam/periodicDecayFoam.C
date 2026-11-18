/*---------------------------------------------------------------------------*\
  periodicDecayFoam  (Paper B, 2026-10-03)

  Mode berkembang-penuh yang meluruh ke arah aliran pada sel periodik:
      c(x) = exp(-lambda z) phi(x),  phi periodik (patch cyclic), phi = 0 di dinding,
  dengan medan kecepatan tunak U (dibaca dari waktu terakhir). Substitusi ke
  u.grad c = D lap c memberi masalah eigen
      -D lap phi + div(u phi) + 2 D lambda d(phi)/dz - D lambda^2 phi = lambda u_z phi.
  lambda terkecil dicari dengan iterasi invers (suku lambda di ruas kiri ditunda).
  Ini analog konsentrasi (atau temperatur) dinding konstan, Sh_T / Nu_T, yaitu kondisi
  yang di Paper A cocok dengan k_g dari front adsorpsi.

  Keluaran: lambda [1/m], <u_z> rata-rata volume gas, Sh = lambda <u_z> d_h^2 / (4 D),
  dan fluks dinding per volume gas untuk cek neraca.

  periodicDecayFoam -region fluid -D <m2/s> -dh <m> [-nIter 200] [-tol 1e-10]
  Patch bernama cyc* (cyclic) dibiarkan; patch lain: fixedValue 0.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    argList::addOption("D", "scalar", "difusivitas [m2/s]");
    argList::addOption("dh", "scalar", "diameter hidraulik acuan [m]");
    argList::addOption("nIter", "label", "iterasi maksimum (200)");
    argList::addOption("tol", "scalar", "toleransi relatif lambda (1e-10)");
    #include "addRegionOption.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const scalar Dv = args.get<scalar>("D");
    const scalar dh = args.get<scalar>("dh");
    const label nIter = args.getOrDefault<label>("nIter", 200);
    const scalar tol = args.getOrDefault<scalar>("tol", 1e-10);

    instantList ts = runTime.times();
    runTime.setTime(ts.last(), ts.size() - 1);
    Info<< "waktu U: " << runTime.timeName() << nl;

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    // Fluks volumetrik: pakai phi KONSERVATIF solver (fluks massa) dibagi rho muka bila tersedia.
    // fvc::flux(U) (interpolasi U) tidak bebas divergensi pada sel terdistorsi dan memberi
    // sumber/sink semu di suku konveksi (REGISTER PB0-SH2: +2,8 % pada gyroid, 0 pada kanal persegi).
    IOobject phiHeader("phi", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE);
    IOobject rhoHeader("rho", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE);
    const bool adaPhi =
        phiHeader.typeHeaderOk<surfaceScalarField>(true)
     && rhoHeader.typeHeaderOk<volScalarField>(true);
    tmp<surfaceScalarField> tphiU;
    if (adaPhi)
    {
        surfaceScalarField phiM(phiHeader, mesh);
        volScalarField rho(rhoHeader, mesh);
        if (phiM.dimensions() == dimMass/dimTime)
        {
            tphiU = tmp<surfaceScalarField>::New("phiU", phiM/fvc::interpolate(rho));
        }
        else
        {
            tphiU = tmp<surfaceScalarField>::New("phiU", phiM);
        }
        Info<< "fluks: phi konservatif solver" << nl;
    }
    else
    {
        tphiU = tmp<surfaceScalarField>::New("phiU", fvc::flux(U));
        Info<< "fluks: fvc::flux(U) (phi/rho tidak ada)" << nl;
    }
    const surfaceScalarField& phiU = tphiU();
    const volScalarField uz(U.component(vector::Z));
    const scalarField& V = mesh.V();

    const scalar Vg = gSum(V);
    const scalar uMean = gSum(uz.primitiveField()*V)/Vg;
    Info<< "V gas = " << Vg << " m3, <u_z> = " << uMean << " m/s, maks div(phiU) = "
        << gMax(mag(fvc::div(phiU))().primitiveField()) << " 1/s" << nl;

    wordList types(mesh.boundary().size());
    forAll(types, pI)
    {
        const word& t = mesh.boundary()[pI].patch().type();
        types[pI] = polyPatch::constraintType(t) ? t : word("fixedValue");
    }
    volScalarField phi
    (
        IOobject("phiDecay", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimless, 1.0), types
    );
    forAll(phi.boundaryField(), pI)
    {
        if (types[pI] == "fixedValue") phi.boundaryFieldRef()[pI] == 0.0;
    }
    volScalarField psi("psi", phi);

    const dimensionedScalar Dd("D", dimArea/dimTime, Dv);
    // tebakan awal: Sh = 3
    scalar lambda = 4*3.0*Dv/(uMean*dh*dh);
    const surfaceScalarField ezSf("ezSf", mesh.Sf() & vector(0, 0, 1));

    dictionary sd;
    sd.add("solver", "PBiCGStab"); sd.add("preconditioner", "DILU");
    sd.add("tolerance", 1e-14); sd.add("relTol", 0.0); sd.add("maxIter", 5000);

    scalar lamOld = lambda;
    for (label it = 0; it < nIter; ++it)
    {
        const surfaceScalarField phiTot
        (
            "phiTot",
            phiU + dimensionedScalar(dimVelocity, 2*Dv*lambda)*ezSf
        );
        fvScalarMatrix A
        (
            fvm::div(phiTot, psi, "div(phi,psi)")
          - fvm::laplacian(Dd, psi, "laplacian(D,psi)")
          - fvm::Sp(dimensionedScalar(dimless/dimTime, Dv*lambda*lambda), psi)
         ==
            dimensionedScalar(dimless/dimLength, lambda)*uz*phi
        );
        A.solve(sd);

        const scalar sPhi = gSum(uz.primitiveField()*phi.primitiveField()*V);
        const scalar sPsi = gSum(uz.primitiveField()*psi.primitiveField()*V);
        lamOld = lambda;
        lambda = lambda*sPhi/sPsi;
        phi.primitiveFieldRef() = psi.primitiveField()*(gSum(uz.primitiveField()*V)/sPsi);
        phi.correctBoundaryConditions();

        const scalar rel = mag(lambda - lamOld)/lambda;
        if (it % 5 == 0 || rel < tol)
        {
            Info<< "iter " << it << " lambda = " << lambda << " 1/m  Sh = "
                << lambda*uMean*dh*dh/(4*Dv) << "  rel " << rel << nl;
        }
        if (rel < tol && it > 3) break;
    }

    // cek neraca: fluks dinding = lambda * INT u_z phi dV (per satuan panjang periode, untuk mode)
    scalar Jw = 0;
    forAll(phi.boundaryField(), pI)
    {
        if (types[pI] != "fixedValue") continue;
        Jw += gSum(Dv*phi.boundaryField()[pI].snGrad()*mesh.magSf().boundaryField()[pI]);
    }
    const scalar src = lambda*gSum(uz.primitiveField()*phi.primitiveField()*V);
    const scalar cCup = gSum(uz.primitiveField()*phi.primitiveField()*V)/gSum(uz.primitiveField()*V);
    Info<< nl << "HASIL lambda = " << lambda << " 1/m" << nl
        << "      Sh (d_h = " << dh << ") = " << lambda*uMean*dh*dh/(4*Dv) << nl
        << "      fluks dinding -INT D dphi/dn dA = " << -Jw << " ; lambda INT u_z phi dV = " << src
        << " (rasio " << -Jw/src << ")" << nl
        << "      c_cup (phi) = " << cCup << nl << endl;

    runTime.writeNow();
    phi.write();
    Info<< "End" << endl;
    return 0;
}
