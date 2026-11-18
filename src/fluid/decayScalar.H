/*---------------------------------------------------------------------------*\
  decayScalar.H  (Paper B, Tahap 1, 2026-10-08)

  Skalar modus-peluruhan untuk Sh/Nu pada aliran periodik TAK TUNAK.

  periodicDecayFoam menyelesaikan modus peluruhan ruang c = exp(-lambda z) phi pada medan aliran
  BEKU; itu hanya sah untuk aliran tunak. Di sini modus yang sama diselesaikan BERSAMA aliran:

      d(rho phi)/dt + div(rho u phi) - div(rho D grad phi)
          = lambda rho u_z phi - 2 rho D lambda dphi/dz + rho D lambda^2 phi ,
      phi = 0 di dinding (penyerap sempurna), periodik di sisi cyclic,

  dengan lambda TETAP (diberikan). Tiap langkah phi dinormalkan ke rata-rata volume 1 dan laju tumbuh
  sesaat sigma = ln(g)/dt dicatat (g = rata-rata volume sebelum normalisasi). Untuk modus peluruhan
  ruang yang stasioner secara statistik, rata-rata waktu sigma harus nol. Beberapa skalar dengan
  lambda berbeda dijalankan serentak; lambda* diperoleh dari akar <sigma>(lambda) = 0, lalu
      Sh = lambda* <u> d_h^2 / (4 D).
  Pada aliran tunak hasilnya harus sama dengan periodicDecayFoam (gerbang verifikasi).

  fvSolution (region fluida), PIMPLE:
      decayScalars { sh1 { D 1.4e-5; lambda 150; } sh2 { D 1.4e-5; lambda 170; } ... }
  fvSchemes: div(phi,phiDec).   Keluaran: postProcessing/decayScalar.<region>.dat
  Tanpa entri decayScalars: tidak berbuat apa-apa.
\*---------------------------------------------------------------------------*/
if (finalIter && pimple.found("decayScalars"))
{
    const dictionary& dsDict = pimple.subDict("decayScalars");
    const wordList dsNames(dsDict.toc());

    if (decFields[i].size() != dsNames.size())
    {
        decFields[i].setSize(dsNames.size());
        wordList types(mesh.boundary().size());
        forAll(types, pI)
        {
            const fvPatch& fp = mesh.boundary()[pI];
            if (polyPatch::constraintType(fp.patch().type())) types[pI] = fp.patch().type();
            else if (isA<mappedPatchBase>(fp.patch())) types[pI] = fixedValueFvPatchScalarField::typeName;
            else types[pI] = zeroGradientFvPatchScalarField::typeName;
        }
        forAll(dsNames, k)
        {
            decFields[i].set
            (
                k,
                new volScalarField
                (
                    IOobject
                    (
                        "phiDec." + dsNames[k], runTime.timeName(), mesh,
                        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
                    ),
                    mesh, dimensionedScalar(dimless, 1.0), types
                )
            );
            volScalarField& f = decFields[i][k];
            forAll(f.boundaryField(), pI)
            {
                if (isA<fixedValueFvPatchScalarField>(f.boundaryField()[pI]))
                {
                    f.boundaryFieldRef()[pI] == 0.0;
                }
            }
            f.oldTime();
        }
    }

    const scalarField& Vc = mesh.V();
    const scalar Vt = gSum(Vc);
    const volScalarField uz("uzDec", U.component(vector::Z));
    const surfaceScalarField rhoSz("rhoSzDec", fvc::interpolate(rho)*(vector(0, 0, 1) & mesh.Sf()));   // berorientasi, seperti phi

    forAll(dsNames, k)
    {
        const dictionary& d = dsDict.subDict(dsNames[k]);
        const dimensionedScalar Dd("D", dimArea/dimTime, d.get<scalar>("D"));
        const dimensionedScalar lam("lambda", dimless/dimLength, d.get<scalar>("lambda"));
        volScalarField& f = decFields[i][k];

        // -2 rho D lambda dphi/dz  ->  fluks konvektif tambahan 2 D lambda rho (e_z . S_f)
        const surfaceScalarField phiEff("phiDecEff", phi + 2.0*Dd*lam*rhoSz);

        fvScalarMatrix fEqn
        (
            fvm::ddt(rho, f)
          + fvm::div(phiEff, f, "div(phi,phiDec)")
          - fvm::laplacian(rho*Dd, f, "laplacian(D,phiDec)")
         ==
            fvm::SuSp(lam*rho*uz + rho*Dd*sqr(lam), f)
        );
        dictionary sd;
        sd.add("solver", "PBiCGStab"); sd.add("preconditioner", "DILU");
        sd.add("tolerance", 1e-12); sd.add("relTol", 1e-6);
        fEqn.solve(sd);

        const scalar g = gSum(f.primitiveField()*Vc)/Vt;
        const scalar sigma = Foam::log(max(g, VSMALL))/runTime.deltaTValue();
        f.primitiveFieldRef() /= g;
        f.correctBoundaryConditions();
        f.oldTime().primitiveFieldRef() /= g;       // skema backward: jaga nisbah antar-tingkat waktu
        f.oldTime().correctBoundaryConditions();

        if (Pstream::master())
        {
            const fileName dir(runTime.globalPath()/"postProcessing");
            mkDir(dir);
            OFstream os
            (
                dir/("decayScalar." + mesh.name() + ".dat"),
                IOstreamOption(), IOstreamOption::APPEND
            );
            os  << runTime.value() << ' ' << dsNames[k] << ' ' << lam.value() << ' '
                << sigma << ' ' << g << nl;
        }
    }
}
