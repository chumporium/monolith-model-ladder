/*---------------------------------------------------------------------------*\
  effCondFoam  (Paper B, Tahap 1, 2026-10-04)

  Konduktivitas efektif DIAM (tanpa aliran) sel periodik dua-fasa, masalah sel homogenisasi:
      T = G.x + theta,  theta periodik,   div( k grad theta ) = - div( k G )
  pada mesh UTUH (sebelum splitMeshRegions) dengan cellZone gas dan padat, sisi cyclic.
  k muka = interpolasi harmonik (tepat untuk dua lapis seri). Hasil dari fluks menembus
  bidang batas yang normalnya searah G:
      k_eff = SUM_muka k_f ( G.S_f + snGrad(theta) |S_f| ) / ( |G| SUM S_f.e )
  Uji: lapis seri -> 1/(f1/k1 + f2/k2); lapis paralel -> f1 k1 + f2 k2.

  effCondFoam -kGas 0.026 -kSolid 1.5 -dir z [-zoneSolid solid] [-nCorr 6]
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "harmonic.H"

int main(int argc, char *argv[])
{
    argList::addOption("kGas", "scalar", "konduktivitas fasa gas [W/m/K]");
    argList::addOption("kSolid", "scalar", "konduktivitas fasa padat [W/m/K]");
    argList::addOption("dir", "word", "arah gradien: x | y | z");
    argList::addOption("zoneSolid", "word", "nama cellZone padat (bawaan solid)");
    argList::addOption("nCorr", "label", "jumlah koreksi non-ortogonal (bawaan 6)");
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    const scalar kG = args.get<scalar>("kGas");
    const scalar kS = args.get<scalar>("kSolid");
    const word dir = args.getOrDefault<word>("dir", "z");
    const word zS = args.getOrDefault<word>("zoneSolid", "solid");
    const label nCorr = args.getOrDefault<label>("nCorr", 6);
    const vector e(dir == "x" ? 1 : 0, dir == "y" ? 1 : 0, dir == "z" ? 1 : 0);
    const dimensionedVector G("G", dimTemperature/dimLength, e);

    const label zI = mesh.cellZones().findZoneID(zS);
    if (zI < 0) FatalErrorInFunction << "cellZone " << zS << " tidak ada" << exit(FatalError);

    volScalarField k
    (
        IOobject("kCond", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimPower/dimLength/dimTemperature, kG), "zeroGradient"
    );
    const labelList& cs = mesh.cellZones()[zI];
    forAll(cs, i) { k[cs[i]] = kS; }
    k.correctBoundaryConditions();
    const scalarField& V = mesh.V();
    const scalar Vt = gSum(V);
    scalar Vs = 0; forAll(cs, i) { Vs += V[cs[i]]; } reduce(Vs, sumOp<scalar>());

    const surfaceScalarField kf("kf", harmonic(mesh).interpolate(k));

    wordList types(mesh.boundary().size());
    forAll(types, pI)
    {
        const word& t = mesh.boundary()[pI].patch().type();
        types[pI] = polyPatch::constraintType(t) ? t : word("zeroGradient");
    }
    volScalarField theta
    (
        IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimTemperature, 0.0), types
    );
    const surfaceScalarField phiG("phiG", kf*(G & mesh.Sf()));

    dictionary sd;
    sd.add("solver", "PCG"); sd.add("preconditioner", "DIC");
    sd.add("tolerance", 1e-14); sd.add("relTol", 0.0); sd.add("maxIter", 20000);

    scalar kEff = 0;
    for (label nc = 0; nc <= nCorr; ++nc)
    {
        fvScalarMatrix A(fvm::laplacian(kf, theta) + fvc::div(phiG));
        A.setReference(0, 0.0);
        A.solve(sd);

        const surfaceScalarField q("q", phiG + kf*fvc::snGrad(theta)*mesh.magSf());
        scalar Q = 0, Ar = 0;
        forAll(mesh.boundary(), pI)
        {
            // patch dipilih dari normal RATA-RATA-nya (sisi cyclic boleh tidak rata: muka miring ikut);
            // luas = luas proyeksi SUM S_f.e (= luas penampang sel periodik)
            const fvPatch& fp = mesh.boundary()[pI];
            const vector Sp(gSum(fp.Sf()));
            const scalar Ap = gSum(fp.magSf());
            if (Ap < VSMALL || (Sp & e)/Ap < 0.5) continue;
            Q += sum(q.boundaryField()[pI]);
            Ar += sum(fp.Sf() & e);
        }
        reduce(Q, sumOp<scalar>()); reduce(Ar, sumOp<scalar>());
        kEff = Q/max(Ar, VSMALL);
        Info<< "koreksi " << nc << ": k_eff = " << kEff << " W/m/K (luas bidang " << Ar << " m2)" << endl;
    }

    const scalar fS = Vs/Vt;
    Info<< nl << "HASIL arah " << dir << ": k_eff = " << kEff << " W/m/K;  k_eff/k_s = " << kEff/kS << nl
        << "      fraksi padat = " << fS << ";  paralel " << fS*kS + (1 - fS)*kG
        << ";  seri " << 1.0/(fS/kS + (1 - fS)/kG) << nl << endl;
    theta.write();
    Info<< "End" << endl;
    return 0;
}
