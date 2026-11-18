/*---------------------------------------------------------------------------*\
  wallShapeFoam  (Paper B, 2026-10-03)

  Faktor bentuk difusi dinding untuk bentuk sembarang (padanan delta^2/3 slab):
      -lap(psi) = 1 di region padat, psi = 0 di antarmuka gas-dinding,
      patch constraint (cyclic, ...) dibiarkan, patch lain zeroGradient.
  <psi> = rata-rata volume psi [m2]. Koefisien LDF kuasi-tunak (profil parabolik) untuk
  serapan seragam:  k_pori = eps_p D_s / <psi>  [1/s], per volume dinding.
  Slab setebal delta (satu sisi terbuka): <psi> = delta^2/3. Silinder R: R^2/8. Bola R: R^2/15.

  wallShapeFoam -region solid [-iface solid_to_fluid]
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    argList::addOption("iface", "word", "nama patch antarmuka (bawaan: semua patch mapped/wall bernama *_to_*)");
    #include "addRegionOption.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word ifaceName = args.getOrDefault<word>("iface", "");
    wordList types(mesh.boundary().size());
    scalar Aif = 0;
    forAll(types, pI)
    {
        const fvPatch& fp = mesh.boundary()[pI];
        const word& t = fp.patch().type();
        const bool isIf = ifaceName.size() ? (fp.name() == ifaceName) : (fp.name().find("_to_") != std::string::npos);
        if (polyPatch::constraintType(t)) types[pI] = t;
        else if (isIf) { types[pI] = "fixedValue"; Aif += gSum(fp.magSf()); }
        else types[pI] = "zeroGradient";
    }
    volScalarField psi
    (
        IOobject("psiWall", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimArea, 0.0), types
    );
    dictionary sd;
    sd.add("solver", "PCG"); sd.add("preconditioner", "DIC");
    sd.add("tolerance", 1e-16); sd.add("relTol", 0.0); sd.add("maxIter", 20000);

    for (label nc = 0; nc < 4; ++nc)      // koreksi non-ortogonal
    {
        fvScalarMatrix A(-fvm::laplacian(psi) == dimensionedScalar(dimless, 1.0));
        A.solve(sd);
    }
    const scalarField& V = mesh.V();
    const scalar Vw = gSum(V);
    const scalar psiM = gSum(psi.primitiveField()*V)/Vw;
    Info<< nl << "HASIL V dinding = " << Vw << " m3, A antarmuka = " << Aif << " m2, V/A = " << Vw/Aif << " m" << nl
        << "      <psi> = " << psiM << " m2" << nl
        << "      <psi>/(V/A)^2 = " << psiM/sqr(Vw/Aif) << "   (slab 1/3; silinder 1/2; bola 3/5)" << nl << endl;
    psi.write();
    Info<< "End" << endl;
    return 0;
}
