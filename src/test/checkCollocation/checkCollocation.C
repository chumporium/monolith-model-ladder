/*---------------------------------------------------------------------------*\
  checkCollocation.C

  Kopling implisit (lduPrimitiveMeshAssembly) memasangkan muka antarmuka
  fluida ke-j dengan muka padatan ke-j MENURUT INDEKS
  (mappedWallPolyPatch::mapCollocatedFaces). Utilitas ini memeriksa apakah
  pasangan indeks itu benar-benar muka yang sama secara geometri, per prosesor.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption("fluidPatch", "word", "default fluid_to_solid");
    argList::addOption("solidPatch", "word", "default solid_to_fluid");
    #include "setRootCase.H"
    #include "createTime.H"

    const word fp(args.getOrDefault<word>("fluidPatch", "fluid_to_solid"));
    const word sp(args.getOrDefault<word>("solidPatch", "solid_to_fluid"));

    fvMesh fluid(IOobject("fluid", runTime.timeName(), runTime, IOobject::MUST_READ));
    fvMesh solid(IOobject("solid", runTime.timeName(), runTime, IOobject::MUST_READ));

    const fvPatch& pf = fluid.boundary()[fluid.boundaryMesh().findPatchID(fp)];
    const fvPatch& ps = solid.boundary()[solid.boundaryMesh().findPatchID(sp)];

    label nBad = 0;
    scalar dMax = 0, aRelMax = 0, hMin = GREAT;
    if (pf.size() != ps.size())
    {
        Pout<< "UKURAN BEDA " << pf.size() << " vs " << ps.size() << endl;
    }
    const label n = min(pf.size(), ps.size());
    for (label j = 0; j < n; ++j)
    {
        const scalar d = mag(pf.Cf()[j] - ps.Cf()[j]);
        const scalar h = Foam::sqrt(pf.magSf()[j]);
        hMin = min(hMin, h);
        dMax = max(dMax, d/h);
        aRelMax = max(aRelMax, mag(pf.magSf()[j] - ps.magSf()[j])/pf.magSf()[j]);
        if (d > 1e-6*h) ++nBad;
    }
    Pout<< "muka " << n << " : pasangan TIDAK kolokasi = " << nBad
        << " ; max |dCf|/h = " << dMax << " ; max |dA|/A = " << aRelMax << endl;

    reduce(nBad, sumOp<label>());
    Info<< "TOTAL pasangan tidak kolokasi: " << nBad << endl;
    return 0;
}
