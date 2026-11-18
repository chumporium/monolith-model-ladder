/*---------------------------------------------------------------------------*\
  LDF.C
\*---------------------------------------------------------------------------*/

#include "LDF.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace kineticClosures
{
    defineTypeNameAndDebug(LDF, 0);
    addToRunTimeSelectionTable(kineticClosure, LDF, dictionary);
}
}


Foam::kineticClosures::LDF::LDF
(
    const fvMesh& mesh,
    const word& specieName,
    const dictionary& dict
)
:
    kineticClosure(mesh, specieName, dict),
    k_(0)
{
    if (dict.found("kLDF"))
    {
        k_ = dict.get<scalar>("kLDF");
    }
    else
    {
        // k = 15 Dc / rc^2
        const scalar Dc = dict.get<scalar>("Dc");
        const scalar rc = dict.get<scalar>("rc");
        k_ = 15.0*Dc/sqr(rc);
    }

    Info<< "    LDF[" << specieName << "]: k = " << k_ << " 1/s" << endl;
}


void Foam::kineticClosures::LDF::update
(
    volScalarField& q,
    const volScalarField& qStar,
    const scalar dt
)
{
    const scalar kdt = k_*dt;
    const scalar denom = 1.0 + kdt;

    // Internal
    {
        scalarField& qc = q.primitiveFieldRef();
        const scalarField& qold = q.oldTime().primitiveField();
        const scalarField& qs = qStar.primitiveField();
        forAll(qc, c)
        {
            qc[c] = max((qold[c] + kdt*qs[c])/denom, scalar(0));
        }
    }
    // Boundaries
    forAll(q.boundaryField(), patchi)
    {
        scalarField& qp = q.boundaryFieldRef()[patchi];
        const scalarField& qold = q.oldTime().boundaryField()[patchi];
        const scalarField& qs = qStar.boundaryField()[patchi];
        forAll(qp, f)
        {
            qp[f] = max((qold[f] + kdt*qs[f])/denom, scalar(0));
        }
    }
}


// ************************************************************************* //
