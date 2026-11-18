/*---------------------------------------------------------------------------*\
  localEquilibrium.C
\*---------------------------------------------------------------------------*/

#include "localEquilibrium.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace kineticClosures
{
    defineTypeNameAndDebug(localEquilibrium, 0);
    addToRunTimeSelectionTable(kineticClosure, localEquilibrium, dictionary);
}
}


Foam::kineticClosures::localEquilibrium::localEquilibrium
(
    const fvMesh& mesh,
    const word& specieName,
    const dictionary& dict
)
:
    kineticClosure(mesh, specieName, dict)
{}


void Foam::kineticClosures::localEquilibrium::update
(
    volScalarField& q,
    const volScalarField& qStar,
    const scalar dt
)
{
    // q = q* aljabar
    q = qStar;
    q.clamp_min(scalar(0));
}


// ************************************************************************* //
