/*---------------------------------------------------------------------------*\
  kineticClosure.C
\*---------------------------------------------------------------------------*/

#include "kineticClosure.H"

namespace Foam
{
    defineTypeNameAndDebug(kineticClosure, 0);
    defineRunTimeSelectionTable(kineticClosure, dictionary);
}


Foam::kineticClosure::kineticClosure
(
    const fvMesh& mesh,
    const word& specieName,
    const dictionary& dict
)
:
    mesh_(mesh),
    specieName_(specieName),
    coeffs_(dict)
{}


Foam::autoPtr<Foam::kineticClosure> Foam::kineticClosure::New
(
    const fvMesh& mesh,
    const word& specieName,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting kinetic closure " << modelType
        << " for species " << specieName << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "kineticClosure",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<kineticClosure>(ctorPtr(mesh, specieName, dict));
}


// ************************************************************************* //
