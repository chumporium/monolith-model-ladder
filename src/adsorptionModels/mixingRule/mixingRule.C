/*---------------------------------------------------------------------------*\
  mixingRule.C
\*---------------------------------------------------------------------------*/

#include "mixingRule.H"

namespace Foam
{
    defineTypeNameAndDebug(mixingRule, 0);
    defineRunTimeSelectionTable(mixingRule, dictionary);
}


Foam::mixingRule::mixingRule(const dictionary& dict)
:
    coeffs_(dict)
{}


Foam::autoPtr<Foam::mixingRule> Foam::mixingRule::New(const dictionary& dict)
{
    const word modelType(dict.getOrDefault<word>("mixingRule", "ideal"));

    Info<< "Selecting mixing rule " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "mixingRule",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<mixingRule>(ctorPtr(dict));
}


// ************************************************************************* //
