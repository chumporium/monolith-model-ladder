/*---------------------------------------------------------------------------*\
  idealMixing.C
\*---------------------------------------------------------------------------*/

#include "idealMixing.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace mixingRules
{
    defineTypeNameAndDebug(idealMixing, 0);
    addToRunTimeSelectionTable(mixingRule, idealMixing, dictionary);
}
}


Foam::mixingRules::idealMixing::idealMixing(const dictionary& dict)
:
    mixingRule(dict)
{}


void Foam::mixingRules::idealMixing::evaluate
(
    const PtrList<isothermModel>& isos,
    const scalarList& p,
    const scalar T,
    scalarList& q,
    scalarSquareMatrix& dqdp,
    scalarList& dqdT
) const
{
    dqdp = Zero;
    forAll(isos, i)
    {
        q[i]       = isos[i].qStar(p[i], T);
        dqdp(i, i) = isos[i].dqStardp(p[i], T);
        dqdT[i]    = isos[i].dqStardT(p[i], T);
    }
}


// ************************************************************************* //
