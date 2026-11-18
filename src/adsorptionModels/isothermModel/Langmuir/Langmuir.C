/*---------------------------------------------------------------------------*\
  Langmuir.C
\*---------------------------------------------------------------------------*/

#include "Langmuir.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace isotherms
{
    defineTypeNameAndDebug(Langmuir, 0);
    addToRunTimeSelectionTable(isothermModel, Langmuir, dictionary);
}
}


Foam::isotherms::Langmuir::Langmuir(const word& specieName, const dictionary& dict)
:
    isothermModel(specieName, dict),
    T0_(dict.get<scalar>("T0")),
    qm0_(dict.get<scalar>("qm0")),
    b0_(dict.get<scalar>("b0")),
    deltaH0_(dict.get<scalar>("deltaH0")),
    chi_(dict.getOrDefault<scalar>("chi", 0))
{}


Foam::scalar Foam::isotherms::Langmuir::qSat(const scalar T) const
{
    return qm0_*Foam::exp(chi_*(1 - max(T, scalar(1))/T0_));
}


Foam::scalar Foam::isotherms::Langmuir::b(const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    return b0_*Foam::exp(deltaH0_/(adsorption::Rgas*T0_)*(T0_/Tc - 1));
}


Foam::scalar Foam::isotherms::Langmuir::qStar(const scalar p, const scalar T) const
{
    const scalar bp = b(T)*max(p, scalar(0));
    return qSat(T)*bp/(1 + bp);
}


Foam::scalar Foam::isotherms::Langmuir::dqStardp(const scalar p, const scalar T) const
{
    const scalar bb = b(T);
    const scalar bp = bb*max(p, scalar(0));
    return qSat(T)*bb/sqr(1 + bp);
}


Foam::scalar Foam::isotherms::Langmuir::dqStardT(const scalar p, const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar dT = max(1e-4*Tc, 1e-3);
    return (qStar(p, Tc + dT) - qStar(p, Tc - dT))/(2*dT);
}


// ************************************************************************* //
