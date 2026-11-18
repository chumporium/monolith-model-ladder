/*---------------------------------------------------------------------------*\
  Sips.C
\*---------------------------------------------------------------------------*/

#include "Sips.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace isotherms
{
    defineTypeNameAndDebug(Sips, 0);
    addToRunTimeSelectionTable(isothermModel, Sips, dictionary);
}
}


Foam::isotherms::Sips::Sips(const word& specieName, const dictionary& dict)
:
    isothermModel(specieName, dict),
    T0_(dict.get<scalar>("T0")),
    qm0_(dict.get<scalar>("qm0")),
    b0_(dict.get<scalar>("b0")),
    n0_(dict.get<scalar>("n0")),
    deltaH0_(dict.get<scalar>("deltaH0")),
    chi_(dict.getOrDefault<scalar>("chi", 0)),
    alphaN_(dict.getOrDefault<scalar>("alpha_n", 0))
{}


Foam::scalar Foam::isotherms::Sips::qStar(const scalar p, const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar bp = bAff(Tc)*max(p, scalar(0));
    const scalar bpn = Foam::pow(bp + SMALL, nExp(Tc));
    return qm(Tc)*bpn/(1 + bpn);
}


Foam::scalar Foam::isotherms::Sips::dqStardp(const scalar p, const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar b  = bAff(Tc);
    const scalar n  = nExp(Tc);
    const scalar pp = max(p, scalar(0));
    const scalar bp = b*pp;
    const scalar bpn = Foam::pow(bp + SMALL, n);
    // d/dp [ bpn/(1+bpn) ] = n bpn / (p (1+bpn)^2)
    const scalar dpp = max(pp, SMALL);
    return qm(Tc)*n*bpn/(dpp*sqr(1 + bpn));
}


Foam::scalar Foam::isotherms::Sips::dqStardT(const scalar p, const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar dT = max(1e-4*Tc, 1e-3);
    return (qStar(p, Tc + dT) - qStar(p, Tc - dT))/(2*dT);
}


// ************************************************************************* //
