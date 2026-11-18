/*---------------------------------------------------------------------------*\
  DualSiteLangmuir.C
\*---------------------------------------------------------------------------*/

#include "DualSiteLangmuir.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace isotherms
{
    defineTypeNameAndDebug(DualSiteLangmuir, 0);
    addToRunTimeSelectionTable(isothermModel, DualSiteLangmuir, dictionary);
}
}


Foam::isotherms::DualSiteLangmuir::DualSiteLangmuir
(
    const word& specieName,
    const dictionary& dict
)
:
    isothermModel(specieName, dict),
    T0_(dict.get<scalar>("T0")),
    qm01_(dict.get<scalar>("qm01")),
    b01_(dict.get<scalar>("b01")),
    dH1_(dict.get<scalar>("deltaH1")),
    chi1_(dict.getOrDefault<scalar>("chi1", 0)),
    qm02_(dict.get<scalar>("qm02")),
    b02_(dict.get<scalar>("b02")),
    dH2_(dict.get<scalar>("deltaH2")),
    chi2_(dict.getOrDefault<scalar>("chi2", 0))
{}


Foam::scalar Foam::isotherms::DualSiteLangmuir::qStar
(
    const scalar p, const scalar T
) const
{
    const scalar Tc = max(T, scalar(1));
    return site(p, Tc, qm01_, b01_, dH1_, chi1_)
         + site(p, Tc, qm02_, b02_, dH2_, chi2_);
}


Foam::scalar Foam::isotherms::DualSiteLangmuir::dqStardp
(
    const scalar p, const scalar T
) const
{
    const scalar Tc = max(T, scalar(1));
    return dSitedp(p, Tc, qm01_, b01_, dH1_, chi1_)
         + dSitedp(p, Tc, qm02_, b02_, dH2_, chi2_);
}


Foam::scalar Foam::isotherms::DualSiteLangmuir::dqStardT
(
    const scalar p, const scalar T
) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar dT = max(1e-4*Tc, 1e-3);
    return (qStar(p, Tc + dT) - qStar(p, Tc - dT))/(2*dT);
}


// ************************************************************************* //
