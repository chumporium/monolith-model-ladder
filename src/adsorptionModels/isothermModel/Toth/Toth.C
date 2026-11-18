/*---------------------------------------------------------------------------*\
  Toth.C
\*---------------------------------------------------------------------------*/

#include "Toth.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace isotherms
{
    defineTypeNameAndDebug(Toth, 0);
    addToRunTimeSelectionTable(isothermModel, Toth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::isotherms::Toth::Toth(const word& specieName, const dictionary& dict)
:
    isothermModel(specieName, dict),
    T0_(dict.get<scalar>("T0")),
    ns0_(dict.get<scalar>("ns0")),
    b0_(dict.get<scalar>("b0")),
    t0_(dict.get<scalar>("t0")),
    deltaH0_(dict.get<scalar>("deltaH0")),
    chi_(dict.getOrDefault<scalar>("chi", 0)),
    alpha_(dict.getOrDefault<scalar>("alpha_t", 0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::isotherms::Toth::qStar(const scalar p, const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar u  = bAff(Tc)*max(p, scalar(0));
    const scalar t  = tExp(Tc);
    const scalar denom = Foam::pow(1 + Foam::pow(u + SMALL, t), 1.0/t);
    return ns(Tc)*u/max(denom, SMALL);
}


Foam::scalar Foam::isotherms::Toth::dqStardp(const scalar p, const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar b  = bAff(Tc);
    const scalar u  = b*max(p, scalar(0));
    const scalar t  = tExp(Tc);
    const scalar denom = Foam::pow(1 + Foam::pow(u + SMALL, t), (t + 1.0)/t);
    return ns(Tc)*b/max(denom, SMALL);
}


Foam::scalar Foam::isotherms::Toth::thermodynamicFactor
(
    const scalar q,
    const scalar T
) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar nsT = ns(Tc);
    if (nsT <= SMALL)
    {
        return scalar(1);
    }

    // theta dikurung DI BAWAH 1: Gamma = 1/(1-theta^t) divergen di theta->1.
    // Divergensi itu FISIS (kapasitas jenuh), tapi tak boleh jadi inf/NaN di
    // solver. THETA_MAX 0.9999 -> Gamma <= ~1e4 utk t~0.25; jauh di atas
    // Gamma=14,7 pada 5,6% (CATATAN 23.42c) sehingga tak mengikat di rentang
    // operasi mana pun yang dipakai proyek ini.
    const scalar THETA_MAX = 0.9999;
    const scalar theta = min(max(q/nsT, scalar(0)), THETA_MAX);
    const scalar t = tExp(Tc);
    const scalar denom = 1.0 - Foam::pow(theta, t);

    // Gamma >= 1 selalu; klem melindungi dari round-off di theta->0.
    return max(1.0/max(denom, SMALL), scalar(1));
}


Foam::scalar Foam::isotherms::Toth::dqStardT(const scalar p, const scalar T) const
{
    // d q*/d T |_p secara numerik (central difference)
    const scalar Tc = max(T, scalar(1));
    const scalar dT = max(1e-4*Tc, 1e-3);
    return (qStar(p, Tc + dT) - qStar(p, Tc - dT))/(2*dT);
}


// ************************************************************************* //
