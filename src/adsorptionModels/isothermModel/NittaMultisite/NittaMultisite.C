/*---------------------------------------------------------------------------*\
  NittaMultisite.C
\*---------------------------------------------------------------------------*/

#include "NittaMultisite.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace isotherms
{
    defineTypeNameAndDebug(NittaMultisite, 0);
    addToRunTimeSelectionTable(isothermModel, NittaMultisite, dictionary);
}
}


Foam::isotherms::NittaMultisite::NittaMultisite
(
    const word& specieName,
    const dictionary& dict
)
:
    isothermModel(specieName, dict),
    qmax_(dict.get<scalar>("qmax")),
    a_(dict.get<scalar>("a")),
    Kref_(dict.get<scalar>("Kref")),
    Tref_(dict.get<scalar>("Tref")),
    deltaH0_(dict.get<scalar>("deltaH0"))
{
    if (!(qmax_ > 0) || !(a_ > 0) || !(Kref_ > 0) || !(Tref_ > 0))
    {
        FatalErrorInFunction
            << "NittaMultisite[" << specieName << "]: qmax, a, Kref, Tref"
            << " harus > 0 (qmax=" << qmax_ << ", a=" << a_ << ", Kref="
            << Kref_ << " 1/Pa, Tref=" << Tref_ << " K)." << exit(FatalError);
    }

    Info<< "    NittaMultisite[" << specieName << "]: qmax=" << qmax_
        << " mol/kg, a=" << a_ << ", Kref=" << Kref_ << " 1/Pa @ "
        << Tref_ << " K (= " << 1e3*Kref_ << " 1/kPa), deltaH0="
        << deltaH0_ << " J/mol" << endl;
}


Foam::scalar Foam::isotherms::NittaMultisite::K(const scalar T) const
{
    const scalar Tc = max(T, scalar(1));
    return Kref_*Foam::exp(deltaH0_/adsorption::Rgas*(1/Tc - 1/Tref_));
}


Foam::scalar Foam::isotherms::NittaMultisite::theta(const scalar x) const
{
    if (!(x > 0))
    {
        return 0;
    }

    // Kurung [lo,hi] selalu memuat akar: F(0) = -x < 0, F(1) = 1 > 0.
    scalar lo = 0, hi = 1;
    scalar th = x/(1 + x);                    // tebakan Langmuir

    for (label it = 0; it < 200; ++it)
    {
        const scalar om = 1 - th;
        const scalar F = th - x*Foam::pow(om, a_);

        if (F > 0) { hi = th; } else { lo = th; }

        const scalar dF = 1 + a_*x*Foam::pow(max(om, VSMALL), a_ - 1);
        scalar thN = th - F/dF;

        // Newton keluar kurung -> bisection
        if (!(thN > lo && thN < hi))
        {
            thN = 0.5*(lo + hi);
        }

        if (mag(thN - th) <= 1e-15 + 1e-13*thN)
        {
            return thN;
        }
        th = thN;
    }

    return th;
}


Foam::scalar Foam::isotherms::NittaMultisite::dthetadx
(
    const scalar x,
    const scalar th
) const
{
    const scalar om = max(1 - th, VSMALL);
    return Foam::pow(om, a_)/(1 + a_*x*Foam::pow(om, a_ - 1));
}


Foam::scalar Foam::isotherms::NittaMultisite::qStar
(
    const scalar p,
    const scalar T
) const
{
    return qmax_*theta(K(T)*max(p, scalar(0)));
}


Foam::scalar Foam::isotherms::NittaMultisite::dqStardp
(
    const scalar p,
    const scalar T
) const
{
    const scalar KT = K(T);
    const scalar x = KT*max(p, scalar(0));
    return qmax_*KT*dthetadx(x, theta(x));
}


Foam::scalar Foam::isotherms::NittaMultisite::dqStardT
(
    const scalar p,
    const scalar T
) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar KT = K(Tc);
    const scalar pp = max(p, scalar(0));
    const scalar x = KT*pp;
    const scalar dKdT = -KT*deltaH0_/(adsorption::Rgas*Tc*Tc);
    return qmax_*pp*dKdT*dthetadx(x, theta(x));
}


Foam::scalar Foam::isotherms::NittaMultisite::thermodynamicFactor
(
    const scalar q,
    const scalar T
) const
{
    const scalar th = min(max(q/qmax_, scalar(0)), 1 - 1e-12);
    return 1 + a_*th/(1 - th);
}


// ************************************************************************* //
