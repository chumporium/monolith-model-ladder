/*---------------------------------------------------------------------------*\
  asymptoticMatched.C
\*---------------------------------------------------------------------------*/

#include "asymptoticMatched.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace isotherms
{
    defineTypeNameAndDebug(asymptoticMatched, 0);
    addToRunTimeSelectionTable(isothermModel, asymptoticMatched, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::isotherms::asymptoticMatched::asymptoticMatched
(
    const word& specieName,
    const dictionary& dict
)
:
    isothermModel(specieName, dict),
    base_(isothermModel::New(specieName, dict.subDict("base"))),
    K_(dict.get<scalar>("K")),
    n_(dict.get<scalar>("n"))
{
    if (K_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "K harus positif; terbaca " << K_ << nl
            << "K adalah kemiringan cabang linear q0 = K c [m3/kg]."
            << exit(FatalIOError);
    }

    // n negatif = MINIMUM halus, yaitu satu-satunya pemakaian yang dimaksud:
    // memotong isoterm cekung di konsentrasi rendah. n positif akan memberi
    // MAKSIMUM halus -- ia menaikkan loading, bukan menurunkannya, dan itu
    // hampir pasti bukan yang diinginkan pemakainya. Ditolak agar salah-pakai
    // berhenti di sini, bukan muncul sebagai hasil yang tampak wajar.
    if (n_ >= 0)
    {
        FatalIOErrorInFunction(dict)
            << "n harus NEGATIF; terbaca " << n_ << nl
            << "n < 0 memberi minimum halus (memotong isoterm di c rendah)." << nl
            << "n > 0 memberi maksimum halus, yang MENAIKKAN loading." << nl
            << "Nilai dapat diturunkan dari rasio pendekatan phi: n = 1/log2(phi)."
            << exit(FatalIOError);
    }

    Info<< "    asymptoticMatched[" << specieName << "]: K=" << K_
        << " m3/kg, n=" << n_ << ", base=" << base_->type() << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

inline Foam::scalar Foam::isotherms::asymptoticMatched::q0
(
    const scalar p,
    const scalar T
) const
{
    return K_*max(p, scalar(0))/(adsorption::Rgas*max(T, scalar(1)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::isotherms::asymptoticMatched::qStar
(
    const scalar p,
    const scalar T
) const
{
    const scalar a = q0(p, T);
    const scalar b = base_->qStar(p, T);

    // Cabang mana pun yang nol membuat minimum halus nol.
    if (a < SMALL || b < SMALL)
    {
        return 0;
    }

    // Difaktorkan terhadap cabang TERKECIL sebelum dipangkatkan.
    // Bentuk langsung (a^n + b^n)^(1/n) dengan n = -6 melimpah begitu salah
    // satu cabang kecil: a^-6 pada a ~ 1e-3 sudah 1e18. Dengan memfaktorkan m,
    // kedua suku jatuh di (0,1] sehingga tak pernah melimpah.
    const scalar m = min(a, b);
    const scalar S = Foam::pow(a/m, n_) + Foam::pow(b/m, n_);

    return m*Foam::pow(S, 1.0/n_);
}


Foam::scalar Foam::isotherms::asymptoticMatched::dqStardp
(
    const scalar p,
    const scalar T
) const
{
    const scalar a = q0(p, T);
    const scalar b = base_->qStar(p, T);

    if (a < SMALL || b < SMALL)
    {
        // Di p -> 0 cabang linear yang mengendalikan, dan kemiringannya K/(RT).
        return K_/(adsorption::Rgas*max(T, scalar(1)));
    }

    const scalar dadp = K_/(adsorption::Rgas*max(T, scalar(1)));
    const scalar dbdp = base_->dqStardp(p, T);

    // dq/dp = S^(1/n - 1) [ (a/m)^(n-1) da/dp + (b/m)^(n-1) db/dp ]
    // Faktor m saling menghapus persis; lihat catatan di qStar().
    const scalar m = min(a, b);
    const scalar S = Foam::pow(a/m, n_) + Foam::pow(b/m, n_);

    return Foam::pow(S, 1.0/n_ - 1.0)
         *(Foam::pow(a/m, n_ - 1.0)*dadp + Foam::pow(b/m, n_ - 1.0)*dbdp);
}


Foam::scalar Foam::isotherms::asymptoticMatched::dqStardT
(
    const scalar p,
    const scalar T
) const
{
    const scalar Tc = max(T, scalar(1));
    const scalar a = q0(p, Tc);
    const scalar b = base_->qStar(p, Tc);

    if (a < SMALL || b < SMALL)
    {
        return 0;
    }

    // q0 = K p/(R T)  ->  dq0/dT = -K p/(R T^2) = -a/T
    const scalar dadT = -a/Tc;
    const scalar dbdT = base_->dqStardT(p, Tc);

    const scalar m = min(a, b);
    const scalar S = Foam::pow(a/m, n_) + Foam::pow(b/m, n_);

    return Foam::pow(S, 1.0/n_ - 1.0)
         *(Foam::pow(a/m, n_ - 1.0)*dadT + Foam::pow(b/m, n_ - 1.0)*dbdT);
}


// ************************************************************************* //
