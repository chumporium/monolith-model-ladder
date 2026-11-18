/*---------------------------------------------------------------------------*\
  isothermModel.C
\*---------------------------------------------------------------------------*/

#include "isothermModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isothermModel, 0);
    defineRunTimeSelectionTable(isothermModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::isothermModel::isothermModel
(
    const word& specieName,
    const dictionary& dict
)
:
    specieName_(specieName),
    coeffs_(dict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::isothermModel> Foam::isothermModel::New
(
    const word& specieName,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting isotherm model " << modelType
        << " for species " << specieName << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "isothermModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<isothermModel>(ctorPtr(specieName, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::isothermModel::reducedSpreadingPressure
(
    const scalar p,
    const scalar T
) const
{
    // psi(p,T) = INT_0^p q(p')/p' dp'
    // Integran q/p' regular di p'->0 (q ~ b p' untuk Henry), tetapi 1/p'
    // singular: pakai substitusi u = ln(p') sehingga dp'/p' = du dan
    // integran menjadi q(e^u). Batas bawah dipotong di p_lo kecil dengan
    // kontribusi Henry analitik INT_0^{p_lo} q/p' dp' ~= (dq/dp|_0) * p_lo.
    if (p <= SMALL)
    {
        return 0;
    }

    const scalar pLo = max(p*1e-8, SMALL);

    // Kontribusi Henry [0, pLo]: q ~ (dq/dp|0) p'  => q/p' ~ const
    scalar psi = dqStardp(scalar(0), T)*pLo;

    // Kuadratur komposit Simpson pada u = ln(p'), [ln(pLo), ln(p)]
    const label n = 200;                 // genap
    const scalar uLo = Foam::log(pLo);
    const scalar uHi = Foam::log(p);
    const scalar du = (uHi - uLo)/n;

    auto integrand = [&](const scalar u) -> scalar
    {
        const scalar pp = Foam::exp(u);
        return qStar(pp, T);            // q dp'/p' = q du
    };

    scalar s = integrand(uLo) + integrand(uHi);
    for (label k = 1; k < n; ++k)
    {
        const scalar u = uLo + k*du;
        s += (k % 2 ? 4.0 : 2.0)*integrand(u);
    }
    psi += du/3.0*s;

    return psi;
}


Foam::scalar Foam::isothermModel::pressureFromPsi
(
    const scalar psi,
    const scalar T
) const
{
    if (psi <= SMALL)
    {
        return 0;
    }

    // Bisection pada p: reducedSpreadingPressure monoton naik dalam p.
    scalar pLo = SMALL;
    scalar pHi = 1.0;                     // 1 Pa awal

    // Perbesar pHi sampai psi(pHi) > psi (maks ~ 1e9 Pa)
    label guard = 0;
    while (reducedSpreadingPressure(pHi, T) < psi && guard++ < 60)
    {
        pHi *= 4.0;
    }

    for (label it = 0; it < 80; ++it)
    {
        const scalar pMid = 0.5*(pLo + pHi);
        const scalar fMid = reducedSpreadingPressure(pMid, T) - psi;
        if (fMid > 0)
        {
            pHi = pMid;
        }
        else
        {
            pLo = pMid;
        }
        if ((pHi - pLo) <= 1e-9*pHi)
        {
            break;
        }
    }

    return 0.5*(pLo + pHi);
}


// ************************************************************************* //
