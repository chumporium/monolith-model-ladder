/*---------------------------------------------------------------------------*\
  extendedLangmuir.C
\*---------------------------------------------------------------------------*/

#include "extendedLangmuir.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace mixingRules
{
    defineTypeNameAndDebug(extendedLangmuir, 0);
    addToRunTimeSelectionTable(mixingRule, extendedLangmuir, dictionary);
}
}


Foam::mixingRules::extendedLangmuir::extendedLangmuir(const dictionary& dict)
:
    mixingRule(dict)
{}


void Foam::mixingRules::extendedLangmuir::loadings
(
    const PtrList<isothermModel>& isos,
    const scalarList& p,
    const scalar T,
    scalarList& q
)
{
    scalar D = 1;
    forAll(isos, k)
    {
        D += isos[k].b(T)*max(p[k], scalar(0));
    }
    forAll(isos, i)
    {
        q[i] = isos[i].qSat(T)*isos[i].b(T)*max(p[i], scalar(0))/D;
    }
}


void Foam::mixingRules::extendedLangmuir::evaluate
(
    const PtrList<isothermModel>& isos,
    const scalarList& p,
    const scalar T,
    scalarList& q,
    scalarSquareMatrix& dqdp,
    scalarList& dqdT
) const
{
    const label N = isos.size();

    // Validasi bentuk Langmuir
    forAll(isos, i)
    {
        if (!isos[i].isLangmuirForm())
        {
            FatalErrorInFunction
                << "extendedLangmuir membutuhkan isoterm murni berbentuk "
                << "Langmuir; species " << isos[i].specieName()
                << " tidak (gunakan IAST untuk isoterm umum)."
                << exit(FatalError);
        }
    }

    scalar D = 1;
    scalarList bp(N);
    forAll(isos, k)
    {
        bp[k] = isos[k].b(T)*max(p[k], scalar(0));
        D += bp[k];
    }

    // Loadings
    forAll(isos, i)
    {
        q[i] = isos[i].qSat(T)*bp[i]/D;
    }

    // d q_i/d p_j = qSat_i b_i ( delta_ij/D - p_i b_j/D^2 )
    for (label i = 0; i < N; ++i)
    {
        const scalar pref = isos[i].qSat(T)*isos[i].b(T);
        for (label j = 0; j < N; ++j)
        {
            const scalar delta = (i == j) ? 1.0 : 0.0;
            dqdp(i, j) =
                pref*(delta/D - max(p[i], scalar(0))*isos[j].b(T)/sqr(D));
        }
    }

    // d q_i/d T |_p numerik
    const scalar dT = max(1e-4*T, 1e-3);
    scalarList qP(N), qM(N);
    loadings(isos, p, T + dT, qP);
    loadings(isos, p, T - dT, qM);
    forAll(isos, i)
    {
        dqdT[i] = (qP[i] - qM[i])/(2*dT);
    }
}


// ************************************************************************* //
