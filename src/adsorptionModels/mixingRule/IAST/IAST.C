/*---------------------------------------------------------------------------*\
  IAST.C
\*---------------------------------------------------------------------------*/

#include "IAST.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace mixingRules
{
    defineTypeNameAndDebug(IAST, 0);
    addToRunTimeSelectionTable(mixingRule, IAST, dictionary);
}
}


Foam::mixingRules::IAST::IAST(const dictionary& dict)
:
    mixingRule(dict),
    maxIter_(dict.getOrDefault<label>("iastMaxIter", 100)),
    tol_(dict.getOrDefault<scalar>("iastTol", 1e-8))
{}


void Foam::mixingRules::IAST::solveIAST
(
    const PtrList<isothermModel>& isos,
    const scalarList& p,
    const scalar T,
    scalarList& q
) const
{
    const label N = isos.size();
    q = Zero;

    // Tekanan total teradsorpsi
    scalar pSum = 0;
    forAll(p, i)
    {
        pSum += max(p[i], scalar(0));
    }
    if (pSum <= SMALL)
    {
        return;
    }

    // g(psi) = sum_i p_i / p_i^0(psi) - 1   (monoton turun dalam psi)
    auto gFun = [&](const scalar psi) -> scalar
    {
        scalar s = 0;
        forAll(isos, i)
        {
            if (p[i] > SMALL)
            {
                const scalar p0 = max(isos[i].pressureFromPsi(psi, T), SMALL);
                s += p[i]/p0;
            }
        }
        return s - 1;
    };

    // Bracket: psiLo kecil (g>0), psiHi diperbesar sampai g<0
    scalar psiLo = SMALL;
    scalar psiHi = 1.0;
    label guard = 0;
    while (gFun(psiHi) > 0 && guard++ < 80)
    {
        psiHi *= 2.0;
    }

    // Bisection
    scalar psi = 0.5*(psiLo + psiHi);
    for (label it = 0; it < maxIter_; ++it)
    {
        psi = 0.5*(psiLo + psiHi);
        const scalar g = gFun(psi);
        if (g > 0)
        {
            psiLo = psi;
        }
        else
        {
            psiHi = psi;
        }
        if (mag(g) < tol_ || (psiHi - psiLo) <= tol_*psiHi)
        {
            break;
        }
    }

    // x_i dan q_t
    scalarList x(N, Zero);
    scalar invqt = 0;
    forAll(isos, i)
    {
        if (p[i] > SMALL)
        {
            const scalar p0 = max(isos[i].pressureFromPsi(psi, T), SMALL);
            x[i] = p[i]/p0;
            const scalar q0 = max(isos[i].qStar(p0, T), SMALL);
            invqt += x[i]/q0;
        }
    }
    const scalar qt = (invqt > SMALL) ? 1.0/invqt : 0.0;

    forAll(isos, i)
    {
        q[i] = x[i]*qt;
    }
}


void Foam::mixingRules::IAST::evaluate
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

    solveIAST(isos, p, T, q);

    // Jacobian numerik (central difference) terhadap p_j
    dqdp = Zero;
    scalarList qP(N), qM(N), pPert(p);
    for (label j = 0; j < N; ++j)
    {
        const scalar dpj = max(1e-6*max(p[j], scalar(1)), SMALL);

        pPert[j] = p[j] + dpj;
        solveIAST(isos, pPert, T, qP);

        pPert[j] = max(p[j] - dpj, scalar(0));
        solveIAST(isos, pPert, T, qM);

        pPert[j] = p[j];

        for (label i = 0; i < N; ++i)
        {
            dqdp(i, j) = (qP[i] - qM[i])/(2*dpj);
        }
    }

    // d q_i/d T |_p numerik
    const scalar dT = max(1e-4*T, 1e-3);
    solveIAST(isos, p, T + dT, qP);
    solveIAST(isos, p, T - dT, qM);
    forAll(isos, i)
    {
        dqdT[i] = (qP[i] - qM[i])/(2*dT);
    }
}


// ************************************************************************* //
