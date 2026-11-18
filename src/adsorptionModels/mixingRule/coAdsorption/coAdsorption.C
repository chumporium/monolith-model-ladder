/*---------------------------------------------------------------------------*\
  coAdsorption.C

  Catatan dict (dictionary-native, mudah dibaca):

      mixingRule  coAdsorption;
      coAdsorptionCoeffs
      {
          enhancement
          {
              CO2 { H2O 2.0; }     // target { promotor gamma; ... }
          }
          pRef { H2O 1000; }       // [Pa]
      }
\*---------------------------------------------------------------------------*/

#include "coAdsorption.H"
#include "addToRunTimeSelectionTable.H"
#include "DynamicList.H"

namespace Foam
{
namespace mixingRules
{
    defineTypeNameAndDebug(coAdsorption, 0);
    addToRunTimeSelectionTable(mixingRule, coAdsorption, dictionary);
}
}


Foam::mixingRules::coAdsorption::coAdsorption(const dictionary& dict)
:
    mixingRule(dict),
    built_(false)
{
    const dictionary& cd = dict.subDict("coAdsorptionCoeffs");

    DynamicList<Tuple2<word, word>> pairs;
    DynamicList<scalar> gammas;

    const dictionary& enh = cd.subDict("enhancement");
    for (const word& target : enh.toc())
    {
        const dictionary& pd = enh.subDict(target);
        for (const word& promoter : pd.toc())
        {
            pairs.append(Tuple2<word, word>(target, promoter));
            gammas.append(pd.get<scalar>(promoter));
        }
    }
    pair_.transfer(pairs);
    gamma_.transfer(gammas);

    const dictionary& pr = cd.subDict("pRef");
    for (const word& k : pr.toc())
    {
        pRef_.insert(k, pr.get<scalar>(k));
    }
}


void Foam::mixingRules::coAdsorption::buildIndex
(
    const PtrList<isothermModel>& isos
) const
{
    const label np = pair_.size();
    targetIdx_.setSize(np, -1);
    promoterIdx_.setSize(np, -1);
    promoterPref_.setSize(np, GREAT);

    forAll(pair_, k)
    {
        forAll(isos, i)
        {
            if (isos[i].specieName() == pair_[k].first())  targetIdx_[k] = i;
            if (isos[i].specieName() == pair_[k].second()) promoterIdx_[k] = i;
        }
        auto it = pRef_.cfind(pair_[k].second());
        promoterPref_[k] = it.good() ? *it : GREAT;
    }
    built_ = true;
}


void Foam::mixingRules::coAdsorption::loadings
(
    const PtrList<isothermModel>& isos,
    const scalarList& p,
    const scalar T,
    scalarList& q
) const
{
    const label N = isos.size();

    // Loading murni
    forAll(isos, i)
    {
        q[i] = isos[i].qStar(p[i], T);
    }

    // Faktor enhancement multiplikatif per target
    scalarList factor(N, scalar(1));
    forAll(pair_, k)
    {
        const label it = targetIdx_[k];
        const label ip = promoterIdx_[k];
        if (it < 0 || ip < 0) continue;

        const scalar pp = max(p[ip], scalar(0));
        const scalar theta = pp/(promoterPref_[k] + pp);
        factor[it] *= (1 + gamma_[k]*theta);
    }

    forAll(isos, i)
    {
        q[i] = max(q[i]*factor[i], scalar(0));
    }
}


void Foam::mixingRules::coAdsorption::evaluate
(
    const PtrList<isothermModel>& isos,
    const scalarList& p,
    const scalar T,
    scalarList& q,
    scalarSquareMatrix& dqdp,
    scalarList& dqdT
) const
{
    if (!built_)
    {
        buildIndex(isos);
    }

    const label N = isos.size();

    loadings(isos, p, T, q);

    // Jacobian numerik terhadap p_j
    dqdp = Zero;
    scalarList qP(N), qM(N), pPert(p);
    for (label j = 0; j < N; ++j)
    {
        const scalar dpj = max(1e-6*max(p[j], scalar(1)), SMALL);

        pPert[j] = p[j] + dpj;
        loadings(isos, pPert, T, qP);

        pPert[j] = max(p[j] - dpj, scalar(0));
        loadings(isos, pPert, T, qM);

        pPert[j] = p[j];

        for (label i = 0; i < N; ++i)
        {
            dqdp(i, j) = (qP[i] - qM[i])/(2*dpj);
        }
    }

    // d q_i/d T |_p numerik
    const scalar dT = max(1e-4*T, 1e-3);
    loadings(isos, p, T + dT, qP);
    loadings(isos, p, T - dT, qM);
    forAll(isos, i)
    {
        dqdT[i] = (qP[i] - qM[i])/(2*dT);
    }
}


// ************************************************************************* //
