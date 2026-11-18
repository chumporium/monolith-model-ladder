/*---------------------------------------------------------------------------*\
  extendedLangmuir.H

  Aturan campuran Extended Langmuir kompetitif (Bab 6.3). Mensyaratkan setiap
  isoterm murni berbentuk Langmuir (isLangmuirForm() == true):

      q_i* = qSat_i b_i p_i / (1 + sum_k b_k p_k)

  Jacobian d q_i/d p_j analitik; d q_i/d T |_p numerik (central difference).
\*---------------------------------------------------------------------------*/

#ifndef mixing_extendedLangmuir_H
#define mixing_extendedLangmuir_H

#include "mixingRule.H"

namespace Foam
{
namespace mixingRules
{

class extendedLangmuir
:
    public mixingRule
{
    //- q_i pada (p,T) tertentu (helper untuk FD terhadap T)
    static void loadings
    (
        const PtrList<isothermModel>& isos,
        const scalarList& p,
        const scalar T,
        scalarList& q
    );

public:

    TypeName("extendedLangmuir");

    explicit extendedLangmuir(const dictionary& dict);

    virtual void evaluate
    (
        const PtrList<isothermModel>& isos,
        const scalarList& p,
        const scalar T,
        scalarList& q,
        scalarSquareMatrix& dqdp,
        scalarList& dqdT
    ) const;
};

} // End namespace mixingRules
} // End namespace Foam

#endif

// ************************************************************************* //
