/*---------------------------------------------------------------------------*\
  IAST.H

  Ideal Adsorbed Solution Theory (Bab 6.3). Aturan campuran kompetitif yang
  menerima isoterm murni APA PUN (termasuk Toth) — memakai reduced spreading
  pressure psi(p,T) = INT_0^p q/p' dp' yang disediakan isothermModel.

  Algoritma per sel (tekanan parsial gas p_i, suhu T):
    1. Cari psi bersama s.h.  sum_i p_i / p_i^0(psi) = 1,
       dengan p_i^0(psi) = isoterm_i.pressureFromPsi(psi, T)   (bisection).
    2. x_i = p_i / p_i^0(psi)                 (fraksi mol fasa teradsorpsi)
    3. 1/q_t = sum_i x_i / q_i^0,  q_i^0 = isoterm_i.qStar(p_i^0, T)
    4. q_i = x_i q_t

  Jacobian d q_i/d p_j dan d q_i/d T |_p dihitung numerik (central difference)
  karena IAST bersifat implisit.
\*---------------------------------------------------------------------------*/

#ifndef mixing_IAST_H
#define mixing_IAST_H

#include "mixingRule.H"

namespace Foam
{
namespace mixingRules
{

class IAST
:
    public mixingRule
{
    //- Toleransi & iterasi solver spreading-pressure
    const label maxIter_;
    const scalar tol_;

    //- Selesaikan loading campuran IAST pada satu sel
    void solveIAST
    (
        const PtrList<isothermModel>& isos,
        const scalarList& p,
        const scalar T,
        scalarList& q
    ) const;

public:

    TypeName("IAST");

    explicit IAST(const dictionary& dict);

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
