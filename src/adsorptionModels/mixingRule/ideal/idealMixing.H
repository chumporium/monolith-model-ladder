/*---------------------------------------------------------------------------*\
  idealMixing.H

  Aturan campuran non-kompetitif / single-component (Bab 6.3). Setiap species
  mengikuti isoterm murninya sendiri tanpa kompetisi situs:

      q_i* = isoterm_i.qStar(p_i, T)

  Dipakai untuk: single-component, atau multikomponen yang teradsorpsi pada
  situs berbeda (tak saling bersaing). Jacobian diagonal.
\*---------------------------------------------------------------------------*/

#ifndef mixing_idealMixing_H
#define mixing_idealMixing_H

#include "mixingRule.H"

namespace Foam
{
namespace mixingRules
{

class idealMixing
:
    public mixingRule
{
public:

    TypeName("ideal");

    explicit idealMixing(const dictionary& dict);

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
