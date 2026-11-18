/*---------------------------------------------------------------------------*\
  coAdsorption.H

  Aturan campuran KO-ADSORPSI SINERGIS (Bab 6.3). Untuk sistem non-kompetitif/
  sinergis di mana satu species menaikkan kapasitas species lain (mis. H2O
  menaikkan kapasitas CO2 pada amine). Bentuk fungsional material-spesifik;
  kelas ini menyediakan satu bentuk pengali yang umum dipakai sebagai TEMPLATE:

      q_i = q_i^murni(p_i,T) * ( 1 + sum_k gamma_ik * theta_k ),
      theta_k = p_k / (pRef_k + p_k)        (cakupan promotor tersaturasi)

  gamma_ik > 0  : sinergis (promotor k menaikkan kapasitas i)
  gamma_ik < 0  : antagonis

  Parameter dibaca dari sub-dictionary coAdsorptionCoeffs:
      enhancement ( (CO2 H2O 2.0) ... );   // (target promotor gamma)
      pRef        ( (H2O 1000) ... );       // [Pa] saturasi promotor

  Jacobian dihitung numerik (central difference). Untuk material nyata, ganti
  bentuk ini dengan korelasi dari data eksperimen lembap/simulasi molekular.
\*---------------------------------------------------------------------------*/

#ifndef mixing_coAdsorption_H
#define mixing_coAdsorption_H

#include "mixingRule.H"
#include "HashTable.H"
#include "Tuple2.H"

namespace Foam
{
namespace mixingRules
{

class coAdsorption
:
    public mixingRule
{
    //- enhancement: (targetName, promoterName, gamma)
    List<Tuple2<word, word>> pair_;
    scalarList gamma_;

    //- pRef per promotor [Pa]
    HashTable<scalar> pRef_;

    //- cache index (dibangun saat panggilan pertama)
    mutable labelList targetIdx_;
    mutable labelList promoterIdx_;
    mutable scalarList promoterPref_;   // sejajar dengan promoterIdx_
    mutable bool built_;

    void buildIndex(const PtrList<isothermModel>& isos) const;

    void loadings
    (
        const PtrList<isothermModel>& isos,
        const scalarList& p,
        const scalar T,
        scalarList& q
    ) const;

public:

    TypeName("coAdsorption");

    explicit coAdsorption(const dictionary& dict);

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
