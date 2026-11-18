/*---------------------------------------------------------------------------*\
  mixingRule.H

  Base class pluggable untuk ATURAN CAMPURAN multikomponen (Bab 6.3 panduan).
  Menggabungkan isoterm-isoterm KOMPONEN MURNI menjadi loading campuran
  q_i*(p_1..p_N, T) beserta Jacobiannya.

  Antarmuka (per sel, satuan natural):
      evaluate(isos, p, T, q, dqdp, dqdT)

        isos  : daftar isoterm murni (satu per species)
        p     : tekanan parsial p_i [Pa]          (ukuran N)
        T     : suhu [K]
        q     : keluaran loading campuran [mol/kg] (ukuran N)
        dqdp  : keluaran d q_i/d p_j  [mol/kg/Pa]  (N x N)
        dqdT  : keluaran d q_i/d T |_p [mol/kg/K]  (ukuran N)

  Chain-rule terhadap konsentrasi pori C_j dan ketergantungan-T melalui p
  (p_j = C_j R T / M_j) ditangani oleh regionAdsorption, BUKAN di sini.

  Implementasi: ideal (non-kompetitif / single), extendedLangmuir, IAST,
  coAdsorption (sinergis, material-spesifik).
\*---------------------------------------------------------------------------*/

#ifndef mixingRule_H
#define mixingRule_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "PtrList.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "isothermModel.H"

namespace Foam
{

class mixingRule
{
protected:

    const dictionary coeffs_;

public:

    TypeName("mixingRule");

    declareRunTimeSelectionTable
    (
        autoPtr,
        mixingRule,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    explicit mixingRule(const dictionary& dict);
    mixingRule(const mixingRule&) = delete;

    static autoPtr<mixingRule> New(const dictionary& dict);

    virtual ~mixingRule() = default;


    //- Hitung loading campuran dan turunannya pada satu sel
    virtual void evaluate
    (
        const PtrList<isothermModel>& isos,
        const scalarList& p,
        const scalar T,
        scalarList& q,
        scalarSquareMatrix& dqdp,
        scalarList& dqdT
    ) const = 0;
};


} // End namespace Foam

#endif

// ************************************************************************* //
