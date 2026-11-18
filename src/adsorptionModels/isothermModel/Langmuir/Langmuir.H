/*---------------------------------------------------------------------------*\
  Langmuir.H

  Isoterm Langmuir bergantung-suhu (komponen murni). Bentuk Langmuir penuh
  sehingga dapat dipakai oleh aturan campuran extended Langmuir.

      q*(p,T) = qm(T) b(T) p / (1 + b(T) p)
      qm(T)   = qm0 exp[ chi (1 - T/T0) ]
      b(T)    = b0  exp[ (dH0/(R T0)) (T0/T - 1) ]    [b0 dalam 1/Pa]

  Entri dictionary: T0 ; qm0 ; b0 ; deltaH0 ; chi (opsional) ;
\*---------------------------------------------------------------------------*/

#ifndef isotherm_Langmuir_H
#define isotherm_Langmuir_H

#include "isothermModel.H"

namespace Foam
{
namespace isotherms
{

class Langmuir
:
    public isothermModel
{
    scalar T0_, qm0_, b0_, deltaH0_, chi_;

public:

    TypeName("Langmuir");

    Langmuir(const word& specieName, const dictionary& dict);

    virtual scalar qStar(const scalar p, const scalar T) const;
    virtual scalar dqStardp(const scalar p, const scalar T) const;
    virtual scalar dqStardT(const scalar p, const scalar T) const;

    virtual bool isLangmuirForm() const { return true; }
    virtual scalar qSat(const scalar T) const;
    virtual scalar b(const scalar T) const;

    //- Panas adsorpsi yang MEMANG dipakai isoterm ini di b(T). Menjadikannya
    //  sumber tunggal bagi `deltaH` source energi (lihat isothermModel.H).
    virtual scalar heatOfAdsorption() const { return deltaH0_; }
};

} // End namespace isotherms
} // End namespace Foam

#endif

// ************************************************************************* //
