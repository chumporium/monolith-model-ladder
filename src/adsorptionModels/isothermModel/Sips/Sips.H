/*---------------------------------------------------------------------------*\
  Sips.H

  Isoterm Sips (Langmuir-Freundlich) bergantung-suhu (komponen murni).

      q*(p,T) = qm(T) (b p)^n / (1 + (b p)^n)
      qm(T)   = qm0 exp[ chi (1 - T/T0) ]
      b(T)    = b0  exp[ (dH0/(R T0)) (T0/T - 1) ]    [b0 dalam 1/Pa]
      n(T)    = n0 + alpha_n (1 - T0/T)

  Entri dictionary: T0 ; qm0 ; b0 ; n0 ; deltaH0 ; chi (opsional) ; alpha_n (opsional) ;
\*---------------------------------------------------------------------------*/

#ifndef isotherm_Sips_H
#define isotherm_Sips_H

#include "isothermModel.H"

namespace Foam
{
namespace isotherms
{

class Sips
:
    public isothermModel
{
    scalar T0_, qm0_, b0_, n0_, deltaH0_, chi_, alphaN_;

    inline scalar qm(const scalar T) const
    {
        return qm0_*Foam::exp(chi_*(1 - T/T0_));
    }
    inline scalar bAff(const scalar T) const
    {
        return b0_*Foam::exp(deltaH0_/(adsorption::Rgas*T0_)*(T0_/T - 1));
    }
    inline scalar nExp(const scalar T) const
    {
        return max(n0_ + alphaN_*(1 - T0_/T), scalar(0.01));
    }

public:

    TypeName("Sips");

    Sips(const word& specieName, const dictionary& dict);

    virtual scalar qStar(const scalar p, const scalar T) const;
    virtual scalar dqStardp(const scalar p, const scalar T) const;
    virtual scalar dqStardT(const scalar p, const scalar T) const;

    //- Panas adsorpsi yang MEMANG dipakai isoterm ini di b(T). Menjadikannya
    //  sumber tunggal bagi `deltaH` source energi (lihat isothermModel.H).
    virtual scalar heatOfAdsorption() const { return deltaH0_; }
};

} // End namespace isotherms
} // End namespace Foam

#endif

// ************************************************************************* //
