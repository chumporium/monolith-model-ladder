/*---------------------------------------------------------------------------*\
  Toth.H

  Isoterm Toth bergantung-suhu (Bab 6.1 panduan). Loading [mol/kg adsorben
  aktif]; faktor fraksi aktif f_a diterapkan di level material (regionAdsorption).

      q*(p,T) = ns(T) u / (1 + u^t)^(1/t),   u = b(T) p
      ns(T)   = ns0 exp[ chi (1 - T/T0) ]
      b(T)    = b0  exp[ (dH0/(R T0)) (T0/T - 1) ]          [b0 dalam 1/Pa]
      t(T)    = t0 + alpha (1 - T0/T)

      dqStar/dp = ns b / (1 + u^t)^((t+1)/t)

  Entri dictionary:
      T0 ; ns0 ; b0 ; t0 ; deltaH0 ; chi ; alpha_t ;
\*---------------------------------------------------------------------------*/

#ifndef isotherm_Toth_H
#define isotherm_Toth_H

#include "isothermModel.H"

namespace Foam
{
namespace isotherms
{

class Toth
:
    public isothermModel
{
    scalar T0_, ns0_, b0_, t0_, deltaH0_, chi_, alpha_;

    inline scalar ns(const scalar T) const
    {
        return ns0_*Foam::exp(chi_*(1 - T/T0_));
    }
    inline scalar bAff(const scalar T) const
    {
        return b0_*Foam::exp(deltaH0_/(adsorption::Rgas*T0_)*(T0_/T - 1));
    }
    inline scalar tExp(const scalar T) const
    {
        return max(t0_ + alpha_*(1 - T0_/T), scalar(0.01));
    }


public:

    TypeName("Toth");

    Toth(const word& specieName, const dictionary& dict);

    virtual scalar qStar(const scalar p, const scalar T) const;
    virtual scalar dqStardp(const scalar p, const scalar T) const;
    virtual scalar dqStardT(const scalar p, const scalar T) const;

    //- Faktor termodinamika Darken, bentuk-dalam-q (tanpa inversi isoterm):
    //      Gamma = 1 + (b p)^t = 1/(1 - (q/ns)^t)
    //  Kedua bentuk terverifikasi identik s/d 1e-14 (CATATAN 23.42c).
    virtual bool hasThermodynamicFactor() const { return true; }
    virtual scalar thermodynamicFactor(const scalar q, const scalar T) const;

    //- Panas adsorpsi yang MEMANG dipakai isoterm ini di b(T). Menjadikannya
    //  sumber tunggal bagi `deltaH` source energi (lihat isothermModel.H).
    virtual scalar heatOfAdsorption() const { return deltaH0_; }
};

} // End namespace isotherms
} // End namespace Foam

#endif

// ************************************************************************* //
