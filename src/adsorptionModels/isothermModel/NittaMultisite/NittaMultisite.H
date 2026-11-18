/*---------------------------------------------------------------------------*\
  NittaMultisite.H

  Isoterm Langmuir MULTISITE (Nitta dkk., J. Chem. Eng. Japan 17:39, 1984),
  bentuk permukaan homogen, komponen murni. Dipakai Grande dkk. 2006 untuk
  propilena/propana pada zeolit 4A (validasi fisisorpsi monolith).

      theta = q/qmax ,   theta = K(T) p (1 - theta)^a
      K(T)  = Kref exp[ (dH/R) (1/T - 1/Tref) ]          [Kref dalam 1/Pa]

  a = jumlah situs yang ditempati satu molekul; a = 1 -> Langmuir biasa.

  theta IMPLISIT. F(theta) = theta - x (1-theta)^a , x = K p, monoton naik
  untuk a > 0 (dF/dtheta = 1 + a x (1-theta)^(a-1) > 0), jadi akarnya tunggal
  di [0,1). Diselesaikan Newton dgn pengaman bisection.

  Turunan analitik:
      dtheta/dx = (1-theta)^a / (1 + a x (1-theta)^(a-1))
      dq/dp     = qmax K dtheta/dx
      dq/dT     = qmax p (dK/dT) dtheta/dx ,   dK/dT = -K dH/(R T^2)
  Faktor termodinamika (dari ln x = ln theta - a ln(1-theta)):
      Gamma = dln p/dln q = 1 + a theta/(1-theta)

  ⚠️ SATUAN. Grande Tabel 1 memberi K dalam kPa^-1 (0,0288 @ 423 K). Solver
  memakai p dalam Pa, jadi Kref = 0,0288e-3 1/Pa. Kunci sengaja `Kref`, BUKAN
  `b0`, supaya konversi satuan tidak diam-diam terbawa dari isoterm lain.

  Entri dictionary:
      type      NittaMultisite;
      qmax      2.857;       // [mol/kg]
      a         2.612;       // [-]
      Kref      2.88e-5;     // [1/Pa] pada Tref
      Tref      423;         // [K]
      deltaH0   34210;       // [J/mol] positif = eksotermik
\*---------------------------------------------------------------------------*/

#ifndef isotherm_NittaMultisite_H
#define isotherm_NittaMultisite_H

#include "isothermModel.H"

namespace Foam
{
namespace isotherms
{

class NittaMultisite
:
    public isothermModel
{
    scalar qmax_, a_, Kref_, Tref_, deltaH0_;

    //- K(T) [1/Pa]
    scalar K(const scalar T) const;

    //- Akar theta dari theta = x (1-theta)^a
    scalar theta(const scalar x) const;

    //- dtheta/dx pada akar theta
    scalar dthetadx(const scalar x, const scalar th) const;

public:

    TypeName("NittaMultisite");

    NittaMultisite(const word& specieName, const dictionary& dict);

    virtual scalar qStar(const scalar p, const scalar T) const;
    virtual scalar dqStardp(const scalar p, const scalar T) const;
    virtual scalar dqStardT(const scalar p, const scalar T) const;

    //- Kapasitas jenuh = qmax (tak bergantung T pada bentuk ini)
    virtual scalar qSat(const scalar T) const { return qmax_; }

    virtual scalar heatOfAdsorption() const { return deltaH0_; }

    virtual bool hasThermodynamicFactor() const { return true; }

    virtual scalar thermodynamicFactor(const scalar q, const scalar T) const;
};

} // End namespace isotherms
} // End namespace Foam

#endif

// ************************************************************************* //
