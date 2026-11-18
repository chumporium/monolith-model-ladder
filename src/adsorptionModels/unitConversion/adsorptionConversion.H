/*---------------------------------------------------------------------------*\
  adsorptionConversion.H

  "Mesin konversi satuan" (Bab 10.5 panduan adsorpsi_v2).

  Representasi internal kanonik transport fluida adalah FRAKSI MASSA Y_i
  (native thermo multikomponen kompresibel OpenFOAM). Model fisika (isoterm,
  kinetik) ditulis dalam SATUAN NATURALNYA: isoterm dalam tekanan parsial p_i
  [Pa] dan loading [mol/kg]; kinetik dalam [mol/kg] dan laju [1/s].

  Header ini adalah satu-satunya tempat konversi dilakukan ("satu sumber
  kebenaran"). Semua fungsi adalah inline & alokasi-bebas agar aman dipakai di
  hot-loop (wajah evaluasi).

  Rumus inti (Bab 10.5):
      X_i  = (Y_i/M_i) / sum_j (Y_j/M_j)          fraksi mol
      p_i  = X_i * p                              tekanan parsial [Pa]
      c_i  = rho * Y_i                            konsentrasi massa [kg/m3]
      p_i  = c_i * R * T / M_i                    (dari konsentrasi massa pori)
      Y_i  = c_i / rho                            (balik ke fraksi massa)
      ppm  = X_i * 1e6                            (molar/volume)
      mdot = M_i * rho_s * dq_i/dt                sink massa species [kg/m3/s]

  Konvensi satuan:
      p  [Pa]      T [K]     M [kg/mol]
      c  [kg/m3]   q [mol/kg]   R = 8.314462618 J/(mol K)
\*---------------------------------------------------------------------------*/

#ifndef adsorptionConversion_H
#define adsorptionConversion_H

#include "scalar.H"
#include "scalarList.H"

namespace Foam
{
namespace adsorption
{

//- Universal gas constant [J/(mol K)]
constexpr scalar Rgas = 8.314462618;

//- Mass concentration [kg/m3] of a species from mass fraction and density
inline scalar concFromY(const scalar rho, const scalar Yi)
{
    return rho*max(Yi, scalar(0));
}

//- Mass fraction from mass concentration [kg/m3] and density (state-dependent)
inline scalar yFromConc(const scalar rho, const scalar ci)
{
    return (rho > SMALL) ? max(ci, scalar(0))/rho : scalar(0);
}

//- Partial pressure [Pa] from pore mass concentration [kg/m3]:
//      p_i = (c_i / M_i) * R * T
inline scalar partialPressureFromConc
(
    const scalar ci,
    const scalar Mi,
    const scalar T
)
{
    return max(ci, scalar(0))*Rgas*max(T, scalar(1))/Mi;
}

//- Pore mass concentration [kg/m3] from partial pressure [Pa]:
//      c_i = p_i * M_i / (R T)
inline scalar concFromPartialPressure
(
    const scalar pi,
    const scalar Mi,
    const scalar T
)
{
    return max(pi, scalar(0))*Mi/(Rgas*max(T, scalar(1)));
}

//- Mole fractions X_i from mass fractions Y_i and molar weights M_i
inline void moleFractions
(
    const scalarList& Y,
    const scalarList& M,
    scalarList& X
)
{
    scalar sum = 0;
    forAll(Y, i)
    {
        X[i] = max(Y[i], scalar(0))/M[i];
        sum += X[i];
    }
    if (sum > SMALL)
    {
        forAll(X, i)
        {
            X[i] /= sum;
        }
    }
}

//- Loading converted from [mol/kg] to [kg/kg]
inline scalar loadingMassBasis(const scalar q_molkg, const scalar Mi)
{
    return q_molkg*Mi;
}

//- ppm (molar/volume) from mole fraction
inline scalar ppmFromX(const scalar Xi)
{
    return Xi*1e6;
}

} // End namespace adsorption
} // End namespace Foam

#endif

// ************************************************************************* //
