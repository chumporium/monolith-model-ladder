/*---------------------------------------------------------------------------*\
  energyCoupledImplicitFvPatchScalarField.H

  Kopling ENERGI padatan-gas implisit (2026-09-22, REGISTER cacat 3).

  Kopling energi partisi (BC suhu mixed + iterasi luar) konvergen sangat
  lambat bila waktu konduksi antarmuka (~0,03 s pada dinding 0,1 mm) jauh di
  bawah dt (~2 s): 30 iterasi luar masih meninggalkan 60% panas adsorpsi tak
  terkirim, dan sisa itu dibuang saat langkah maju. Jalur coupledEnergyField
  bawaan tak dapat dipakai: koefisien silang mixedEnergy::manipulateMatrix
  ditulis dalam satuan h region sendiri tanpa konversi Cp, dan lower() hanya
  ditulis bila matriks sudah asimetris (ENr1trefMono divergen, T 204-707 K).

  BC ini dipasang pada field BANTU `hCpl` (salinan he tiap region), bukan
  pada he: tipe BC he diturunkan otomatis dari BC suhu dan tak dapat diganti.
  Solver menyalin matriks EEqn/hEqn ke hCpl, merakit keduanya, lalu menyalin
  hasil kembali ke he (createEnergyCoupling.H).

  Muka antarmuka, suhu dilinearkan per sel sekitar iterasi kini (T*, h*):
      T = T* + (h - h*)/Cp                         (eksak utk hConst, Y tetap)
      K = kappaEff * deltaCoeffs,  G = K_P K_Q/(K_P + K_Q) |Sf|     [W/K]
      q (Q -> P) = G (T_Q - T_P)
  Baris P:  diag += G/Cp_P,  off(P,Q) -= G/Cp_Q,
            source += G [ (T*_Q - h*_Q/Cp_Q) - (T*_P - h*_P/Cp_P) ]
  Kedua sisi menulis baris masing-masing dgn G dan linearisasi YANG SAMA,
  jadi q(Q->P) = -q(P->Q) tepat di SETIAP iterasi, bukan hanya saat konvergen.

  Masukan dari region masing-masing (dikelola solver, terdaftar di db):
      T         suhu (thermo)
      kappaCpl  konduktivitas efektif [W/m/K]: gas kappaEff; padatan betav*kappa
      CpCpl     Cp per sel [J/kg/K]
      hCpl      h* = nilai field ini sebelum solve (= he)

  Batas: pasangan muka kolokasi berindeks sama dan di prosesor yang sama
  (sama dgn speciesCoupledMixed implisit; dekomposisi aksial).

  Contoh (dipasang solver, bukan ditulis user):
      fluid_to_solid { type energyCoupledImplicit; }
\*---------------------------------------------------------------------------*/

#ifndef energyCoupledImplicitFvPatchScalarField_H
#define energyCoupledImplicitFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

namespace Foam
{

class energyCoupledImplicitFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    TypeName("energyCoupledImplicit");

    energyCoupledImplicitFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    energyCoupledImplicitFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    energyCoupledImplicitFvPatchScalarField
    (
        const energyCoupledImplicitFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    energyCoupledImplicitFvPatchScalarField
    (
        const energyCoupledImplicitFvPatchScalarField&
    );

    energyCoupledImplicitFvPatchScalarField
    (
        const energyCoupledImplicitFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new energyCoupledImplicitFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new energyCoupledImplicitFvPatchScalarField(*this, iF)
        );
    }

    //- Konduktansi seri G [W/K] per muka (dipakai juga oleh neraca energi)
    tmp<scalarField> conductance() const;

    //- Tulis koefisien antar sel (useImplicit true)
    virtual void manipulateMatrix
    (
        fvMatrix<scalar>& matrix,
        const label iMatrix,
        const direction cmpt
    );
};

} // End namespace Foam

#endif
