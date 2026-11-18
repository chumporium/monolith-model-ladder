/*---------------------------------------------------------------------------*\
  danckwertsInletTemperatureFvPatchScalarField.H

  Syarat batas inlet DANCKWERTS untuk SUHU (2026-09-23) -- pasangan energi
  dari danckwertsInlet (species).

      fluks entalpi TOTAL masuk (konveksi + konduksi) = fluks umpan:
      (-phi) Cp_b (T_b - T_feed) + kappa_b (T_b - T_P) delta |Sf| = 0

  Dirichlet (fixedValue T = T_feed) mengizinkan konduksi keluar lewat inlet
  selama ada gradien. Pe termal Grande = u L rho cp / kappa ~ 1,3: pada F2b
  (G56_423, 0-3000 s) konduksi inlet membuang 30% panas adsorpsi. Danckwerts
  adalah BC energi baku kolom terdispersi (Ruthven), tanpa parameter bebas.

  Bentuk mixed: T_b = f T_feed + (1-f) T_P, refGrad = 0,
      f = (-phi) Cp_b / ((-phi) Cp_b + kappa_b delta |Sf|)
  Muka keluar (phi > 0): f = 0 (zeroGradient). Cp_b dan kappa_b dibaca dari
  thermo region (basicThermo::lookupThermo): tak ada sifat yang ditulis dua
  kali. kappa laminar (kappaEff turbulen tidak dimasukkan).

  Contoh:
      inlet { type danckwertsInletTemperature; feedValue 423; phi phi;
              value uniform 423; }
\*---------------------------------------------------------------------------*/

#ifndef danckwertsInletTemperatureFvPatchScalarField_H
#define danckwertsInletTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class danckwertsInletTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    scalar feedValue_;
    word phiName_;

    //- Campuran T dihitung di sini (T_b = f T_umpan + (1-f) T_P) dan dipasang dgn fraksi nilai 1 (V2-SF1e).
    //  Tanpa ini mixed-T dikonversi ke mixed-h, yang mencampur h sel (komposisi Y_P) dan h umpan
    //  (komposisi Y_b) -> suhu muka bergantung pada acuan entalpi. Bawaan false (kasus lama bit-identik).
    bool explicitBlend_;

public:

    TypeName("danckwertsInletTemperature");

    danckwertsInletTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    danckwertsInletTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    danckwertsInletTemperatureFvPatchScalarField
    (
        const danckwertsInletTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    danckwertsInletTemperatureFvPatchScalarField
    (
        const danckwertsInletTemperatureFvPatchScalarField&
    );

    danckwertsInletTemperatureFvPatchScalarField
    (
        const danckwertsInletTemperatureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new danckwertsInletTemperatureFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new danckwertsInletTemperatureFvPatchScalarField(*this, iF)
        );
    }

    scalar feedValue() const noexcept { return feedValue_; }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
