/*---------------------------------------------------------------------------*\
  danckwertsInletFvPatchScalarField.H

  Syarat batas inlet DANCKWERTS untuk fraksi massa species (2026-09-18).

      fluks TOTAL masuk (konveksi + difusi) = fluks umpan:
      (-phi) Y_b + rho_b D (Y_b - Y_P) delta |Sf| = (-phi) Y_feed

  Dirichlet (fixedValue Y = Y_feed) memompa fluks difusif TAMBAHAN selama ada
  gradien di muka inlet. Pada Pe_L = uL/D ~ 3-5 (Grande 2006) itu 18% massa
  terserap dan breakthrough 18-33% terlalu awal. Danckwerts adalah BC fisik
  kolom terdispersi.

  Bentuk mixed: Y_b = f Y_feed + (1-f) Y_P, refGrad = 0,
      f = (-phi) / ((-phi) + rho_b D delta |Sf|)
  -> fluks masuk diskret = (-phi) Y_feed tepat. Muka keluar (phi > 0): f = 0
  (zeroGradient). D = difusivitas molekular [m2/s] (HARUS sama dgn DmolList
  species ini; bagian turbulen tidak dimasukkan -> hanya untuk aliran laminar).

  Contoh:
      inlet { type danckwertsInlet; feedValue 0.136411; D 1.068e-4;
              rho rho; phi phi; value uniform 0.136411; }
\*---------------------------------------------------------------------------*/

#ifndef danckwertsInletFvPatchScalarField_H
#define danckwertsInletFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class danckwertsInletFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    scalar feedValue_;
    scalar D_;
    word rhoName_;
    word phiName_;

public:

    TypeName("danckwertsInlet");

    danckwertsInletFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    danckwertsInletFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    danckwertsInletFvPatchScalarField
    (
        const danckwertsInletFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    danckwertsInletFvPatchScalarField
    (
        const danckwertsInletFvPatchScalarField&
    );

    danckwertsInletFvPatchScalarField
    (
        const danckwertsInletFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new danckwertsInletFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new danckwertsInletFvPatchScalarField(*this, iF)
        );
    }

    scalar feedValue() const noexcept { return feedValue_; }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
