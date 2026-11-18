/*---------------------------------------------------------------------------*\
  speciesCoupledMixedFvPatchScalarField.H

  BC mixed/Robin kopling species di interface fluid-solid (Bab 5.3 / 10.7).
  Menjamin kontinuitas NILAI dan FLUKS konsentrasi massa C [kg/m3] antara:
    - sisi fluid : variabel FRAKSI MASSA Y  (phase fluid),  C_f = rho_f Y
    - sisi solid : variabel KONSENTRASI pori C (phase solid)

  Karena kedua sisi memakai variabel berbeda, BC memanggil "mesin konversi"
  (rho) untuk menyamakannya pada state lokal:

      refValue (fluid)  = C_solid,internal / rho_fluid          (C -> Y)
      refValue (solid)  = Y_fluid,internal * rho_fluid          (Y -> C)
      valueFraction     = KDelta_nbr / (KDelta_nbr + KDelta_this)
      KDelta            = Deff * deltaCoeffs   (ruang konsentrasi)

  Faktor rho saling hapus pada valueFraction (kedua sisi diskalakan rho yang
  sama di interface tipis), sehingga valueFraction = rasio konduktansi
  konsentrasi murni.

  Deff (ruang konsentrasi):
    - fluid : D molekular species          (D)
    - solid : eps_p * Ds                   (D = eps_p, DField = Ds.<sp>)

  Entri dictionary:
      type        speciesCoupledMixed;
      phase       fluid;          // atau solid
      nbrField    C.CO2;          // nama field tetangga
      D           1.6e-5;         // fluid SAJA: D molekular [m2/s]
      DField      none;           // opsional nama field (solid: Ds.CO2 untuk Bosanquet)
      rho         rho;            // nama field densitas fluid (default "rho")
      value       uniform 0;

  !! SISI SOLID: `D` (= eps_p) SUDAH TIDAK DIPAKAI -- JANGAN DITULIS.
     eps_p bersumber TUNGGAL dari constant/<region>/adsorptionProperties
     (entri `epsilonP`). Sebelumnya ia hidup di DUA berkas tanpa pemeriksaan
     silang (CATATAN 23.80[1]): mengubah material di satu tempat saja memberi
     persamaan pori dan fluks antarmuka dgn eps_p BERBEDA -- ketidakcocokan yang
     baru terlihat setelah satu run, lewat neraca massa antarmuka.
     Perilaku sekarang:
       `D` absen            -> eps_p dari adsorptionProperties.            (BENAR)
       `D` ada & sama       -> peringatan USANG, nilai dict yang dipakai.
                               (case lama tetap berjalan identik -- provenans)
       `D` ada & MENYIMPANG -> FatalIOError, dgn kedua nilai dicetak.
     `write()` tidak lagi mengeluarkan `D` di sisi solid, supaya duplikatnya
     tidak lahir kembali di tiap direktori waktu.
\*---------------------------------------------------------------------------*/

#ifndef speciesCoupledMixedFvPatchScalarField_H
#define speciesCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class speciesCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    //- "fluid" (field = Y) atau "solid" (field = C)
    word phase_;

    //- Nama field tetangga
    word nbrFieldName_;

    //- Koefisien difusi skalar [m2/s] (ruang konsentrasi)
    scalar D_;

    //- Nama field difusi opsional (mis. Ds.<sp>); "none" jika tak dipakai
    word DFieldName_;

    //- Nama field densitas fluid untuk konversi
    word rhoName_;


    // Private Member Functions

        //- Sisi solid: ambil eps_p dari constant/<region>/adsorptionProperties
        //  dan periksa silang entri `D` yang usang. Lihat catatan di kepala
        //  berkas ini.
        static scalar readEpsilonP(const fvPatch&, const dictionary&);


public:

    TypeName("speciesCoupledMixed");

    // Constructors

        speciesCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        speciesCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        speciesCoupledMixedFvPatchScalarField
        (
            const speciesCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        speciesCoupledMixedFvPatchScalarField
        (
            const speciesCoupledMixedFvPatchScalarField&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new speciesCoupledMixedFvPatchScalarField(*this)
            );
        }

        speciesCoupledMixedFvPatchScalarField
        (
            const speciesCoupledMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new speciesCoupledMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Koefisien difusi efektif pada patch (ruang konsentrasi) [m2/s]
        tmp<scalarField> effectiveD() const;

        virtual void updateCoeffs();

        //- Kopling IMPLISIT (useImplicit true): tulis koefisien antar sel
        //  fluida-padatan langsung ke matriks rakitan (fvMatrix::addFvMatrix).
        //  Lihat catatan di .C.
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>& matrix,
            const label iMatrix,
            const direction cmpt
        );

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
