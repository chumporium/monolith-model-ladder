/*---------------------------------------------------------------------------*\
  regionAdsorption.H

  Koordinator adsorpsi untuk SATU solid region. Mengikat komponen pluggable
  (isothermModel x mixingRule x kineticClosure) + mesin konversi, lalu
  menyelesaikan persamaan species pori dan menyediakan source energi S_ads.

  Persamaan species pori (Bab 5.2 / 10.6), per adsorbat i:
      eps_p dC_i/dt + rho_s d(M_i q_i)/dt = div(eps_p Ds_i grad C_i)

  Mode local-equilibrium (closure aljabar): q_i = q_i*, storage transien
  dirakit via Jacobian:
      effStorage_i = eps_p + rho_s (dqStar_i/dp_i) R T         (diagonal)
      + source eksplisit silang (dqStar_i/dC_j) dan termal (dqStar_i/dT).

  Mode finite-rate (LDF/dualSite): q_i diintegrasi closure, species memakai
  source eksplisit -rho_s M_i dq_i/dt.

  Source energi solid (Bab 5.2):
      S_ads = sum_i deltaH_i rho_s dq_i/dt           [W/m3]

  Variabel:
      C_i [kg/m3] konsentrasi massa pori   (field, MUST_READ)
      q_i [mol/kg] loading                 (field, AUTO_WRITE)
  Isoterm dievaluasi dalam tekanan parsial p_i = C_i R T / M_i (mesin konversi).
\*---------------------------------------------------------------------------*/

#ifndef regionAdsorption_H
#define regionAdsorption_H

#include "fvCFD.H"
#include "solidThermo.H"
#include "isothermModel.H"
#include "mixingRule.H"
#include "kineticClosure.H"
#include "adsorptionConversion.H"

namespace Foam
{

// Forward declaration
class fvMesh;

/*---------------------------------------------------------------------------*                  Class speciesImplicitCoupler (antarmuka)
\*---------------------------------------------------------------------------*/
//- Kait kopling species IMPLISIT fluida-padatan (2026-09-17).
//  Bila diberikan ke regionAdsorption::correct, persamaan C_i yang sudah
//  dilinearisasi Newton TIDAK diselesaikan sendiri, melainkan diserahkan ke
//  kait ini untuk dirakit bersama Y_i fluida dalam satu matriks
//  (speciesCoupledMixed::manipulateMatrix). Tanpa kait -> jalur lama identik.
class speciesImplicitCoupler
{
public:

    virtual ~speciesImplicitCoupler() = default;

    //- Species padatan i dirakit implisit?
    virtual bool coupled(const label i) const = 0;

    //- Rakit C_i (padatan) + Y_i (fluida) lalu selesaikan bersama
    virtual void solve(const label i, fvScalarMatrix& CiEqn) = 0;

    //- Redaman line-search yang sama untuk Y_i fluida
    virtual void damp(const label i, const scalar alpha) = 0;
};


/*---------------------------------------------------------------------------*\
                     Class regionAdsorption Declaration
\*---------------------------------------------------------------------------*/

class regionAdsorption
{
    // Private data

        const fvMesh& mesh_;
        const solidThermo& thermo_;

        IOdictionary dict_;

        //- Daftar species adsorbat
        wordList species_;

        //- Berat molekul [kg/mol]
        scalarList M_;

        //- Panas adsorpsi [J/mol] (positif = eksotermik)
        scalarList deltaH_;

        //- Difusivitas pori referensi [m2/s]
        scalarList DsRef_;

        //- Bosanquet (Ds bergantung tekanan DAN suhu) per species
        boolList   bosanquet_;
        scalarList Dbulk_;     // D_bulk pada (pRef,TRef) [m2/s]
        scalarList pRefB_;     // tekanan referensi [Pa]
        scalarList TRefB_;     // suhu referensi [K]
        scalarList Dknud_;     // D_Knudsen [m2/s] (dipakai; terhitung bila rPore ada)
        boolList   DknudGiven_; // `Dknudsen` benar-benar DITULIS di dict?
                                // Tanpa ini, penjaga mencetak nilai BAWAAN
                                // seolah-olah itu angka yang diketik orang.
        scalarList rPore_;     // jari-jari pori [m]; <=0 = tidak diberikan

        //- Properti pori/padatan
        //  DUA porositas, karena dinding wash-coat punya DUA skala pori:
        //    epsilonP_    makropori -- rongga yang MENGANGKUT (fluks difusi)
        //    epsilonPtot_ pori TOTAL -- rongga yang MENYIMPAN gas (makro +
        //                 pori di dalam partikel wash-coat)
        //  Sampai 2026-09-04 keduanya dilebur jadi satu angka. Untuk monolith
        //  CACE angkanya 48,0% vs 53,1% -- simpanan gas pori kurang 5 poin
        //  volume kalau dipakai nilai transport. Tak ditulis -> epsilonPtot_
        //  = epsilonP_, jadi seluruh case lama BIT-IDENTIK.
        scalar epsilonP_;
        scalar epsilonPtot_;
        scalar rhoS_;
        scalar fActive_;       // f_a fraksi massa adsorben aktif

        //- Matikan panas adsorpsi: deltaH dipaksa 0 -> S_ads == 0.
        //  Untuk menjalankan kasus ISOTERMAL secara sengaja (mis. verifikasi
        //  silang terhadap model yang memang isotermal). Ia juga satu-satunya
        //  cara sah membuat `deltaH` menyimpang dari isoterm tanpa FatalError.
        bool isothermalAdsorption_;

        //- Kontrol iterasi Newton storage-konservatif (jalur algebraic)
        label  nNewton_;       // batas iterasi Newton per panggilan correct()
        scalar newtonTol_;     // toleransi relatif update C (konvergensi)

        //- Komponen pluggable
        PtrList<isothermModel> isotherms_;
        autoPtr<mixingRule>    mixing_;
        PtrList<kineticClosure> kinetics_;

        //- Field
        PtrList<volScalarField> C_;       // konsentrasi pori [kg/m3]
        PtrList<volScalarField> q_;       // loading [mol/kg]
        PtrList<volScalarField> qStar_;   // loading kesetimbangan [mol/kg]
        PtrList<volScalarField> Ds_;      // difusivitas pori [m2/s]

        //- Simpanan efektif  eps_p + rho_s M dq*/dC   [-]
        //  Disimpan supaya bilangan kopling Lambda dapat dihitung oleh
        //  pengendali langkah waktu. Persamaan padatan berbentuk
        //      effStor dC/dt = div(eps_p Ds grad C)
        //  jadi difusivitas EFEKTIF untuk C adalah eps_p Ds / effStor, dan
        //  itulah yang menentukan batas dt kopling -- bukan Ds telanjang.
        PtrList<volScalarField> effStor_;

        volScalarField Sads_;             // source energi [W/m3]

        //- Fluks volumetrik Darcy di pori [m3/s] (Paper B, tahap V2). nullptr = tanpa adveksi
        //  (semua kasus lama). Diset solver lewat setPoreFlux() bila `poreAdvection yes`.
        const surfaceScalarField* poreFlux_ = nullptr;
        //- Koefisien pembuangan [1/s] di sel antarmuka: adveksi KELUAR lewat antarmuka memakai C sel
        //  (upwind), bukan nilai batas BC kopling. nullptr = tak dipakai.
        const volScalarField* poreSink_ = nullptr;
        //- Sumber eksplisit per species [kg/m3/s] di sel antarmuka: gas kanal yang MASUK pori membawa
        //  komposisi sel kanal (hulu). nullptr / tak diset = tak dipakai.
        const PtrList<volScalarField>* poreSrc_ = nullptr;


    // Private Member Functions

        //- Perbarui Ds(p,T) bila Bosanquet aktif
        void updateDs();

        //- P1: rho_s dict vs thermo. Satu besaran, dua berkas.
        void checkSolidDensity() const;

        //- Tetapkan deltaH_[i] dari SUMBER TUNGGAL, dgn pemeriksaan silang.
        //  Kelas jebakan yang sama dgn rho_s dan eps_p: panas adsorpsi hidup
        //  di `deltaH` (source energi) dan `isotherm.deltaH0` (ket. suhu b(T)).
        //  Panggil SETELAH isotherms_[i] terbentuk. Lihat regionAdsorption.C.
        void resolveDeltaH(const label i, const dictionary& sd);

        //- P1: D_Knudsen terhitung vs DsRef pada (pRef,TRef). Melaporkan
        //  tortuositas TERSIRAT -- pengukuran, bukan sekadar kebersihan kode.
        void checkBosanquet();

        //- Inisialisasi q dari isoterm bila tidak dibaca dari disk
        void initialiseLoading();

        regionAdsorption(const regionAdsorption&) = delete;
        void operator=(const regionAdsorption&) = delete;


public:

    // Constructors
    regionAdsorption(const fvMesh& mesh, const solidThermo& thermo);


    // Member Functions

        //- Simpan oldTime field state (q dan state closure) — panggil sekali
        //  per langkah waktu sebelum loop PIMPLE
        void storeOldTime();

        //- Selesaikan species pori + isoterm + closure + hitung S_ads.
        //  nCouple = sub-iterasi q*<->T (Bab 10.1, nSolidCouple).
        void correct
        (
            const label nCouple,
            speciesImplicitCoupler* coupler = nullptr
        );

        //- Source energi adsorpsi [W/m3]
        const volScalarField& Sads() const noexcept { return Sads_; }

        //- Akses field/daftar
        const wordList& species() const noexcept { return species_; }
        const volScalarField& C(const label i) const { return C_[i]; }
        const volScalarField& q(const label i) const { return q_[i]; }

        //- Akses tambahan untuk diagnostik neraca massa interface (P1)
        const volScalarField& Ds(const label i) const { return Ds_[i]; }
        const volScalarField& effStor(const label i) const { return effStor_[i]; }
        scalar epsilonP() const noexcept { return epsilonP_; }
        scalar epsilonPtot() const noexcept { return epsilonPtot_; }
        scalar rhoS() const noexcept { return rhoS_; }
        scalar M(const label i) const noexcept { return M_[i]; }

        //- Adveksi species oleh fluks Darcy pori: persamaan C_i mendapat suku + div(phiD C_i)
        void setPoreFlux(const surfaceScalarField* phiD, const volScalarField* sink = nullptr) noexcept
        {
            poreFlux_ = phiD; poreSink_ = sink;
        }
        void setPoreSource(const PtrList<volScalarField>* src) noexcept { poreSrc_ = src; }
        bool poreAdvection() const noexcept { return poreFlux_ != nullptr; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
