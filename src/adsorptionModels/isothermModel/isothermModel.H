/*---------------------------------------------------------------------------*\
  isothermModel.H

  Base class pluggable untuk isoterm KOMPONEN MURNI (Bab 6 panduan).

  Variabel natural isoterm adalah TEKANAN PARSIAL p [Pa] dan suhu T [K];
  loading dikembalikan dalam [mol/kg]. Konversi dari/ke fraksi massa dilakukan
  di luar (mesin konversi, Bab 10.5) sehingga isoterm ditulis persis seperti
  publikasi.

  Antarmuka wajib (per komponen murni):
      qStar(p,T)      loading kesetimbangan          [mol/kg]
      dqStardp(p,T)   turunan qStar terhadap p (T tetap)   [mol/kg/Pa]
      dqStardT(p,T)   turunan qStar terhadap T (p tetap)   [mol/kg/K]

  Antarmuka opsional (untuk aturan campuran):
      isLangmuirForm()             true jika q* = qSat*b*p/(1+b*p)
      qSat(T), b(T)                kapasitas & afinitas (extended Langmuir)
      reducedSpreadingPressure()   psi(p,T) = INT_0^p q(p')/p' dp'   (IAST)
      pressureFromPsi()            inversi psi -> p                  (IAST)

  reducedSpreadingPressure & pressureFromPsi memiliki implementasi NUMERIK
  default di base (kuadratur + bisection) sehingga IAST bekerja untuk SEMUA
  isoterm murni (termasuk Toth), sesuai panduan Bab 6.3. Subkelas boleh
  meng-override dengan bentuk analitik.
\*---------------------------------------------------------------------------*/

#ifndef isothermModel_H
#define isothermModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"
#include "adsorptionConversion.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class isothermModel Declaration
\*---------------------------------------------------------------------------*/

class isothermModel
{
protected:

    //- Nama species pemilik isoterm
    const word specieName_;

    //- Salinan sub-dictionary koefisien
    const dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("isothermModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        isothermModel,
        dictionary,
        (
            const word& specieName,
            const dictionary& dict
        ),
        (specieName, dict)
    );


    // Constructors

        isothermModel(const word& specieName, const dictionary& dict);

        //- Disallow default bitwise copy construction
        isothermModel(const isothermModel&) = delete;


    // Selector
    static autoPtr<isothermModel> New
    (
        const word& specieName,
        const dictionary& dict
    );


    //- Destructor
    virtual ~isothermModel() = default;


    // Member Functions

        const word& specieName() const noexcept { return specieName_; }

        // --- Antarmuka wajib ---

            //- Loading kesetimbangan q*(p,T) [mol/kg]
            virtual scalar qStar(const scalar p, const scalar T) const = 0;

            //- d q*/d p |_T [mol/kg/Pa]
            virtual scalar dqStardp(const scalar p, const scalar T) const = 0;

            //- d q*/d T |_p [mol/kg/K]
            virtual scalar dqStardT(const scalar p, const scalar T) const = 0;


        // --- Antarmuka opsional (campuran) ---

            //- Apakah berbentuk Langmuir murni (untuk extended Langmuir)
            virtual bool isLangmuirForm() const { return false; }

            //- Kapasitas jenuh [mol/kg] (Langmuir-form). Default: q*->inf approx
            virtual scalar qSat(const scalar T) const
            {
                return qStar(GREAT, T);
            }

            //- Afinitas b(T) [1/Pa] (Langmuir-form). Default dari kemiringan awal
            virtual scalar b(const scalar T) const
            {
                const scalar qs = qSat(T);
                return (qs > SMALL) ? dqStardp(scalar(0), T)/qs : scalar(0);
            }

        // --- Panas adsorpsi: SUMBER TUNGGAL ---
        //
        //  Panas adsorpsi isosterik [J/mol], positif = eksotermik. Ia besaran
        //  fisis yang SAMA dengan `deltaH` pada source energi solid, tetapi
        //  sebelumnya keduanya hidup di dua entri dict terpisah tanpa
        //  pemeriksaan silang. Kelas jebakan yang sama dgn rho_s dan eps_p.
        //
        //  Isoterm yang memakai satu nilai untuk b(T) mengembalikannya di sini,
        //  sehingga regionAdsorption dapat memakainya sebagai default dan
        //  memeriksa silang `deltaH` yang ditulis tangan.
        //
        //  Nilai 0 berarti "isoterm tidak menyatakan satu nilai" -- dipakai
        //  DualSiteLangmuir, yang punya DUA energi situs dan tidak dapat
        //  meringkasnya jadi satu tanpa berbohong. Untuk isoterm semacam itu
        //  `deltaH` tetap WAJIB ditulis di dict.

            //- Panas adsorpsi isosterik [J/mol]; 0 = tidak menyatakan
            virtual scalar heatOfAdsorption() const { return scalar(0); }

        // --- Faktor termodinamika (koreksi Darken) ---
        //
        //  Gamma == dln p/dln q  pada T tetap. Dipakai kineticClosure yang
        //  memodelkan difusi fase teradsorpsi dgn difusivitas TERKOREKSI D0:
        //      D_Fick(q) = D0 * Gamma(q)
        //  Gamma >= 1 selalu; Gamma -> 1 di limit encer.
        //
        //  ⚠️ Argumen q dalam [mol/kg adsorben AKTIF] — skala isoterm sendiri,
        //  BUKAN mol/kg material. Pemanggil wajib membagi fActive lebih dulu
        //  (regionAdsorption menyimpan q sbg fActive*q_isoterm).
        //
        //  Sengaja dirumuskan dalam q, bukan p: loading per kulit sudah ada di
        //  closure, sedangkan p menuntut inversi isoterm per kulit. Untuk Toth
        //  kedua bentuk analitik dan terverifikasi identik s/d 1e-14
        //  (CATATAN 23.42c).

            //- Apakah subkelas menyediakan Gamma analitik
            virtual bool hasThermodynamicFactor() const { return false; }

            //- Gamma(q,T) [-]; default 1 (= tanpa koreksi)
            virtual scalar thermodynamicFactor
            (
                const scalar q,
                const scalar T
            ) const
            {
                return scalar(1);
            }


            //- Reduced spreading pressure psi(p,T) = INT_0^p q/p' dp' [mol/kg]
            //  Default: kuadratur Gauss-Legendre adaptif (substitusi log).
            virtual scalar reducedSpreadingPressure
            (
                const scalar p,
                const scalar T
            ) const;

            //- Inversi: cari p sehingga reducedSpreadingPressure(p,T) = psi.
            //  Default: bisection.
            virtual scalar pressureFromPsi
            (
                const scalar psi,
                const scalar T
            ) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //
