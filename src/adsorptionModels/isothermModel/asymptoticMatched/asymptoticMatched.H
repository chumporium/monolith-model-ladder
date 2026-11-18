/*---------------------------------------------------------------------------*\
  asymptoticMatched.H

  PEMBUNGKUS isoterm: mencocokkan isoterm fisis dengan cabang LINEAR di
  konsentrasi rendah lewat teknik asymptotic matching.

      q*(p,T) = [ (q0)^n + (qp)^n ]^(1/n)        dgn  q0 = K c = K p/(R T)

  qp adalah isoterm fisis (sub-dict `base`), K kemiringan cabang linear
  [m3/kg], dan n eksponen pencocokan. Untuk n NEGATIF rumus ini adalah
  "minimum halus": ia mengikuti cabang yang lebih KECIL, dengan peralihan mulus
  di sekitar titik potong kedua cabang.

  KENAPA ADA
  ----------
  Isoterm cekung (Langmuir, Toth) membuat gelombang adsorpsi MENAJAM SENDIRI:
  konsentrasi tinggi bergerak lebih cepat daripada yang rendah, sehingga bagian
  belakang mengejar bagian depan sampai terbentuk front mirip-kejut. Skema
  berbasis fungsi basis kontinu (mis. elemen hingga Galerkin) tidak dapat
  merepresentasikan front seperti itu tanpa overshoot besar, kecuali jaringnya
  sangat halus.

  Antunes Morgado & Mazzotti, Comput. Chem. Eng. 205 (2026) 109404, memakai
  bentuk ini untuk meredam penajaman itu (SI Pers. 10-12). Model di sini
  dibuat agar solver ini dapat menjalankan model MEREKA apa adanya, sehingga
  verifikasi silang Tahap 1.2 membandingkan dua PENYELESAIAN dari model yang
  sama -- bukan dua model yang berbeda.

  /!\ INI BUKAN FISIKA. Ia peranti numerik yang MENGUBAH isoterm: di bawah
      ~25 ppm ia menurunkan loading sampai puluhan persen (79% pada 1 ppm untuk
      parameter CACE). Jangan memakainya untuk hasil yang diklaim fisis kecuali
      memang sedang mereproduksi karya yang memakainya.

  DICT
  ----
      isotherm
      {
          type    asymptoticMatched;
          K       262;        // [m3/kg] kemiringan cabang linear q0 = K c
          n       -6;         // [-] eksponen; WAJIB negatif (minimum halus)
          base                // isoterm fisis qp -- isotherm apa pun
          {
              type    Toth;
              T0      298;
              ns0     1.5482;
              b0      41.890;
              t0      0.1703;
              deltaH0 70000;
          }
      }

  Nilai n dapat diturunkan dari rasio pendekatan phi yang dikehendaki di titik
  potong:  n = 1/log2(phi)  (SI Pers. 11). phi = 0,9 memberi n = -6,58, yang
  dibulatkan CACE ke -6.
\*---------------------------------------------------------------------------*/
#ifndef isotherm_asymptoticMatched_H
#define isotherm_asymptoticMatched_H

#include "isothermModel.H"
#include "autoPtr.H"

namespace Foam
{
namespace isotherms
{

/*---------------------------------------------------------------------------*\
                    Class asymptoticMatched Declaration
\*---------------------------------------------------------------------------*/

class asymptoticMatched
:
    public isothermModel
{
    // Private data

        //- Isoterm fisis yang dibungkus
        autoPtr<isothermModel> base_;

        //- Kemiringan cabang linear q0 = K c   [m3/kg]
        scalar K_;

        //- Eksponen pencocokan (negatif)
        scalar n_;


    // Private Member Functions

        //- q0 = K c = K p/(R T) dan turunannya
        inline scalar q0(const scalar p, const scalar T) const;

public:

    //- Runtime type information
    TypeName("asymptoticMatched");

    // Constructors
        asymptoticMatched(const word& specieName, const dictionary& dict);

    //- Destructor
    virtual ~asymptoticMatched() = default;


    // Member Functions

        virtual scalar qStar(const scalar p, const scalar T) const;
        virtual scalar dqStardp(const scalar p, const scalar T) const;
        virtual scalar dqStardT(const scalar p, const scalar T) const;

        //- Panas adsorpsi diteruskan dari isoterm fisis: pencocokan asimtotik
        //  adalah peranti numerik dan tidak mengubah energetika adsorpsi.
        virtual scalar heatOfAdsorption() const
        {
            return base_->heatOfAdsorption();
        }

        //- Bukan bentuk Langmuir: cabang linear merusak bentuk itu, sehingga
        //  extendedLangmuir TIDAK boleh menganggapnya Langmuir.
        virtual bool isLangmuirForm() const { return false; }

        //- Akses isoterm fisis yang dibungkus (untuk diagnostik/uji)
        const isothermModel& base() const { return base_(); }
};

} // End namespace isotherms
} // End namespace Foam

#endif

// ************************************************************************* //
