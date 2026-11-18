/*---------------------------------------------------------------------------*\
  intraGrainDiffusion.H

  Closure DIFUSI INTRA-GRAIN (HSDM, homogeneous solid diffusion model).
  Resistansi sub-grid di lapisan amine DIRESOLUSI sebagai difusi 1-D fase
  teradsorpsi di dalam domain amine, dengan KESETIMBANGAN ISOTERM di permukaan
  domain (q_permukaan = q*). Loading rata-rata domain q-bar dikembalikan ke
  persamaan species makro sebagai source -rho_s M (q-bar - q-bar_old)/dt
  (jalur finite-rate, sama seperti LDF/dualSite).

  Geometri slab (default): film amine tebal delta yang melapisi dinding pori.
  Satu muka (sisi gas pori) setimbang dgn q*; muka lain (menempel rangka padat)
  no-flux. Adsorbat berdifusi menembus film dengan difusivitas padatan Dc.

      dq/dt = Dc d2q/dx2 ,  0<x<delta ;  q(delta)=q*(C,T) ;  dq/dx|_0 = 0

  Ekor breakthrough LAHIR dari profil difusi yang belum rata; satu-satunya
  parameter fisik baru = Dc (difusivitas padatan), dgn skala waktu tau=delta^2/Dc.
  TIDAK ada laju k atau eta yang difit.

  ⚠️ delta dan Dc TIDAK teridentifikasi terpisah. Keduanya hanya masuk lewat
  h=delta/nShells, sehingga r = Dc*dt*nShells^2/delta^2 = dt*nShells^2/tau:
  HANYA tau yang menyentuh fisika. delta berperan sbg pembawa satuan. Karena itu
  JANGAN kutip Dc sebagai difusivitas material (CATATAN 23.26a).

  Ekuivalen LDF momen-pertama (konvensi Glueckauf, sama spt k_s=15*Ds/r_c^2 pada
  paper acuan pers. 25/26): waktu serapan rata-rata slab = tau/3 TEPAT, jadi
      k_LDF = 3*Dc/delta^2 = 3/tau
  Jangan tertukar dgn pi^2/(4*tau) ~ 2.47/tau -- itu laju peluruhan mode
  terlambat, bukan momen pertama. Turunan lengkap: CATATAN 23.26(b).

  Diskretisasi: nShells kulit planar seragam; Euler implisit (tak-bersyarat
  stabil) -> sistem tridiagonal per sel. Koefisien identik tiap sel (Dc,dt,h
  konstan) -> faktorisasi Thomas SEKALI per langkah, back-substitusi per sel.

  ---------------------------------------------------------------------------
  KOREKSI DARKEN (opsional, `darkenCorrection yes`) -- CATATAN 23.42/23.43
  ---------------------------------------------------------------------------
  Difusivitas Fick jadi bergantung loading lewat faktor termodinamika isoterm:

      D(q) = D0 * Gamma(q) ,   Gamma = dln p/dln q = 1/(1-(q/ns)^t)  [Toth]

  D0 = difusivitas TERKOREKSI (Maxwell-Stefan), menggantikan Dc. JUMLAH
  PARAMETER TETAP SATU. Gamma datang gratis dari isoterm yang sudah ada.

  ⚠️ D0 BUKAN difusivitas-diri. D_self <= D0, berimpit hanya di limit encer,
  dan theta=0,41 sudah tercapai di 400 ppm -> kita tidak di sana. D0 paling
  jujur disebut PARAMETER EFEKTIF HASIL FIT, sama seperti Dc (CATATAN 23.43a).

  BENTUK KONSERVATIF WAJIB. Dgn D bervariasi persamaannya

      dq/dt = d/dx ( D0 Gamma(q) dq/dx )

  BUKAN D d2q/dx2. Memakai bentuk lama dgn D variabel TIDAK konservatif dan
  akan membocorkan massa -- persis kelas galat yang neraca dua sisi dibangun
  untuk menangkap (23.37/23.40). Diskretisasi di sini berbasis MUKA:
      F_{i+1/2} = D_{i+1/2} (q_{i+1}-q_i)/h
  sehingga fluks yang keluar sel i identik dgn yang masuk sel i+1 (teleskopik).

  ⚠️ D_{i+1/2} = rata-rata HARMONIK 2 D_i D_{i+1}/(D_i + D_{i+1}) -- DIUBAH
  2026-08-14 dari aritmetik (CATATAN 23.94). Kedua kulit adalah resistansi SERI;
  aritmetik didominasi nilai BESAR padahal fluks dibatasi nilai KECIL, jadi ia
  MELEBIHKAN konduktansi. Galat A/H = (1+r)^2/(4r), r = kontras antar kulit:
      Darken Gamma 1..5      r=1,26 -> +1,3%   (karena itu aritmetik memadai dulu)
      self-limiting a=10,9   r=1,89 -> +10,4%  @400 ppm
      self-limiting          r=3,23 -> +38,5%  @5,6%
  Konservasi TIDAK bergantung pilihan ini (nilai muka yang sama dipakai kedua
  sel); yang bergantung adalah KETELITIANNYA.

  Gamma di-LAG pada oldTime (linearisasi Picard) supaya sistem tetap linear
  tridiagonal. Konsisten dgn cara qStar sendiri dievaluasi pada langkah luar;
  nOuterCorrectors yang menutup sisa ketaklinearannya.

  BIAYA: koefisien tak lagi identik antar sel -> faktorisasi Thomas per sel
  (~18.000/langkah, bukan sekali). Closure naik ~2,5x; dampak total kecil krn
  solver tekanan yang dominan.

  Entri dict:
      type      intraGrainDiffusion;
      Dc        1e-15;   // [m2/s] difusivitas padatan (TANPA koreksi Darken)
      delta     1e-6;    // [m]    tebal film/jari-jari domain
      nShells   8;       // jumlah kulit (default 8)

      darkenCorrection yes;   // opsional, default no
      D0        2e-16;        // [m2/s] WAJIB bila darkenCorrection yes;
                              //        menggantikan Dc, BUKAN menskalakannya
  ⚠️ `D0` sengaja dibuat kunci TERPISAH dan WAJIB, bukan menafsirkan ulang
  `Dc`: keduanya berbeda ~Gamma (5x di 400 ppm), dan menafsir ulang angka
  hasil fit secara diam-diam adalah cara termudah menghasilkan run yang salah
  tanpa jejak.

  Valid untuk: kemisorpsi amine dgn difusi lapisan lambat (Bab 7.3/7.4 panduan).

  ---------------------------------------------------------------------------
  GEOMETRI BOLA (opsional, `geometry sphere;`) -- ditambahkan 2026-09-16
  ---------------------------------------------------------------------------
  Untuk difusi MIKROPORI di kristal zeolit (validasi Grande dkk. 2006, 4A):
      dq/dt = (Dc/r^2) d/dr ( r^2 dq/dr ) ,  q(R)=q*(C,T) ,  dq/dr|_0 = 0
  `delta` dibaca sebagai JARI-JARI kristal R. Default tetap `slab`, dan jalur
  slab TIDAK disentuh -> kasus lama bit-identik.

  Volume hingga konservatif pada kulit radial seragam h = R/N:
      bobot kulit  w_i   = ((i+1)^3 - i^3)/3      (dalam h^3)
      luas muka    A_f   = f^2                    (dalam h^2), f = i atau i+1
      rL_i = r i^2/w_i ,  rR_i = r (i+1)^2/w_i ,  r = Dc dt/h^2
      muka permukaan berjarak h/2 -> rR_{N-1} = 2 r N^2/w_{N-1}
  Fluks muka f ditulis identik dari kedua sisi (w_i rR_i = w_{i+1} rL_{i+1}
  = r f^2) -> teleskopik. q-bar = rata-rata BERBOBOT VOLUME, bukan aritmetik.

  Akurasi (momen pertama uptake terhadap analitik R^2/(15 Dc), Euler implisit
  memberi momen pertama semi-diskret TEPAT berapa pun dt):
      N=8 +2,6% | N=16 +0,65% | N=32 +0,16%   (orde 2)
  -> pakai nShells >= 16 untuk bola.
  Ekuivalen LDF Glueckauf: k = 15 Dc/R^2.

  ⛔ Bola + Darken BELUM diturunkan dan ditolak dgn FatalError.

  Entri dict tambahan:
      geometry  sphere;   // slab (default) | sphere

  ---------------------------------------------------------------------------
  BOLA + SELF-LIMITING + SITUS PERMUKAAN -- ditambahkan 2026-09-26
  ---------------------------------------------------------------------------
  Struktur model Grossmann & Mazzotti, IECR 64 (2025) 2339, pers. 8-15/37:
      q      = q1 + q2
      q1     = psi q*(C,T)                         situs permukaan, setimbang seketika
      dq2/dt = (1/r^2) d/dr ( r^2 D(q2) dq2/dr ),  q2(R) = (1-psi) q* ,  dq2/dr|_0 = 0
      D(q2)  = D0 exp(-alpha (q2/qSatSL)^beta)     self-limiting (cross-linking amina)
  Kulit menyimpan profil q2 (bulk). Perakitan muka bola konservatif dgn D muka
  HARMONIK (kernel statik assembleSphereVar); pada D seragam tereduksi TEPAT ke
  assembleSphere. Faktor implisit = psi + (1-psi) g_bulk.
      surfaceFraction  0.26;   // psi, 0 <= psi < 1 ; default 0 (hanya dgn sphere+selfLimiting)
      qSatSL           0.038;  // [mol/kg, basis q TERSIMPAN solver] pasangan alpha dari
                               // sumber yg SAMA; bila tak diberi: theta = q_iso/qSat(T) (lama)
  Tanpa kedua kunci baru: jalur lama bit-identik.
\*---------------------------------------------------------------------------*/

#ifndef closure_intraGrainDiffusion_H
#define closure_intraGrainDiffusion_H

#include "kineticClosure.H"
#include "PtrList.H"

namespace Foam
{
namespace kineticClosures
{

class intraGrainDiffusion
:
    public kineticClosure
{
    // Parameter fisik. Dc_ = difusivitas basis: Dc (tanpa Darken) atau D0
    // (dengan Darken). Hanya SATU yang pernah difit.
    scalar Dc_;
    scalar delta_;
    label  nShells_;

    // Geometri domain: false = slab (default, jalur lama), true = bola
    bool sphere_;

    // --- Koreksi Darken (opsional) ---
    bool darken_;

    // --- Difusi SELF-LIMITING (opsional, EKSKLUSIF thd Darken) ---
    //   D(theta) = D0 exp(-alpha theta^beta) ,  theta = q_isoterm/qSat(T)
    // Mekanisme: reaksi CO2-amina menciptakan cross-linking yang MENGHAMBAT
    // difusi lebih lanjut (Vallace dkk. CEJ 2023; Grossmann & Mazzotti I&ECR
    // 2025, 64, 2339). Arahnya BERLAWANAN Darken -- lihat CATATAN 23.90/23.95.
    bool   selfLim_;
    scalar alphaSL_;   // [-] sensitivitas cakupan, tak berdimensi -> DAPAT dipindah
    scalar betaSL_;    // [-] default 1 (Tabel 4 paper: beta ~ 1, tidak difit)

    // --- Situs permukaan & basis theta eksplisit (2026-09-26) ---
    scalar psi_;       // fraksi situs permukaan (setimbang seketika); 0 = tanpa
    scalar qSatSL_;    // [mol/kg, basis q tersimpan] theta = q/qSatSL; <=0 = basis lama

    // Konteks dari regionAdsorption lewat setIsothermContext(). Pointer,
    // bukan referensi, krn diisi SETELAH konstruksi.
    const isothermModel* iso_;
    const volScalarField* T_;
    scalar fActive_;

    // Gamma rata-rata volume langkah terakhir, dipakai dqImplicitFactor()
    // yang harus mengembalikan SATU skalar. Direduksi lintas-rank supaya
    // linearisasi tak bergantung dekomposisi. 1 sebelum update() pertama.
    scalar gammaMean_;

    // State sub-grid: profil loading per kulit [mol/kg], kulit 0=backing
    // (no-flux), kulit nShells-1 = permukaan (setimbang q*).
    PtrList<volScalarField> shells_;

    // Selesaikan slab 1-D implisit untuk satu kolom kulit.
    //   old[]  : loading kulit pd oldTime (sumber)
    //   qs     : q* permukaan
    //   cprime,invm : koefisien Thomas ter-faktor (konstan, dihitung sekali)
    //   r      : Dc dt / h^2
    //   xout[] : loading kulit baru ; return q-bar (rata-rata kulit)
    scalar solveColumn
    (
        const scalarList& old,
        const scalar qs,
        const scalarList& cprime,
        const scalarList& invm,
        const scalar r,
        scalarList& dwork,
        scalarList& xout
    ) const;

    // Versi D BERVARIASI (koreksi Darken). Merakit tridiagonal berbasis muka
    // lalu memfaktorkan Thomas UNTUK KOLOM INI SAJA.
    //   Tcell  : suhu sel [K] (Gamma bergantung T lewat ns(T), t(T))
    //   gSum   : akumulator Gamma (ditambah rata-rata kulit kolom ini)
    // Mengembalikan q-bar.
    //- Faktor keadaan pengali Dc_: Gamma(q) untuk Darken, exp(-alpha th^b)
    //  untuk self-limiting. SATU tempat, supaya kedua hukum memakai jalur
    //  perakitan konservatif yang sama.
    scalar stateFactor(const scalar qStored, const scalar Tcell) const;

    scalar solveColumnDarken
    (
        const scalarList& old,
        const scalar qs,
        const scalar Tcell,
        const scalar dt,
        scalarList& Dsh,
        scalarList& adiag,
        scalarList& bdiag,
        scalarList& cdiag,
        scalarList& dwork,
        scalarList& cprime,
        scalarList& xout,
        scalar& gSum
    ) const;

    //- Bola + D(q) (self-limiting): satu kolom. qsBulk = (1-psi) q*.
    //  Mengembalikan q2-bar (rata-rata berbobot volume).
    scalar solveColumnSphereVar
    (
        const scalarList& old,
        const scalar qsBulk,
        const scalar Tcell,
        const scalar dt,
        scalarList& Dsh,
        scalarList& a,
        scalarList& b,
        scalarList& c,
        scalarList& rhs,
        scalarList& w,
        scalarList& cp,
        scalarList& dw,
        scalarList& x,
        scalar& gSum
    ) const;

    //- Jalur D konstan geometri bola (dipanggil update() bila sphere_)
    void updateSphere
    (
        volScalarField& q,
        const volScalarField& qStar,
        const scalar dt
    );

public:

    TypeName("intraGrainDiffusion");

    // --- Kernel numerik bola: statik, tanpa mesh (dipakai juga uji satuan) ---

        //- Rakit tridiagonal kulit radial untuk D konstan.
        //  r = Dc dt/h^2. Keluaran a (sub), b (diag), c (super), sRhs
        //  (koefisien q* di baris permukaan), w (bobot volume kulit).
        static void assembleSphere
        (
            const label N,
            const scalar r,
            scalarList& a,
            scalarList& b,
            scalarList& c,
            scalar& sRhs,
            scalarList& w
        );

        //- Rakit tridiagonal kulit radial untuk D per kulit (Dsh, N nilai) dan
        //  D permukaan Dsurf; D muka = rata-rata harmonik. dth2 = dt/h^2.
        //  Pada Dsh = Dsurf = Dc seragam identik assembleSphere(r = Dc dth2).
        static void assembleSphereVar
        (
            const label N,
            const scalar dth2,
            const scalarList& Dsh,
            const scalar Dsurf,
            scalarList& a,
            scalarList& b,
            scalarList& c,
            scalar& sRhs,
            scalarList& w
        );

        //- Thomas umum lalu rata-rata berbobot w. clamp -> x >= 0.
        static scalar solveTridiagMean
        (
            const scalarList& a,
            const scalarList& b,
            const scalarList& c,
            const scalarList& rhs,
            const scalarList& w,
            const bool clamp,
            scalarList& cp,
            scalarList& dw,
            scalarList& x
        );

    intraGrainDiffusion
    (
        const fvMesh& mesh,
        const word& specieName,
        const dictionary& dict
    );

    virtual bool isAlgebraic() const { return false; }

    virtual void update
    (
        volScalarField& q,
        const volScalarField& qStar,
        const scalar dt
    );

    virtual void storeOldTime();

    virtual void setIsothermContext
    (
        const isothermModel& iso,
        const volScalarField& T,
        const scalar fActive
    );

    //- g = d(q-bar)/d(q*) = respons q-bar thd q* permukaan utk SATU langkah
    //  difusi implisit dari kondisi awal nol (operator linear). Skalar (r sama
    //  tiap sel). Utk suku implisit di regionAdsorption jalur finite-rate.
    virtual scalar dqImplicitFactor(const scalar dt) const;
};

} // End namespace kineticClosures
} // End namespace Foam

#endif

// ************************************************************************* //
