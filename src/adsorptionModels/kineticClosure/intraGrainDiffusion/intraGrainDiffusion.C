/*---------------------------------------------------------------------------*\
  intraGrainDiffusion.C
\*---------------------------------------------------------------------------*/

#include "intraGrainDiffusion.H"
#include "isothermModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"

namespace Foam
{
namespace kineticClosures
{
    defineTypeNameAndDebug(intraGrainDiffusion, 0);
    addToRunTimeSelectionTable(kineticClosure, intraGrainDiffusion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

Foam::kineticClosures::intraGrainDiffusion::intraGrainDiffusion
(
    const fvMesh& mesh,
    const word& specieName,
    const dictionary& dict
)
:
    kineticClosure(mesh, specieName, dict),
    // Sakelar dibaca dua kali: urutan inisialisasi anggota mengikuti urutan
    // DEKLARASI, jadi darken_ belum ada nilainya saat Dc_ diinisialisasi.
    // Membaca ulang dari dict lebih jelas daripada menyusun ulang anggota.
    Dc_
    (
        (   dict.getOrDefault<bool>("darkenCorrection", false)
         || dict.getOrDefault<bool>("selfLimiting", false) )
      ? dict.get<scalar>("D0")
      : dict.get<scalar>("Dc")
    ),
    delta_(dict.get<scalar>("delta")),
    nShells_(dict.getOrDefault<label>("nShells", 8)),
    sphere_(dict.getOrDefault<word>("geometry", "slab") == "sphere"),
    darken_(dict.getOrDefault<bool>("darkenCorrection", false)),
    selfLim_(dict.getOrDefault<bool>("selfLimiting", false)),
    alphaSL_
    (
        dict.getOrDefault<bool>("selfLimiting", false)
      ? dict.get<scalar>("alphaSL")
      : 0.0
    ),
    betaSL_(dict.getOrDefault<scalar>("betaSL", 1.0)),
    psi_(dict.getOrDefault<scalar>("surfaceFraction", 0.0)),
    qSatSL_(dict.getOrDefault<scalar>("qSatSL", -1.0)),
    iso_(nullptr),
    T_(nullptr),
    fActive_(1.0),
    gammaMean_(1.0),
    shells_(nShells_)
{
    if (nShells_ < 2)
    {
        FatalErrorInFunction
            << "nShells harus >= 2 (diberikan " << nShells_ << ")"
            << exit(FatalError);
    }

    // Darken dan self-limiting adalah mekanisme BERLAWANAN (D naik lawan D
    // turun terhadap loading). Menyalakan keduanya bukan "kombinasi", ia
    // kekeliruan konfigurasi yang hasilnya akan tampak sah.
    if (darken_ && selfLim_)
    {
        FatalErrorInFunction
            << "darkenCorrection DAN selfLimiting sama-sama aktif." << nl
            << "Keduanya mekanisme BERLAWANAN: Darken menaikkan D dgn loading,"
            << " self-limiting menurunkannya (CATATAN 23.90/23.95)." << nl
            << "Pilih SATU." << exit(FatalError);
    }

    {
        const word geom = dict.getOrDefault<word>("geometry", "slab");
        if (geom != "slab" && geom != "sphere")
        {
            FatalErrorInFunction
                << "geometry '" << geom << "' tidak dikenal. Pilihan: slab"
                << " (default) | sphere." << exit(FatalError);
        }
    }

    // Bola + Darken belum diturunkan. Bola + self-limiting: assembleSphereVar.
    if (sphere_ && darken_)
    {
        FatalErrorInFunction
            << "geometry sphere belum mendukung darkenCorrection." << nl
            << "Perakitan konservatif bola D(q) hanya diuji untuk selfLimiting."
            << exit(FatalError);
    }

    // Situs permukaan hanya didefinisikan pada struktur model Grossmann &
    // Mazzotti (bola + self-limiting). Menyalakannya di jalur lain akan
    // diam-diam diabaikan -> ditolak.
    if (psi_ < 0 || psi_ >= 1)
    {
        FatalErrorInFunction
            << "surfaceFraction = " << psi_ << " di luar [0, 1)." << exit(FatalError);
    }
    if (psi_ > 0 && !(sphere_ && selfLim_))
    {
        FatalErrorInFunction
            << "surfaceFraction > 0 hanya didukung bersama geometry sphere +"
            << " selfLimiting." << exit(FatalError);
    }
    if (qSatSL_ > 0 && !selfLim_)
    {
        FatalErrorInFunction
            << "qSatSL diberikan tetapi selfLimiting tidak aktif." << exit(FatalError);
    }

    if (selfLim_ && alphaSL_ <= 0)
    {
        FatalErrorInFunction
            << "selfLimiting aktif tetapi alphaSL = " << alphaSL_
            << " (harus > 0)." << nl
            << "alphaSL adalah sensitivitas cakupan TAK BERDIMENSI; ia yang"
            << " dapat dipindahkan antar material, bukan D0." << nl
            << exit(FatalError);
    }

    const dimensionSet dimLoad(dimMoles/dimMass);

    forAll(shells_, k)
    {
        shells_.set
        (
            k,
            new volScalarField
            (
                IOobject
                (
                    "qShell" + Foam::name(k) + "." + specieName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("qShell", dimLoad, 0)
            )
        );
        shells_[k].storeOldTime();
    }

    Info<< "    intraGrainDiffusion[" << specieName << "]: "
        << (darken_ ? "D0=" : "Dc=") << Dc_
        << " m2/s, delta=" << delta_ << " m, nShells=" << nShells_
        << " ; tau=delta^2/" << (darken_ ? "D0" : "Dc") << "="
        << delta_*delta_/max(Dc_, VSMALL) << " s"
        << " (k_LDF Glueckauf = 3D/delta^2 = "
        << 3.0*Dc_/max(delta_*delta_, VSMALL) << " 1/s)" << endl;

    if (sphere_)
    {
        Info<< "        geometri BOLA: delta = jari-jari R; k_LDF Glueckauf"
            << " yang setara = 15 Dc/R^2 = "
            << 15.0*Dc_/max(delta_*delta_, VSMALL) << " 1/s"
            << " (angka k_LDF slab di atas TIDAK berlaku)." << endl;

        if (nShells_ < 16)
        {
            Info<< "        ⚠️ nShells=" << nShells_ << " < 16: momen pertama"
                << " uptake bola meleset +2,6% pada N=8 (+0,65% pada N=16)."
                << endl;
        }
    }

    if (darken_)
    {
        Info<< "        koreksi DARKEN AKTIF: D(q) = D0*Gamma(q), Gamma dari"
            << " isoterm; tau di atas adalah tau0 pd Gamma=1 (limit encer),"
            << " jadi tau EFEKTIF lebih kecil sebesar Gamma." << endl;
    }

    if (selfLim_)
    {
        Info<< "        difusi SELF-LIMITING AKTIF: D(theta) = D0*exp(-"
            << alphaSL_ << "*theta^" << betaSL_ << "), theta = q/qSat." << nl
            << "        tau di atas = tau(theta=0); tau EFEKTIF LEBIH BESAR"
            << " sebesar exp(+alpha theta^beta). Resistansi pengendali ada di"
            << " KULIT PERMUKAAN (23.95)." << endl;

        if (qSatSL_ > 0)
        {
            Info<< "        basis theta EKSPLISIT: theta = q_tersimpan/qSatSL, qSatSL = "
                << qSatSL_ << " mol/kg (basis q solver; pasangan alpha dari sumber"
                << " yang sama, TIDAK dikonversi)." << endl;
        }
        if (sphere_)
        {
            Info<< "        bola + self-limiting: situs permukaan psi = " << psi_
                << " (setimbang seketika), bulk (1-psi) berdifusi; q = psi q* +"
                << " q2-bar." << endl;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void Foam::kineticClosures::intraGrainDiffusion::storeOldTime()
{
    forAll(shells_, k)
    {
        shells_[k].storeOldTime();
    }
}


Foam::scalar Foam::kineticClosures::intraGrainDiffusion::solveColumn
(
    const scalarList& old,
    const scalar qs,
    const scalarList& cprime,
    const scalarList& invm,
    const scalar r,
    scalarList& dwork,
    scalarList& xout
) const
{
    const label N = nShells_;

    // RHS: d_i = old_i ; tambah suku permukaan di kulit terluar (Dirichlet q*)
    for (label i = 0; i < N; ++i) { dwork[i] = old[i]; }
    dwork[N - 1] += 2.0*r*qs;

    // Forward sweep (a_i = -r untuk i>=1)
    dwork[0] *= invm[0];
    for (label i = 1; i < N; ++i)
    {
        dwork[i] = (dwork[i] + r*dwork[i - 1])*invm[i];
    }

    // Back substitution + rata-rata (grid seragam -> mean aritmetik)
    scalar sum = 0;
    xout[N - 1] = max(dwork[N - 1], scalar(0));
    sum += xout[N - 1];
    for (label i = N - 2; i >= 0; --i)
    {
        xout[i] = max(dwork[i] - cprime[i]*xout[i + 1], scalar(0));
        sum += xout[i];
    }

    return sum/scalar(N);
}


void Foam::kineticClosures::intraGrainDiffusion::setIsothermContext
(
    const isothermModel& iso,
    const volScalarField& T,
    const scalar fActive
)
{
    iso_ = &iso;
    T_ = &T;
    fActive_ = max(fActive, SMALL);

    if (darken_ && !iso.hasThermodynamicFactor())
    {
        FatalErrorInFunction
            << "darkenCorrection AKTIF tetapi isoterm '" << iso.type()
            << "' tidak menyediakan faktor termodinamika Gamma." << nl
            << "Diam-diam memakai Gamma=1 akan menghasilkan run yang tampak"
            << " sah padahal identik dgn tanpa koreksi. Pakai isoterm Toth,"
            << " atau matikan darkenCorrection." << exit(FatalError);
    }

    if (selfLim_)
    {
        const scalar qs = iso.qSat(298.15);
        if (!(qs > SMALL) || qs > GREAT/2)
        {
            FatalErrorInFunction
                << "selfLimiting AKTIF tetapi isoterm '" << iso.type()
                << "' tidak memberi qSat yang berhingga (qSat(298 K) = "
                << qs << ")." << nl
                << "theta = q/qSat tak terdefinisi -> alphaSL tak bermakna."
                << exit(FatalError);
        }
    }
}


Foam::scalar Foam::kineticClosures::intraGrainDiffusion::stateFactor
(
    const scalar qStored,
    const scalar Tcell
) const
{
    const scalar qIso = qStored/fActive_;      // mol/kg AKTIF (basis isoterm)

    if (!selfLim_)
    {
        return iso_->thermodynamicFactor(qIso, Tcell);
    }

    // Basis eksplisit (qSatSL): theta = q_tersimpan/qSatSL, tanpa konversi.
    if (qSatSL_ > 0)
    {
        const scalar thE = qStored/qSatSL_;
        if (thE > 1.0 + 1e-6)
        {
            FatalErrorInFunction
                << "theta = q/qSatSL = " << thE << " MELEBIHI 1 (q = " << qStored
                << ", qSatSL = " << qSatSL_ << ")." << nl
                << "alpha dan qSatSL WAJIB berpasangan dari sumber yang sama."
                << exit(FatalError);
        }
        const scalar thbE =
            (betaSL_ == 1.0) ? thE : Foam::pow(max(thE, 0.0), betaSL_);
        return Foam::exp(-alphaSL_*thbE);
    }

    // theta = q_isoterm / qSat(T). Basisnya WAJIB isoterm yang sedang dipakai
    // -- alpha dan q_s hanya bermakna berpasangan (CATATAN 23.91[2], 23.92[2]).
    const scalar qSat = iso_->qSat(Tcell);

    if (qSat <= SMALL)
    {
        FatalErrorInFunction
            << "qSat = " << qSat << " <= 0 pada T = " << Tcell
            << " K; theta tak terdefinisi untuk selfLimiting." << exit(FatalError);
    }

    const scalar th = qIso/qSat;

    // ⛔ theta > 1 = loading MELEBIHI kapasitas jenuh isotermnya sendiri. Itu
    // tidak dapat terjadi bila alpha dan q_s berasal dari basis fit yang SAMA;
    // kalau muncul, ia menandakan pencampuran basis -- persis galat yang
    // CATATAN 23.92[1] buat dan 23.93 koreksi. Dikodekan supaya tak terulang
    // diam-diam.
    if (th > 1.0 + 1e-6)
    {
        FatalErrorInFunction
            << "theta = q/qSat = " << th << " MELEBIHI 1 (q_iso = " << qIso
            << " mol/kg aktif, qSat = " << qSat << ")." << nl
            << "exp(-alpha theta) kehilangan arti di sini." << nl
            << "Penyebab yang sudah pernah terjadi (23.92[1]): alpha diambil"
            << " dari satu basis fit dan q dari basis lain. alpha dan q_s WAJIB"
            << " berpasangan dari sumber yang sama." << nl
            << exit(FatalError);
    }

    const scalar thb = (betaSL_ == 1.0) ? th : Foam::pow(max(th, 0.0), betaSL_);
    return Foam::exp(-alphaSL_*thb);
}


Foam::scalar Foam::kineticClosures::intraGrainDiffusion::solveColumnDarken
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
) const
{
    const label N = nShells_;
    const scalar h = delta_/scalar(N);
    const scalar dth2 = dt/(h*h);

    // --- D per kulit, faktor keadaan di-LAG pd oldTime (Picard) ---
    //   Darken        : D = D0 Gamma(q)          -- Gamma NAIK dgn loading
    //   self-limiting : D = D0 exp(-alpha th^b)  -- D TURUN dgn loading
    // Keduanya lewat jalur yang sama; hanya hukum D-nya berbeda.
    scalar gAcc = 0;
    for (label i = 0; i < N; ++i)
    {
        const scalar f = stateFactor(old[i], Tcell);
        Dsh[i] = Dc_*f;
        gAcc += f;
    }
    gSum += gAcc/scalar(N);

    // D di permukaan: loading di sana = q* (Dirichlet).
    const scalar Dsurf = Dc_*stateFactor(qs, Tcell);

    // --- Koefisien BERBASIS MUKA (konservatif) ---
    //   rL_i = D_{i-1/2} dt/h^2 , rR_i = D_{i+1/2} dt/h^2
    //   D_{i+1/2} = rata-rata HARMONIK (lihat catatan di bawah)
    //   Muka permukaan berjarak h/2 dari pusat sel terakhir -> faktor 2.
    // Pada D konstan ini tereduksi TEPAT ke b = 1+r | 1+2r | 1+3r dgn
    // RHS += 2r q*, yaitu skema lama (harmonik = aritmetik saat D_i = D_{i+1}).
    //
    // ⚠️ KENAPA HARMONIK (CATATAN 23.94). Kedua kulit adalah resistansi SERI,
    // jadi konduktansi muka yang benar adalah harmonik. Rata-rata ARITMETIK
    // didominasi nilai BESAR padahal fluks dibatasi nilai KECIL -> ia
    // MELEBIHKAN konduktansi, yaitu laju serapan terlalu cepat.
    // Galatnya A/H = (1+r)^2/(4r) dgn r = kontras D antar kulit bertetangga:
    //   Darken (Gamma 1..5 @400 ppm)   r=1,26 -> +1,3%   (dapat diabaikan)
    //   Darken (Gamma 1..14,7 @5,6%)   r=1,47 -> +3,7%
    //   self-limiting (alpha=10,9)     r=1,89 -> ** +10,4% ** @400 ppm
    //   self-limiting                  r=3,23 -> ** +38,5% ** @5,6%
    // Aritmetik memadai untuk Darken; TIDAK memadai untuk self-limiting, dan
    // salahnya paling parah persis di 5,6% yang jadi sumbu kedua Paper 1.
    // Konservasi TIDAK bergantung cara D_{i+1/2} dihitung, asalkan nilai yang
    // SAMA dipakai kedua sel yang berbagi muka -- dan itu tetap dijaga di sini.
    for (label i = 0; i < N; ++i)
    {
        // harmonik: 2 D_a D_b / (D_a + D_b), aman thd D -> 0
        auto hmean = [](const scalar Da, const scalar Db)
        {
            return 2.0*Da*Db/max(Da + Db, VSMALL);
        };

        const scalar rL = (i == 0) ? 0.0 : hmean(Dsh[i - 1], Dsh[i])*dth2;
        const scalar rR = (i == N - 1)
                        ? 2.0*hmean(Dsh[N - 1], Dsurf)*dth2
                        : hmean(Dsh[i], Dsh[i + 1])*dth2;

        adiag[i] = -rL;
        bdiag[i] = 1.0 + rL + rR;
        cdiag[i] = (i == N - 1) ? 0.0 : -rR;
        dwork[i] = old[i];

        if (i == N - 1)
        {
            dwork[i] += rR*qs;
        }
    }

    // --- Thomas UNTUK KOLOM INI (koefisien tak lagi seragam) ---
    scalar m = bdiag[0];
    cprime[0] = cdiag[0]/m;
    dwork[0] = dwork[0]/m;
    for (label i = 1; i < N; ++i)
    {
        m = bdiag[i] - adiag[i]*cprime[i - 1];
        cprime[i] = (i < N - 1) ? cdiag[i]/m : 0.0;
        dwork[i] = (dwork[i] - adiag[i]*dwork[i - 1])/m;
    }

    scalar sum = 0;
    xout[N - 1] = max(dwork[N - 1], scalar(0));
    sum += xout[N - 1];
    for (label i = N - 2; i >= 0; --i)
    {
        xout[i] = max(dwork[i] - cprime[i]*xout[i + 1], scalar(0));
        sum += xout[i];
    }

    return sum/scalar(N);
}


void Foam::kineticClosures::intraGrainDiffusion::update
(
    volScalarField& q,
    const volScalarField& qStar,
    const scalar dt
)
{
    if (darken_ || selfLim_)
    {
        if (!iso_ || !T_)
        {
            FatalErrorInFunction
                << "darkenCorrection aktif tetapi setIsothermContext() belum"
                << " dipanggil. regionAdsorption seharusnya memanggilnya tepat"
                << " setelah closure dibuat." << exit(FatalError);
        }

        const label N = nShells_;
        scalarList old(N), Dsh(N), adiag(N), bdiag(N), cdiag(N),
                   dwork(N), cprime(N), xout(N);

        scalar gSum = 0;
        label nCol = 0;

        if (sphere_)
        {
            // Bola + self-limiting (+ situs permukaan psi). Kulit = profil q2.
            scalarList w(N), dw2(N);
            {
                scalarField& qc = q.primitiveFieldRef();
                const scalarField& qsF = qStar.primitiveField();
                const scalarField& Tf = T_->primitiveField();
                forAll(qc, c)
                {
                    for (label k = 0; k < N; ++k)
                    {
                        old[k] = shells_[k].oldTime().primitiveField()[c];
                    }
                    const scalar q2 = solveColumnSphereVar
                    (
                        old, (1.0 - psi_)*qsF[c], Tf[c], dt, Dsh,
                        adiag, bdiag, cdiag, dwork, w, cprime, dw2, xout, gSum
                    );
                    qc[c] = psi_*qsF[c] + q2;
                    ++nCol;
                    for (label k = 0; k < N; ++k)
                    {
                        shells_[k].primitiveFieldRef()[c] = xout[k];
                    }
                }
            }
            forAll(q.boundaryField(), patchi)
            {
                scalarField& qp = q.boundaryFieldRef()[patchi];
                const scalarField& qsF = qStar.boundaryField()[patchi];
                const scalarField& Tf = T_->boundaryField()[patchi];
                forAll(qp, f)
                {
                    for (label k = 0; k < N; ++k)
                    {
                        old[k] = shells_[k].oldTime().boundaryField()[patchi][f];
                    }
                    scalar dummy = 0;
                    const scalar q2 = solveColumnSphereVar
                    (
                        old, (1.0 - psi_)*qsF[f], Tf[f], dt, Dsh,
                        adiag, bdiag, cdiag, dwork, w, cprime, dw2, xout, dummy
                    );
                    qp[f] = psi_*qsF[f] + q2;
                    for (label k = 0; k < N; ++k)
                    {
                        shells_[k].boundaryFieldRef()[patchi][f] = xout[k];
                    }
                }
            }
            // Faktor keadaan rata-rata (< 1 utk self-limiting, TANPA clamp ke 1):
            // hanya koefisien linearisasi dqImplicitFactor().
            reduce(gSum, sumOp<scalar>());
            reduce(nCol, sumOp<label>());
            gammaMean_ = (nCol > 0) ? max(gSum/scalar(nCol), SMALL) : scalar(1);
            return;
        }

        {
            scalarField& qc = q.primitiveFieldRef();
            const scalarField& qsF = qStar.primitiveField();
            const scalarField& Tf = T_->primitiveField();
            forAll(qc, c)
            {
                for (label k = 0; k < N; ++k)
                {
                    old[k] = shells_[k].oldTime().primitiveField()[c];
                }
                qc[c] = solveColumnDarken
                (
                    old, qsF[c], Tf[c], dt,
                    Dsh, adiag, bdiag, cdiag, dwork, cprime, xout, gSum
                );
                ++nCol;
                for (label k = 0; k < N; ++k)
                {
                    shells_[k].primitiveFieldRef()[c] = xout[k];
                }
            }
        }

        forAll(q.boundaryField(), patchi)
        {
            scalarField& qp = q.boundaryFieldRef()[patchi];
            const scalarField& qsF = qStar.boundaryField()[patchi];
            const scalarField& Tf = T_->boundaryField()[patchi];
            forAll(qp, f)
            {
                for (label k = 0; k < N; ++k)
                {
                    old[k] = shells_[k].oldTime().boundaryField()[patchi][f];
                }
                scalar dummy = 0;   // muka batas tak masuk rata-rata Gamma
                qp[f] = solveColumnDarken
                (
                    old, qsF[f], Tf[f], dt,
                    Dsh, adiag, bdiag, cdiag, dwork, cprime, xout, dummy
                );
                for (label k = 0; k < N; ++k)
                {
                    shells_[k].boundaryFieldRef()[patchi][f] = xout[k];
                }
            }
        }

        // Rata-rata Gamma DIREDUKSI lintas-rank: dqImplicitFactor() harus
        // memberi angka yang sama di semua rank, kalau tidak linearisasi jadi
        // bergantung dekomposisi dan reproduktibilitas paralel rusak.
        reduce(gSum, sumOp<scalar>());
        reduce(nCol, sumOp<label>());
        gammaMean_ = (nCol > 0) ? max(gSum/scalar(nCol), scalar(1)) : scalar(1);

        return;
    }


    if (sphere_)
    {
        updateSphere(q, qStar, dt);
        return;
    }

    const label N = nShells_;
    const scalar h = delta_/scalar(N);
    const scalar r = Dc_*dt/(h*h);

    // --- Faktorisasi Thomas (koefisien identik tiap sel) dihitung SEKALI ---
    //   diag b_i : 1+r (backing) | 1+2r (interior) | 1+3r (permukaan)
    //   sub a_i = -r (i>=1) ; super c_i = -r (i<N-1)
    scalarList bdiag(N), cprime(N), invm(N);
    for (label i = 0; i < N; ++i)
    {
        bdiag[i] = (i == 0) ? (1.0 + r)
                 : (i == N - 1) ? (1.0 + 3.0*r)
                 : (1.0 + 2.0*r);
    }
    invm[0] = 1.0/bdiag[0];
    cprime[0] = (-r)*invm[0];
    for (label i = 1; i < N; ++i)
    {
        const scalar m = bdiag[i] - (-r)*cprime[i - 1];   // a_i=-r
        invm[i] = 1.0/m;
        cprime[i] = (i < N - 1 ? (-r) : 0.0)*invm[i];
    }

    // Buffer per kolom
    scalarList old(N), dwork(N), xout(N);

    // --- Internal field ---
    {
        scalarField& qc = q.primitiveFieldRef();
        const scalarField& qs = qStar.primitiveField();
        forAll(qc, c)
        {
            for (label k = 0; k < N; ++k)
            {
                old[k] = shells_[k].oldTime().primitiveField()[c];
            }
            qc[c] = solveColumn(old, qs[c], cprime, invm, r, dwork, xout);
            for (label k = 0; k < N; ++k)
            {
                shells_[k].primitiveFieldRef()[c] = xout[k];
            }
        }
    }

    // --- Boundary fields ---
    forAll(q.boundaryField(), patchi)
    {
        scalarField& qp = q.boundaryFieldRef()[patchi];
        const scalarField& qs = qStar.boundaryField()[patchi];
        forAll(qp, f)
        {
            for (label k = 0; k < N; ++k)
            {
                old[k] = shells_[k].oldTime().boundaryField()[patchi][f];
            }
            qp[f] = solveColumn(old, qs[f], cprime, invm, r, dwork, xout);
            for (label k = 0; k < N; ++k)
            {
                shells_[k].boundaryFieldRef()[patchi][f] = xout[k];
            }
        }
    }
}


Foam::scalar Foam::kineticClosures::intraGrainDiffusion::dqImplicitFactor
(
    const scalar dt
) const
{
    if (sphere_)
    {
        // Respons q-bar (berbobot volume) thd q* satuan, dari nol, tanpa clamp.
        // Self-limiting: D = Dc_ * faktor keadaan rata-rata langkah terakhir
        // (koefisien linearisasi saja). Situs permukaan: dq1/dq* = psi eksak.
        const label N = nShells_;
        const scalar h = delta_/scalar(N);
        const scalar r = Dc_*(selfLim_ ? gammaMean_ : scalar(1))*dt/(h*h);

        scalarList a(N), b(N), c(N), w(N), rhs(N, 0.0), cp(N), dw(N), x(N);
        scalar sRhs = 0;
        assembleSphere(N, r, a, b, c, sRhs, w);
        rhs[N - 1] = sRhs;
        const scalar gBulk = solveTridiagMean(a, b, c, rhs, w, false, cp, dw, x);
        return psi_ + (1.0 - psi_)*gBulk;
    }

    // g = d(q-bar)/d(q*): operator solve LINEAR pd RHS -> respons q-bar thd q*
    // = solveColumn(old=0, qs=1) (dari nol, clamp tak aktif). Faktorisasi Thomas
    // + back-subst identik update(), tanpa clamp (unit-response selalu positif).
    const label N = nShells_;
    const scalar h = delta_/scalar(N);
    // Dgn Darken, D bervariasi per sel sedangkan antarmuka menuntut SATU
    // skalar. Dipakai Gamma rata-rata volume dari langkah terakhir. Ini hanya
    // koefisien LINEARISASI: skemanya deferred-correction, sukunya saling
    // batal saat konvergen, jadi nilai hampiran mempengaruhi laju konvergensi
    // dan bukan jawaban yang dikonvergensi.
    const scalar r = Dc_*(darken_ ? gammaMean_ : scalar(1))*dt/(h*h);

    scalarList bdiag(N), cprime(N), invm(N), dwork(N);
    for (label i = 0; i < N; ++i)
    {
        bdiag[i] = (i == 0) ? (1.0 + r)
                 : (i == N - 1) ? (1.0 + 3.0*r)
                 : (1.0 + 2.0*r);
    }
    invm[0] = 1.0/bdiag[0];
    cprime[0] = (-r)*invm[0];
    for (label i = 1; i < N; ++i)
    {
        const scalar m = bdiag[i] - (-r)*cprime[i - 1];
        invm[i] = 1.0/m;
        cprime[i] = (i < N - 1 ? (-r) : 0.0)*invm[i];
    }

    // RHS unit q*: hanya suku permukaan 2r*qs (qs=1), old=0
    for (label i = 0; i < N; ++i) { dwork[i] = 0.0; }
    dwork[N - 1] = 2.0*r;
    dwork[0] *= invm[0];
    for (label i = 1; i < N; ++i)
    {
        dwork[i] = (dwork[i] + r*dwork[i - 1])*invm[i];
    }
    // Back-substitusi + rata-rata (tanpa clamp)
    scalar sum = dwork[N - 1];
    scalar xip = dwork[N - 1];
    for (label i = N - 2; i >= 0; --i)
    {
        const scalar xi = dwork[i] - cprime[i]*xip;
        sum += xi;
        xip = xi;
    }
    return sum/scalar(N);
}


// * * * * * * * * * * * * * * * Geometri bola  * * * * * * * * * * * * * * * //

void Foam::kineticClosures::intraGrainDiffusion::assembleSphere
(
    const label N,
    const scalar r,
    scalarList& a,
    scalarList& b,
    scalarList& c,
    scalar& sRhs,
    scalarList& w
)
{
    for (label i = 0; i < N; ++i)
    {
        const scalar ii = scalar(i);
        const scalar wi = ((ii + 1)*(ii + 1)*(ii + 1) - ii*ii*ii)/3.0;
        const scalar rL = (i == 0) ? 0.0 : r*ii*ii/wi;
        const scalar rR =
            (i == N - 1)
          ? 2.0*r*scalar(N)*scalar(N)/wi        // muka permukaan, jarak h/2
          : r*(ii + 1)*(ii + 1)/wi;

        w[i] = wi;
        a[i] = -rL;
        b[i] = 1.0 + rL + rR;
        c[i] = (i == N - 1) ? 0.0 : -rR;

        if (i == N - 1)
        {
            sRhs = rR;
        }
    }
}


void Foam::kineticClosures::intraGrainDiffusion::assembleSphereVar
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
)
{
    // Konduktansi muka SERI -> harmonik (CATATAN 23.94); aman thd D -> 0.
    auto hmean = [](const scalar Da, const scalar Db)
    {
        return 2.0*Da*Db/max(Da + Db, VSMALL);
    };

    for (label i = 0; i < N; ++i)
    {
        const scalar ii = scalar(i);
        const scalar wi = ((ii + 1)*(ii + 1)*(ii + 1) - ii*ii*ii)/3.0;
        const scalar rL =
            (i == 0) ? 0.0 : hmean(Dsh[i - 1], Dsh[i])*dth2*ii*ii/wi;
        const scalar rR =
            (i == N - 1)
          ? 2.0*hmean(Dsh[N - 1], Dsurf)*dth2*scalar(N)*scalar(N)/wi
          : hmean(Dsh[i], Dsh[i + 1])*dth2*(ii + 1)*(ii + 1)/wi;

        w[i] = wi;
        a[i] = -rL;
        b[i] = 1.0 + rL + rR;
        c[i] = (i == N - 1) ? 0.0 : -rR;

        if (i == N - 1)
        {
            sRhs = rR;
        }
    }
}


Foam::scalar Foam::kineticClosures::intraGrainDiffusion::solveColumnSphereVar
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
) const
{
    const label N = nShells_;
    const scalar h = delta_/scalar(N);

    // D per kulit dari loading bulk LOKAL, di-LAG pada oldTime (Picard).
    scalar gAcc = 0;
    for (label i = 0; i < N; ++i)
    {
        const scalar f = stateFactor(old[i], Tcell);
        Dsh[i] = Dc_*f;
        gAcc += f;
    }
    gSum += gAcc/scalar(N);

    const scalar Dsurf = Dc_*stateFactor(qsBulk, Tcell);

    scalar sRhs = 0;
    assembleSphereVar(N, dt/(h*h), Dsh, Dsurf, a, b, c, sRhs, w);

    for (label i = 0; i < N; ++i) { rhs[i] = old[i]; }
    rhs[N - 1] += sRhs*qsBulk;

    return solveTridiagMean(a, b, c, rhs, w, true, cp, dw, x);
}


Foam::scalar Foam::kineticClosures::intraGrainDiffusion::solveTridiagMean
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
)
{
    const label N = b.size();

    scalar m = b[0];
    cp[0] = c[0]/m;
    dw[0] = rhs[0]/m;
    for (label i = 1; i < N; ++i)
    {
        m = b[i] - a[i]*cp[i - 1];
        cp[i] = c[i]/m;
        dw[i] = (rhs[i] - a[i]*dw[i - 1])/m;
    }

    scalar sumW = 0, sumWX = 0;

    x[N - 1] = clamp ? max(dw[N - 1], scalar(0)) : dw[N - 1];
    sumW += w[N - 1];
    sumWX += w[N - 1]*x[N - 1];

    for (label i = N - 2; i >= 0; --i)
    {
        const scalar xi = dw[i] - cp[i]*x[i + 1];
        x[i] = clamp ? max(xi, scalar(0)) : xi;
        sumW += w[i];
        sumWX += w[i]*x[i];
    }

    return sumWX/sumW;
}


void Foam::kineticClosures::intraGrainDiffusion::updateSphere
(
    volScalarField& q,
    const volScalarField& qStar,
    const scalar dt
)
{
    const label N = nShells_;
    const scalar h = delta_/scalar(N);
    const scalar r = Dc_*dt/(h*h);

    // Koefisien identik tiap sel (D konstan) -> dirakit SEKALI per langkah.
    scalarList a(N), b(N), c(N), w(N), rhs(N), cp(N), dw(N), x(N);
    scalar sRhs = 0;
    assembleSphere(N, r, a, b, c, sRhs, w);

    // --- Internal field ---
    {
        scalarField& qc = q.primitiveFieldRef();
        const scalarField& qs = qStar.primitiveField();
        forAll(qc, ci)
        {
            for (label k = 0; k < N; ++k)
            {
                rhs[k] = shells_[k].oldTime().primitiveField()[ci];
            }
            rhs[N - 1] += sRhs*qs[ci];

            qc[ci] = solveTridiagMean(a, b, c, rhs, w, true, cp, dw, x);

            for (label k = 0; k < N; ++k)
            {
                shells_[k].primitiveFieldRef()[ci] = x[k];
            }
        }
    }

    // --- Boundary fields ---
    forAll(q.boundaryField(), patchi)
    {
        scalarField& qp = q.boundaryFieldRef()[patchi];
        const scalarField& qs = qStar.boundaryField()[patchi];
        forAll(qp, f)
        {
            for (label k = 0; k < N; ++k)
            {
                rhs[k] = shells_[k].oldTime().boundaryField()[patchi][f];
            }
            rhs[N - 1] += sRhs*qs[f];

            qp[f] = solveTridiagMean(a, b, c, rhs, w, true, cp, dw, x);

            for (label k = 0; k < N; ++k)
            {
                shells_[k].boundaryFieldRef()[patchi][f] = x[k];
            }
        }
    }
}


// ************************************************************************* //
