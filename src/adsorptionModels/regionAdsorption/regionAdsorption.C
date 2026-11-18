/*---------------------------------------------------------------------------*\
  regionAdsorption.C
\*---------------------------------------------------------------------------*/

#include "regionAdsorption.H"
#include "fvm.H"
#include "fvc.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Helper (file-local) * * * * * * * * * * * * * //

namespace Foam
{

// Evaluasi isoterm+campuran pada satu titik (T, C[]) menghasilkan koefisien
// yang dipakai persamaan species & energi. Lihat regionAdsorption.H.
//   qStarOut[i]  = f_a q_i*                         [mol/kg]
//   effStorOut[i]= eps_p + rho_s M_i dq_i/dC_i      [-]
//   STcoefOut[i] = rho_s M_i dq_i/dT|_C             [kg/m3/K]
//   crossOut(i,j)= rho_s M_i dq_i/dC_j              [-]   (i!=j)
static inline void evalPoint
(
    const PtrList<isothermModel>& isos,
    const mixingRule& mixing,
    const scalarList& M,
    const scalar rhoS,
    const scalar epsPstor,   // pori TOTAL: evalPoint hanya memakainya utk effStor
    const scalar fa,
    const scalar Tin,
    const scalarList& Cvals,
    scalarList& qStarOut,
    scalarList& effStorOut,
    scalarList& STcoefOut,
    scalarSquareMatrix& crossOut,
    // buffers (dialokasikan pemanggil, ukuran N)
    scalarList& p,
    scalarList& qc,
    scalarList& dqdT_p,
    scalarSquareMatrix& dqdp
)
{
    const label N = isos.size();
    const scalar T = max(Tin, scalar(1));
    const scalar R = adsorption::Rgas;

    forAll(isos, j)
    {
        p[j] = adsorption::partialPressureFromConc(Cvals[j], M[j], T);
    }

    mixing.evaluate(isos, p, T, qc, dqdp, dqdT_p);

    for (label i = 0; i < N; ++i)
    {
        qStarOut[i] = fa*qc[i];

        // dq_i/dT|_C = fa ( dq_i/dT|_p + sum_k dq_i/dp_k * C_k R / M_k )
        scalar dqdT_tot = dqdT_p[i];
        for (label k = 0; k < N; ++k)
        {
            dqdT_tot += dqdp(i, k)*Cvals[k]*R/M[k];
        }
        dqdT_tot *= fa;
        STcoefOut[i] = rhoS*M[i]*dqdT_tot;

        for (label j = 0; j < N; ++j)
        {
            // dq_i/dC_j = fa dq_i/dp_j * R T / M_j
            const scalar dqdC = fa*dqdp(i, j)*R*T/M[j];
            crossOut(i, j) = rhoS*M[i]*dqdC;
        }
        effStorOut[i] = epsPstor + crossOut(i, i);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::regionAdsorption::regionAdsorption
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    mesh_(mesh),
    thermo_(thermo),
    dict_
    (
        IOobject
        (
            "adsorptionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    species_(dict_.get<wordList>("species")),
    M_(species_.size(), 0.0),
    deltaH_(species_.size(), 0.0),
    DsRef_(species_.size(), 0.0),
    bosanquet_(species_.size(), false),
    Dbulk_(species_.size(), 0.0),
    pRefB_(species_.size(), 0.0),
    TRefB_(species_.size(), 0.0),
    Dknud_(species_.size(), 0.0),
    DknudGiven_(species_.size(), false),
    rPore_(species_.size(), 0.0),
    epsilonP_(dict_.get<scalar>("epsilonP")),
    epsilonPtot_(dict_.getOrDefault<scalar>("epsilonPtotal", epsilonP_)),
    rhoS_(dict_.get<scalar>("rhoS")),
    fActive_(dict_.getOrDefault<scalar>("fActive", dict_.getOrDefault<scalar>("fMono", 1.0))),
    isothermalAdsorption_(dict_.getOrDefault<bool>("isothermalAdsorption", false)),
    nNewton_(dict_.getOrDefault<label>("nNewton", 50)),
    newtonTol_(dict_.getOrDefault<scalar>("newtonTol", 1e-8)),
    isotherms_(species_.size()),
    mixing_(mixingRule::New(dict_)),
    kinetics_(species_.size()),
    C_(species_.size()),
    q_(species_.size()),
    qStar_(species_.size()),
    Ds_(species_.size()),
    effStor_(species_.size()),
    Sads_
    (
        IOobject
        (
            "Sads",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimEnergy/dimVolume/dimTime, 0)
    )
{
    Info<< "Constructing adsorption for region " << mesh.name()
        << " : species " << flatOutput(species_)
        << " (mixingRule " << mixing_->type() << ")" << nl
        << "  epsilonP=" << epsilonP_
        << " (transport) epsilonPtotal=" << epsilonPtot_ << " (simpanan)"
        << " rhoS=" << rhoS_
        << " fActive=" << fActive_ << endl;

    // Pori total tak mungkin lebih kecil daripada makropori: yang kedua adalah
    // bagian dari yang pertama. Kalau terbalik, dictnya keliru dan diamkan saja
    // akan memberi simpanan gas pori yang terlalu kecil tanpa jejak.
    if (epsilonPtot_ < epsilonP_)
    {
        FatalIOErrorInFunction(dict_)
            << "epsilonPtotal (" << epsilonPtot_ << ") < epsilonP ("
            << epsilonP_ << ")." << nl
            << "    epsilonP     = makropori, rongga yang MENGANGKUT" << nl
            << "    epsilonPtotal= pori TOTAL, rongga yang MENYIMPAN" << nl
            << "    Yang kedua memuat yang pertama, jadi ia tak boleh lebih kecil."
            << exit(FatalIOError);
    }

    const dimensionSet dimLoad(dimMoles/dimMass);

    forAll(species_, i)
    {
        const word& sp = species_[i];
        const dictionary& sd = dict_.subDict(sp);

        M_[i]      = sd.get<scalar>("molWeight")/1000.0;   // g/mol -> kg/mol
        DsRef_[i]  = sd.get<scalar>("Ds");
        // deltaH_[i] ditetapkan SETELAH isoterm terbentuk -- lihat resolveDeltaH()

        if (sd.found("bosanquet"))
        {
            const dictionary& bd = sd.subDict("bosanquet");
            bosanquet_[i] = bd.getOrDefault<bool>("active", false);
            Dbulk_[i]     = bd.getOrDefault<scalar>("Dbulk", DsRef_[i]);
            pRefB_[i]     = bd.getOrDefault<scalar>("pRef", 101325.0);
            TRefB_[i]     = bd.getOrDefault<scalar>("TRef", 298.0);
            DknudGiven_[i] = bd.found("Dknudsen");
            Dknud_[i]     = bd.getOrDefault<scalar>("Dknudsen", DsRef_[i]);
            rPore_[i]     = bd.getOrDefault<scalar>("rPore", 0.0);
        }

        isotherms_.set(i, isothermModel::New(sp, sd.subDict("isotherm")));

        // Panas adsorpsi: sumber tunggal + pemeriksaan silang. WAJIB setelah
        // isoterm terbentuk, karena isotermlah sumbernya.
        resolveDeltaH(i, sd);

        kinetics_.set(i, kineticClosure::New(mesh, sp, sd.subDict("kinetics")));

        // Closure yang memakai koreksi Darken butuh isoterm (utk Gamma), suhu,
        // dan fActive (krn q disimpan sbg fActive*q_isoterm). Closure lain
        // mengabaikan panggilan ini.
        kinetics_[i].setIsothermContext(isotherms_[i], thermo_.T(), fActive_);

        // --- field C.<sp> [kg/m3] ---
        C_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "C." + sp,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        // --- field Ds.<sp> [m2/s] ---
        Ds_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Ds." + sp,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE      // P1: dapat didiagnosis dari disk
                ),
                mesh,
                dimensionedScalar("Ds", dimArea/dimTime, DsRef_[i])
            )
        );

        // --- field effStor.<sp> [-] : eps_p + rho_s M dq*/dC ---
        // Nilai awal epsilonPtot_ (belum ada muatan), diisi ulang tiap solve().
        // NO_WRITE: ia besaran turunan, bukan keadaan; menulisnya tiap waktu
        // hanya membesarkan arsip tanpa menambah keterangan.
        effStor_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "effStor." + sp,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("eS", dimless, epsilonPtot_)
            )
        );

        // --- field qStar.<sp> [mol/kg] ---
        qStar_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "qStar." + sp,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE      // P1: dapat didiagnosis dari disk
                ),
                mesh,
                dimensionedScalar("qStar", dimLoad, 0)
            )
        );

        // --- field q.<sp> [mol/kg] (READ_IF_PRESENT) ---
        IOobject qIO
        (
            "q." + sp,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        q_.set
        (
            i,
            new volScalarField
            (
                qIO,
                mesh,
                dimensionedScalar("q", dimLoad, 0)
            )
        );

        q_[i].oldTime();
    }

    // P1 -- penjaga dijalankan SEBELUM langkah waktu pertama, supaya
    // ketidakcocokan material ketahuan saat konstruksi, bukan setelah satu run
    // lewat neraca massa antarmuka.
    checkSolidDensity();
    checkBosanquet();

    initialiseLoading();
}


// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::regionAdsorption::resolveDeltaH(const label i, const dictionary& sd)
{
    // Panas adsorpsi hidup di DUA entri yang menggambarkan besaran fisis SAMA:
    //   `deltaH`             -> source energi solid  S_ads = dH rho_s dq/dt
    //   `isotherm.deltaH0`   -> ketergantungan suhu  b(T) = b0 exp[dH/(R T0)(T0/T - 1)]
    // Kelas jebakan yang sama dgn rho_s (checkSolidDensity) dan eps_p. Selama
    // materialnya satu dan tetap ia tak menggigit; begitu pergantian material
    // jadi rutin, memperbarui SATU tempat saja memberi isoterm dan neraca
    // energi dgn dH berbeda -- bergeser diam-diam, tampak fisik, dan hanya
    // ketahuan lewat suhu gelombang adsorpsi.
    //
    // Aturan di sini:
    //   isothermalAdsorption  -> dH = 0, S_ads mati (jalur sengaja)
    //   `deltaH` ditulis      -> dipakai, TAPI diperiksa silang thd isoterm
    //   `deltaH` tak ditulis  -> diambil dari isoterm (sumber tunggal)
    //   keduanya tak ada      -> FatalIOError, bukan diam-diam nol
    const word& sp = species_[i];
    const scalar dHiso = isotherms_[i].heatOfAdsorption();
    const bool   given = sd.found("deltaH");
    const scalar tol   = 1e-3;      // 0,1% -- ambang sama dgn checkSolidDensity()

    if (isothermalAdsorption_)
    {
        deltaH_[i] = 0;

        Info<< "  [P1] " << sp << ": isothermalAdsorption AKTIF -> deltaH = 0,"
            << " S_ads dimatikan." << nl
            << "       Kasus ini TIDAK memodelkan panas adsorpsi. Jangan"
            << " kutip suhunya sebagai hasil fisis." << endl;
        return;
    }

    if (given)
    {
        deltaH_[i] = sd.get<scalar>("deltaH");

        if (mag(dHiso) > SMALL)
        {
            const scalar dev = mag(deltaH_[i] - dHiso)/max(mag(dHiso), SMALL);

            if (dev > tol)
            {
                FatalIOErrorInFunction(sd)
                    << "Panas adsorpsi TIDAK COCOK untuk species " << sp
                    << " di region " << mesh_.name() << ':' << nl
                    << "    `deltaH` (source energi)      : "
                    << deltaH_[i] << " J/mol" << nl
                    << "    `isotherm.deltaH0` (b(T))     : "
                    << dHiso << " J/mol" << nl
                    << "    simpangan relatif             : " << dev
                    << "  (ambang " << tol << ')' << nl
                    << nl
                    << "Keduanya adalah panas adsorpsi isosterik yang SAMA."
                    << " Kalau berbeda, isoterm dan neraca energi memakai"
                    << " termodinamika yang berbeda." << nl
                    << nl
                    << "PERBAIKAN: hapus `deltaH` -- ia otomatis mengambil"
                    << " nilai dari isoterm. Tulis hanya bila memang perlu"
                    << " berbeda, dan itu menuntut `isothermalAdsorption yes`."
                    << nl
                    << exit(FatalIOError);
            }

            Info<< "  [P1] " << sp << ": deltaH cocok, dict " << deltaH_[i]
                << " vs isoterm " << dHiso << " J/mol" << endl;
        }
        else
        {
            // Isoterm tak menyatakan satu nilai (mis. DualSiteLangmuir).
            // Dict jadi satu-satunya sumber; tak ada yang bisa disilangkan.
            Info<< "  [P1] " << sp << ": deltaH = " << deltaH_[i]
                << " J/mol dari dict; isoterm '" << isotherms_[i].type()
                << "' tidak menyatakan nilai tunggal -> TAK ADA pemeriksaan"
                << " silang." << endl;
        }
        return;
    }

    if (mag(dHiso) > SMALL)
    {
        deltaH_[i] = dHiso;

        Info<< "  [P1] " << sp << ": deltaH bersumber dari isotherm.deltaH0 = "
            << deltaH_[i] << " J/mol (entri `deltaH` tidak ditulis -- benar)"
            << endl;
        return;
    }

    FatalIOErrorInFunction(sd)
        << "Panas adsorpsi tidak diketahui untuk species " << sp
        << " di region " << mesh_.name() << ':' << nl
        << "    entri `deltaH` tidak ditulis, dan isoterm '"
        << isotherms_[i].type() << "' tidak menyatakan nilai tunggal." << nl
        << nl
        << "Isoterm bersitus-jamak (mis. DualSiteLangmuir) punya lebih dari"
        << " satu energi situs, sehingga panas adsorpsi efektifnya bergantung"
        << " keadaan. Solver memakai SATU nilai per species, jadi pilihannya"
        << " harus dinyatakan tangan dan terekam di dict." << nl
        << nl
        << "PERBAIKAN: tulis `deltaH <J/mol>;` di sub-dict species " << sp
        << ", atau pakai isoterm bersitus-tunggal."
        << nl
        << exit(FatalIOError);
}


void Foam::regionAdsorption::checkSolidDensity() const
{
    // P1. rho_s hidup di DUA berkas: `rhoS` di adsorptionProperties dan `rho`
    // di thermophysicalProperties. Kelas jebakan yang sama dgn eps_p (CATATAN
    // 23.80[1]); bedanya eps_p sudah dijadikan bersumber-tunggal, sedangkan
    // rho_s tidak bisa -- thermo memang harus punya densitasnya sendiri.
    // Karena itu ia DIJAGA, bukan dihapus.
    const volScalarField rhoT(thermo_.rho());

    const scalar rhoMin = min(rhoT).value();
    const scalar rhoMax = max(rhoT).value();
    const scalar tol    = 1e-3;      // 0,1% -- keduanya seharusnya angka SAMA

    const scalar dMin = mag(rhoMin - rhoS_)/max(mag(rhoS_), SMALL);
    const scalar dMax = mag(rhoMax - rhoS_)/max(mag(rhoS_), SMALL);

    if (max(dMin, dMax) > tol)
    {
        FatalErrorInFunction
            << "rho_s TIDAK COCOK antara dua berkas di region "
            << mesh_.name() << ':' << nl
            << "    `rhoS` di adsorptionProperties        : " << rhoS_ << nl
            << "    thermo rho (thermophysicalProperties) : "
            << rhoMin << " .. " << rhoMax << nl
            << "    simpangan relatif terbesar            : "
            << max(dMin, dMax) << "  (ambang " << tol << ')' << nl
            << nl
            << "Keduanya adalah densitas dinding tampak yang SAMA. Kalau"
            << " berbeda, kapasitas adsorpsi dan inersia termal memakai"
            << " densitas yang berbeda." << nl
            << "PERBAIKAN: samakan `rhoS` dan `equationOfState { rho ... }`."
            << nl
            << exit(FatalError);
    }

    Info<< "  [P1] rho_s cocok: dict " << rhoS_ << " vs thermo "
        << rhoMin << " .. " << rhoMax << endl;
}


void Foam::regionAdsorption::checkBosanquet()
{
    // P1. D_Knudsen dari GEOMETRI PORI, bukan dari angka yang diketik:
    //     D_K = (2/3) r_pore sqrt(8 R T / (pi M))
    // Lalu Bosanquet pada (pRef,TRef) dibandingkan dgn `Ds` yang dipakai.
    // Selisihnya BUKAN sekadar cacat: ia MENGUKUR tortuositas tersirat
    //     tau_tersirat = D_Bosanquet(pRef,TRef) / Ds
    // Itu yang menaikkan penjaga ini dari kebersihan kode menjadi pengukuran
    // (CATATAN 23.78[7]): `Ds` proyek ini mewarisi tau = 2 dari mullite 2024,
    // dan dgn geometri pori material lain angka itu tidak lagi berlaku.
    forAll(species_, i)
    {
        if (rPore_[i] <= 0)
        {
            if (bosanquet_[i])
            {
                Info<< "  [P1] " << species_[i] << ": Bosanquet AKTIF tanpa"
                    << " `rPore` -> D_Knudsen dipakai apa adanya ("
                    << Dknud_[i]
                    << (DknudGiven_[i] ? " m2/s dari dict" : " m2/s = NILAI BAWAAN (=Ds), bukan dari dict")
                    << "), tortuositas tersirat TIDAK"
                    << " terukur." << endl;
            }
            continue;
        }

        const scalar Tref = TRefB_[i] > 0 ? TRefB_[i] : 298.0;

        // v_bar = sqrt(8 R T / (pi M)) ; R = 8.31446 J/(mol K), M [kg/mol]
        const scalar vBar =
            Foam::sqrt(8.0*8.314462618*Tref/(constant::mathematical::pi*M_[i]));

        const scalar DkCalc = (2.0/3.0)*rPore_[i]*vBar;
        const scalar Dbos   = 1.0/(1.0/max(Dbulk_[i], SMALL)
                                 + 1.0/max(DkCalc, SMALL));
        const scalar tauImp = Dbos/max(DsRef_[i], SMALL);

        // `word` tidak boleh memuat spasi (dipotong / assert), jadi teksnya
        // dirakit sebagai string lebih dulu, bukan lewat ternary di dalam <<.
        const Foam::string dknudTxt =
            DknudGiven_[i]
          ? Foam::string(Foam::name(Dknud_[i]))
          : Foam::string("tidak ditulis");

        Info<< "  [P1] " << species_[i] << " Bosanquet @ (" << pRefB_[i]
            << " Pa, " << Tref << " K):" << nl
            << "       r_pore     = " << rPore_[i] << " m" << nl
            << "       D_K hitung = " << DkCalc << " m2/s   (dict: "
            << dknudTxt.c_str()      // .c_str(): Foam::string di-stream BER-KUTIP
            << ')' << nl
            << "       D_Bosanquet= " << Dbos << " m2/s" << nl
            << "       Ds dipakai = " << DsRef_[i] << " m2/s" << nl
            << "       ** tortuositas TERSIRAT = " << tauImp << " **" << endl;

        // Ambang: tau tersirat harus berada di rentang fisik. tau < 1 mustahil
        // (jalur pori tak bisa lebih pendek dari garis lurus); tau > 10 berarti
        // Ds dan geometri pori bercerita tentang material yang berbeda.
        if (bosanquet_[i] && (tauImp < 1.0 || tauImp > 10.0))
        {
            FatalErrorInFunction
                << "tortuositas tersirat " << tauImp << " di luar rentang"
                << " fisik [1,10] untuk " << species_[i] << '.' << nl
                << "Bosanquet AKTIF, jadi Ds(p,T) dihitung dari geometri pori"
                << " yang tidak konsisten dgn `Ds` = " << DsRef_[i] << nl
                << "PERBAIKAN: periksa `rPore`, `Dbulk`, dan `Ds` -- salah satu"
                << " milik material lain." << nl
                << exit(FatalError);
        }

        // D_K terhitung MENGGANTIKAN angka dict saat Bosanquet aktif.
        if (bosanquet_[i])
        {
            Dknud_[i] = DkCalc;
        }
    }
}


void Foam::regionAdsorption::updateDs()
{
    // Bosanquet: 1/Ds = 1/D_gas(p,T) + 1/D_Knudsen(T)
    //   D_gas ~ T^1.75 / p   (Fuller; p dari pers. gas, T dari solver)
    //   D_K   ~ sqrt(T)      (teori kinetik)
    // Keduanya dinormalkan pada (pRef,TRef). Saat `bosanquet.active false`
    // (kasus adsorpsi 1 atm), Ds TETAP konstan = DsRef -- jadi seluruh hasil
    // kampanye tidak tersentuh oleh penambahan ketergantungan T ini.
    const volScalarField* pPtr = mesh_.findObject<volScalarField>("p");
    const volScalarField& T = thermo_.T();

    forAll(species_, i)
    {
        if (!bosanquet_[i] || !pPtr)
        {
            Ds_[i] == dimensionedScalar("Ds", dimArea/dimTime, DsRef_[i]);
            continue;
        }

        const volScalarField& p = *pPtr;
        const scalar Tref = TRefB_[i] > 0 ? TRefB_[i] : 298.0;
        const scalar DbosRef =
            1.0/(1.0/max(Dbulk_[i], SMALL) + 1.0/max(Dknud_[i], SMALL));

        scalarField& dsc = Ds_[i].primitiveFieldRef();
        const scalarField& pc = p.primitiveField();
        const scalarField& Tc = T.primitiveField();

        forAll(dsc, c)
        {
            const scalar pp = max(pc[c], SMALL);
            const scalar tt = max(Tc[c], SMALL);
            const scalar tr = tt/Tref;

            const scalar Dgas = Dbulk_[i]*(pRefB_[i]/pp)*Foam::pow(tr, 1.75);
            const scalar Dk   = max(Dknud_[i], SMALL)*Foam::sqrt(tr);

            // Paper B (2026-10-03): dinormalkan ke Ds pada (pRef, TRef), yaitu
            // dibagi tortuositas tersirat tau = D_Bos(pRef,TRef)/Ds (checkBosanquet).
            // Dulu D_Bos dipakai langsung -> pada pRef Ds = D_Bos, bukan `Ds` dict
            // (mis. 1,39e-5 vs 8e-6 untuk mullite). Tak ada run lama dgn Bosanquet aktif.
            dsc[c] = DsRef_[i]*(1.0/(1.0/Dgas + 1.0/Dk))/DbosRef;
        }
        Ds_[i].correctBoundaryConditions();
    }
}


void Foam::regionAdsorption::initialiseLoading()
{
    const label N = species_.size();
    const volScalarField& T = thermo_.T();

    // Buffer
    scalarList Cvals(N), qStarOut(N), effStorOut(N), STcoefOut(N);
    scalarList p(N), qc(N), dqdT_p(N);
    scalarSquareMatrix crossOut(N, Zero), dqdp(N, Zero);

    forAll(species_, i)
    {
        // Hanya isi bila q tidak terbaca dari disk (nilai awal nol)
        if (gMax(mag(q_[i].primitiveField())()) > SMALL)
        {
            continue;
        }

        const scalarField& Tc = T.primitiveField();
        scalarField& qiC = q_[i].primitiveFieldRef();

        forAll(qiC, c)
        {
            forAll(species_, j)
            {
                Cvals[j] = C_[j].primitiveField()[c];
            }
            evalPoint
            (
                isotherms_, *mixing_, M_, rhoS_, epsilonPtot_, fActive_,
                Tc[c], Cvals, qStarOut, effStorOut, STcoefOut, crossOut,
                p, qc, dqdT_p, dqdp
            );
            qiC[c] = qStarOut[i];
        }
        q_[i].correctBoundaryConditions();
        q_[i].oldTime() = q_[i];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionAdsorption::storeOldTime()
{
    forAll(species_, i)
    {
        q_[i].storeOldTime();
        kinetics_[i].storeOldTime();
    }
}


void Foam::regionAdsorption::correct
(
    const label nCouple,
    speciesImplicitCoupler* coupler
)
{
    const label N = species_.size();
    const volScalarField& T = thermo_.T();
    const scalar dt = mesh_.time().deltaTValue();

    const dimensionedScalar dtDim("dt", dimTime, max(dt, SMALL));
    const dimensionedScalar rhoSDim("rhoS", dimDensity, rhoS_);
    const dimensionedScalar epsPDim("epsP", dimless, epsilonP_);       // transport
    const dimensionedScalar epsSDim("epsS", dimless, epsilonPtot_);    // simpanan

    updateDs();

    // Klasifikasi closure tiap species
    boolList algebraic(N);
    forAll(species_, i)
    {
        algebraic[i] = kinetics_[i].isAlgebraic();
    }

    // Buffer skalar per sel (evalPoint)
    scalarList Cvals(N), qStarOut(N), effStorOut(N), STcoefOut(N);
    scalarList pbuf(N), qcbuf(N), dqdT_pbuf(N);
    scalarSquareMatrix crossOut(N, Zero), dqdp(N, Zero);

    // Field koefisien silang cross_ij = rho_s M_i dq*_i/dC_j  [-]
    List<PtrList<volScalarField>> cross(N);
    forAll(species_, i)
    {
        cross[i].setSize(N);
        forAll(species_, j)
        {
            cross[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        "cross." + species_[i] + "." + species_[j],
                        mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE, false
                    ),
                    mesh_,
                    dimensionedScalar("c", dimless, 0)
                )
            );
        }
    }

    // Salinan C pada awal tiap iterasi Newton (C^k) -> residu + line-search
    PtrList<volScalarField> Ckeep(N);

    // ------------------------------------------------------------------
    //  Newton fully-coupled pada STORAGE KEKAL EKSAK (Bab 7.1 / 10.6):
    //      d/dt[ eps_p C_i + rho_s M_i q*_i(C,T) ] = div(eps_p Ds_i grad C_i)
    //  Jacobian J_ij = eps_p delta_ij + cross_ij (cross_ij = rho_s M_i dq*_i/dC_j).
    //  Linearisasi full-form per iterasi k, species i (algebraic):
    //      eps_p ddt(C_i) + Sp(cross_ii/dt) C_i - lap(C_i)
    //        = cross_ii C_i^k/dt - rho_s M_i (q*_i^k - q_i_old)/dt - offdiag
    //  Saat C->C^k suku cross_ii batal => persamaan KEKAL EKSAK:
    //      eps_p ddt(C) + rho_s M (q*(C)-q_old)/dt = lap (= Js).
    //  Konservasi TIDAK bergantung jumlah iterasi (storage = selisih sekan
    //  q*(C)-q_old, bukan tangent*dC). Iterasi hanya mengejar nonlinearitas
    //  q*(C); Newton -> konvergensi kuadratik (beda dgn Picard ter-relaksasi
    //  lama yg macet). Termal masuk lewat q*(C,T) (T kini) vs q_old -> tak ada
    //  double-count STcoef. Off-diagonal multikomponen terkopel Gauss-Seidel;
    //  single-species otomatis implisit penuh.
    // ------------------------------------------------------------------
    const label maxIter = max(nCouple, nNewton_);
    label iter = 0;
    bool newtonConverged = false;
    for (iter = 0; iter < maxIter; ++iter)
    {
        forAll(species_, i)
        {
            if (Ckeep.set(i)) { Ckeep[i] == C_[i]; }
            else { Ckeep.set(i, new volScalarField(C_[i])); }
        }

        // Evaluasi isoterm di C^k -> qStar_, cross_ij (internal + boundary)
        {
            const scalarField& Tc = T.primitiveField();
            forAll(Tc, c)
            {
                forAll(species_, j) { Cvals[j] = C_[j].primitiveField()[c]; }
                evalPoint
                (
                    isotherms_, *mixing_, M_, rhoS_, epsilonPtot_, fActive_,
                    Tc[c], Cvals, qStarOut, effStorOut, STcoefOut, crossOut,
                    pbuf, qcbuf, dqdT_pbuf, dqdp
                );
                forAll(species_, i)
                {
                    qStar_[i].primitiveFieldRef()[c] = qStarOut[i];
                    effStor_[i].primitiveFieldRef()[c] = effStorOut[i];
                    forAll(species_, j)
                    {
                        cross[i][j].primitiveFieldRef()[c] = crossOut(i, j);
                    }
                }
            }
        }
        forAll(mesh_.boundary(), patchi)
        {
            const scalarField& Tp = T.boundaryField()[patchi];
            forAll(Tp, f)
            {
                forAll(species_, j)
                {
                    Cvals[j] = C_[j].boundaryField()[patchi][f];
                }
                evalPoint
                (
                    isotherms_, *mixing_, M_, rhoS_, epsilonPtot_, fActive_,
                    Tp[f], Cvals, qStarOut, effStorOut, STcoefOut, crossOut,
                    pbuf, qcbuf, dqdT_pbuf, dqdp
                );
                forAll(species_, i)
                {
                    qStar_[i].boundaryFieldRef()[patchi][f] = qStarOut[i];
                    forAll(species_, j)
                    {
                        cross[i][j].boundaryFieldRef()[patchi][f] =
                            crossOut(i, j);
                    }
                }
            }
        }

        // Integrasi closure finite-rate (q jadi state)
        forAll(species_, i)
        {
            if (!algebraic[i])
            {
                kinetics_[i].update(q_[i], qStar_[i], dt);
            }
        }

        // Solve C tiap species
        forAll(species_, i)
        {
            const dimensionedScalar MiDim("M", dimMass/dimMoles, M_[i]);

            volScalarField epsDs
            (
                IOobject
                (
                    "epsDs." + species_[i],
                    mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE, false
                ),
                epsPDim*Ds_[i]
            );

            if (algebraic[i])
            {
                // Kopling silang Newton (Gauss-Seidel): sum_{j!=i} cross_ij dC_j/dt
                volScalarField offdiag
                (
                    IOobject
                    (
                        "offdiag." + species_[i],
                        mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE, false
                    ),
                    mesh_,
                    dimensionedScalar("o", dimDensity/dimTime, 0)
                );
                forAll(species_, j)
                {
                    if (j != i)
                    {
                        offdiag += cross[i][j]*(C_[j] - Ckeep[j])/dtDim;
                    }
                }

                fvScalarMatrix CiEqn
                (
                    epsSDim*fvm::ddt(C_[i])
                  + fvm::Sp(cross[i][i]/dtDim, C_[i])
                  - fvm::laplacian(epsDs, C_[i])
                 ==
                    cross[i][i]*Ckeep[i]/dtDim
                  - rhoSDim*MiDim*(qStar_[i] - q_[i].oldTime())/dtDim
                  - offdiag
                );
                if (poreFlux_)
                {
                    CiEqn += fvm::div(*poreFlux_, C_[i], "div(phiDarcy,C)");
                    if (poreSink_)
                    {
                        CiEqn += fvm::Sp(*poreSink_, C_[i]);
                    }
                    if (poreSrc_ && poreSrc_->set(i))
                    {
                        CiEqn -= (*poreSrc_)[i];
                    }
                }
                if (coupler && coupler->coupled(i))
                {
                    coupler->solve(i, CiEqn);
                }
                else
                {
                    CiEqn.solve();
                }
            }
            else
            {
                // Finite-rate: storage eps_p + SINK sekan -rho_s M (q-q_old)/dt,
                // dgn ketergantungan q_ pd C DILINEARISASI IMPLISIT (fix konservasi
                // 2026-07-01 opsi-b, CATATAN BAGIAN 14). dS/dC = g*cross_ii, g =
                // dq_/dqStar (dqImplicitFactor). Suku fvm::Sp(g*cross_ii/dt) +
                // deferred g*cross_ii*Ckeep/dt saling BATAL saat konvergen (C->Ckeep)
                // -> sink sekan EKSAK (konservatif); saat iterasi -> diagonal positif
                // cegah C<0 -> TAK butuh redaman line-search -> tak bocor pd closure
                // KAKU (HSDM/LDF k-tinggi). g=0 (default) -> jalur eksplisit lama.
                const scalar g = kinetics_[i].dqImplicitFactor(dt);
                fvScalarMatrix CiEqn
                (
                    epsSDim*fvm::ddt(C_[i])
                  + fvm::Sp(g*cross[i][i]/dtDim, C_[i])
                  - fvm::laplacian(epsDs, C_[i])
                 ==
                    g*cross[i][i]*Ckeep[i]/dtDim
                  - rhoSDim*MiDim*(q_[i] - q_[i].oldTime())/dtDim
                );
                if (poreFlux_)
                {
                    CiEqn += fvm::div(*poreFlux_, C_[i], "div(phiDarcy,C)");
                    if (poreSink_)
                    {
                        CiEqn += fvm::Sp(*poreSink_, C_[i]);
                    }
                    if (poreSrc_ && poreSrc_->set(i))
                    {
                        CiEqn -= (*poreSrc_)[i];
                    }
                }
                if (coupler && coupler->coupled(i))
                {
                    coupler->solve(i, CiEqn);
                }
                else
                {
                    CiEqn.solve();
                }
            }
        }

        // Line-search global (fraction-to-boundary): C >= 0 tanpa buang massa.
        //
        // AMBANG DERAU (2026-09-17, bug GR). Dulu setiap dC < -SMALL ikut
        // membatasi alpha. Pada sel dgn Ck <= 0 (C belum tersentuh front) derau
        // pembulatan solver linier memberi dC ~ -2e-15 -> amax = 0 -> alpha = 0
        // -> SELURUH medan C tak diperbarui, dan kriteria lama (alpha*|dC|)
        // menganggapnya konvergen. Grande 2006: alpha = 0 di 90 dari ~270
        // iterasi, 99,9% massa yg masuk dinding hilang, twoWayCoupling lalu
        // membalik aliran kanal.
        //
        // Aturan: sel ikut membatasi alpha HANYA bila Ck > epsC dan dC < -epsC,
        // epsC = newtonTol*max|C| (lantai derau yg sama dgn kriteria konvergen).
        // Sel dgn Ck <= epsC tak dapat dilindungi fraction-to-boundary (amax = 0
        // utk langkah negatif sekecil apa pun) dan nilainya sudah setara derau:
        // uji 2026-09-17 menunjukkan pembatasnya sel di z = 0,11 m (belum
        // tersentuh front) dgn Ck = -8e-13, dC = -6e-12. Langkah penuh di sel itu
        // tetap memenuhi persamaan C (tanpa clamp) -> konservatif; C negatif
        // yg tersisa berorde epsC.
        scalar alpha = 1.0;
        const scalar tau = 0.01;
        forAll(species_, i)
        {
            const scalarField& Ck = Ckeep[i].primitiveField();
            const scalarField& Cn = C_[i].primitiveField();
            const scalar epsC = max
            (
                newtonTol_*max(gMax(mag(Ck)()), gMax(mag(Cn)())),
                SMALL
            );
            forAll(Cn, c)
            {
                const scalar dC = Cn[c] - Ck[c];
                if (Ck[c] > epsC && dC < -epsC)
                {
                    const scalar amax = (tau - 1.0)*max(Ck[c], scalar(0))/dC;
                    if (amax < alpha) { alpha = amax; }
                }
            }
        }
        reduce(alpha, minOp<scalar>());
        alpha = max(alpha, scalar(0));

        // Terapkan langkah (teredam bila perlu) + ukur update relatif.
        // CATATAN: closure finite-rate KAKU (HSDM / LDF k-tinggi) BOCOR massa di
        // sini (A1 89%, B* 94-96%) krn exit Newton dini iters=2. Fix#1 (2026-07-01):
        // maxRel pakai |dCi| tak-teredam -> GAGAL (B2 gap ~85% tetap) + regresi
        // 50-iters -> DIREVERT. Lihat CATATAN_EVALUASI_ITERASI.txt BAGIAN 13.5.
        // Fix benar (belum dikerjakan): sink finite-rate IMPLISIT dlm C (tiru jalur
        // algebraic fvm::Sp(cross_ii/dt)+linearisasi dq/dC), bukan kriteria stop.
        scalar maxRel = 0;
        forAll(species_, i)
        {
            volScalarField dCi(C_[i] - Ckeep[i]);
            if (alpha < 1.0)
            {
                C_[i] == Ckeep[i] + alpha*dCi;
                if (coupler && coupler->coupled(i))
                {
                    coupler->damp(i, alpha);
                }
            }
            C_[i].correctBoundaryConditions();

            const scalar num = alpha*gMax(mag(dCi.primitiveField())());
            const scalar den = max(gMax(mag(C_[i].primitiveField())()), SMALL);
            maxRel = max(maxRel, num/den);
        }

        // Keluar HANYA bila langkah terakhir penuh (alpha = 1): langkah teredam
        // bukan solusi persamaan C, berapa pun kecil alpha*|dC|-nya. Kriteria
        // relatif sendiri TIDAK diperketat (Fix#1 2026-07-01 yg memakai |dC|
        // tak-teredam macet di 50 iterasi; CATATAN_EVALUASI_ITERASI 13.5 opsi a).
        if (iter >= 1 && alpha >= 1.0 && maxRel < newtonTol_)
        {
            ++iter;
            newtonConverged = true;
            break;
        }
    }

    // Loading akhir: q = q*(C_final) untuk species algebraic (konsisten m)
    {
        const scalarField& Tc = T.primitiveField();
        forAll(Tc, c)
        {
            forAll(species_, j) { Cvals[j] = C_[j].primitiveField()[c]; }
            evalPoint
            (
                isotherms_, *mixing_, M_, rhoS_, epsilonPtot_, fActive_,
                Tc[c], Cvals, qStarOut, effStorOut, STcoefOut, crossOut,
                pbuf, qcbuf, dqdT_pbuf, dqdp
            );
            forAll(species_, i)
            {
                qStar_[i].primitiveFieldRef()[c] = qStarOut[i];
            }
        }
        forAll(mesh_.boundary(), patchi)
        {
            const scalarField& Tp = T.boundaryField()[patchi];
            forAll(Tp, f)
            {
                forAll(species_, j)
                {
                    Cvals[j] = C_[j].boundaryField()[patchi][f];
                }
                evalPoint
                (
                    isotherms_, *mixing_, M_, rhoS_, epsilonPtot_, fActive_,
                    Tp[f], Cvals, qStarOut, effStorOut, STcoefOut, crossOut,
                    pbuf, qcbuf, dqdT_pbuf, dqdp
                );
                forAll(species_, i)
                {
                    qStar_[i].boundaryFieldRef()[patchi][f] = qStarOut[i];
                }
            }
        }
        forAll(species_, i)
        {
            if (algebraic[i])
            {
                q_[i] == qStar_[i];
            }
        }
    }

    Info<< "  solidSpecies[" << mesh_.name() << "]: Newton iters=" << iter
        << endl;

    if (!newtonConverged)
    {
        static label nWarnNewton = 0;
        if (nWarnNewton < 20)
        {
            ++nWarnNewton;
            WarningInFunction
                << "Newton padatan " << mesh_.name() << " TIDAK konvergen dalam "
                << maxIter << " iterasi pada t = " << mesh_.time().timeName()
                << " (langkah terakhir teredam atau update relatif > newtonTol "
                << newtonTol_ << "). Persamaan C tak terpenuhi -> neraca massa"
                << " padatan tak dijamin. (peringatan " << nWarnNewton << "/20)"
                << endl;
        }
    }

    // --- Source energi adsorpsi S_ads = sum_i dH_i rho_s dq_i/dt ---
    Sads_ == dimensionedScalar("0", dimEnergy/dimVolume/dimTime, 0);
    forAll(species_, i)
    {
        const dimensionedScalar dHDim("dH", dimEnergy/dimMoles, deltaH_[i]);
        Sads_ += dHDim*rhoSDim*(q_[i] - q_[i].oldTime())/dtDim;
    }

    Info<< "Adsorption[" << mesh_.name() << "]:";
    forAll(species_, i)
    {
        Info<< " " << species_[i]
            << " C=" << gMin(C_[i].primitiveField()) << "/"
            << gMax(C_[i].primitiveField())
            << " q=" << gMin(q_[i].primitiveField()) << "/"
            << gMax(q_[i].primitiveField()) << " mol/kg;";
    }
    Info<< " Sads=" << gMin(Sads_.primitiveField()) << "/"
        << gMax(Sads_.primitiveField()) << " W/m3" << endl;
}


// ************************************************************************* //
