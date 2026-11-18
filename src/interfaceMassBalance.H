/*---------------------------------------------------------------------------*\
  interfaceMassBalance.H   (DIAGNOSTIK P1 — tanpa fisika baru)

  Mengukur konservasi massa species DI INTERFACE secara LANGSUNG, tiap langkah
  waktu, untuk menjawab: apakah BC speciesCoupledMixed benar-benar bocor, atau
  kebocoran berasal dari kopling interface yang belum konvergen?

  Fluks difusif CO2 (per interface, dijumlah seluruh region & face):
      Jf = sum ( -rho Dmol snGrad(Y) ) |Sf|     [kg/s]   sisi FLUID (gas keluar)
      Js = sum (  eps_p Ds snGrad(C) ) |Sf|     [kg/s]   sisi SOLID (masuk pori)

  Konstruksi mixed-BC bersifat konservatif PADA KONVERGENSI -> Jf == Js (mesin).
  leak = Jf - Js mengukur ketidak-konservasian sesaat (kg/s). Integralkan
  (trapesium, offline) -> total massa hilang; bandingkan dgn t50 dt-dependent.

  Inventori untuk NERACA TERTUTUP (sisi gas bebas yang DIABAIKAN skrip lama):
      channelGas = INT rho Y dV        (fluid)   [kg]
      poreGas    = eps_p INT C dV      (solid)   [kg]
      adsorbed   = rho_s M INT q dV    (solid)   [kg]

  Neraca tertutup (post): INT(inlet-outlet)dt
        = d(channelGas+poreGas+adsorbed) + kebocoran_numerik
  -> kalau Jf==Js TAPI neraca masih bocor => sumbernya BUKAN interface BC.

  Output: postProcessing/interfaceBalance.dat (master).
  Hanya membaca field; tidak mengubah solusi.
\*---------------------------------------------------------------------------*/
{
    scalar Jf = 0;          // fluks CO2 keluar gas   [kg/s]
    scalar Js = 0;          // fluks CO2 masuk solid  [kg/s]
    scalar channelGas = 0;  // gas bebas channel      [kg]
    scalar poreGas = 0;     // gas bebas pori         [kg]
    scalar adsorbed = 0;    // teradsorpsi (q)        [kg]

    // --- Sisi FLUID: fluks difusif yang meninggalkan fasa gas ---
    forAll(fluidRegions, ri)
    {
        const fvMesh& fmesh = fluidRegions[ri];
        const rhoReactionThermo& th = thermoFluid[ri];
        const auto& comp = th.composition();
        const PtrList<volScalarField>& Yf = comp.Y();
        const volScalarField& rhof = rhoFluid[ri];

        forAll(Yf, j)
        {
            const volScalarField& Yi = Yf[j];
            const scalar Dmol = DmolList[ri][j];

            bool isAdsorbate = false;
            forAll(fmesh.boundary(), patchi)
            {
                const fvPatchScalarField& Yp = Yi.boundaryField()[patchi];
                if (Yp.type() != "speciesCoupledMixed")
                {
                    continue;
                }
                isAdsorbate = true;

                const fvPatch& pp = fmesh.boundary()[patchi];
                const scalarField rhop(rhof.boundaryField()[patchi]);
                const scalarField snG(Yp.snGrad());     // (Yb-Yint)*deltaCoeffs
                Jf += gSum(-rhop*Dmol*snG*pp.magSf());
            }

            if (isAdsorbate)
            {
                const scalarField& Vf = fmesh.V();
                channelGas += gSum
                (
                    rhof.primitiveField()*Yi.primitiveField()*Vf
                );
            }
        }
    }

    // --- Sisi SOLID: fluks difusif yang masuk pori + inventori ---
    forAll(solidRegions, ri)
    {
        const fvMesh& smesh = solidRegions[ri];
        const regionAdsorption& ads = adsorption[ri];
        const scalar epsP = ads.epsilonP();       // transport -> Js
        const scalar epsS = ads.epsilonPtot();    // simpanan  -> poreGas
        const scalar rhoS = ads.rhoS();

        forAll(ads.species(), j)
        {
            const volScalarField& Cj = ads.C(j);
            const volScalarField& Dsj = ads.Ds(j);
            const volScalarField& qj = ads.q(j);
            const scalar Mj = ads.M(j);

            forAll(smesh.boundary(), patchi)
            {
                const fvPatchScalarField& Cp = Cj.boundaryField()[patchi];
                if (Cp.type() != "speciesCoupledMixed")
                {
                    continue;
                }
                const fvPatch& pp = smesh.boundary()[patchi];
                const scalarField Dsp(Dsj.boundaryField()[patchi]);
                const scalarField snG(Cp.snGrad());     // (Cb-Cint)*deltaCoeffs
                Js += gSum(epsP*Dsp*snG*pp.magSf());
            }

            const scalarField& Vs = smesh.V();
            poreGas  += epsS*gSum(Cj.primitiveField()*Vs);
            adsorbed += rhoS*Mj*gSum(qj.primitiveField()*Vs);
        }
    }

    // Laju massa adsorbat umpan -- skala gerbang di bawah (kolektif: semua rank).
    #include "interfaceFeedFlux.H"

    // --- NERACA GLOBAL (2026-09-17, bug GR) ---------------------------------
    // Jf == Js hanya menjamin antarmuka. Run Grande lolos gerbang itu sambil
    // kehilangan 99,9% massa di DALAM padatan (Newton palsu-konvergen) dan
    // disuplai difusi di inlet fixedValue yg tak terhitung `feed`. Yang
    // menangkapnya: fluks bersih masuk lewat SEMUA batas luar (konveksi +
    // difusi, rho Dmol sama dgn Jf) dibanding perubahan inventori total.
    // Jumlah lokal, reduce sekali (aturan MPI di interfaceFeedFlux.H).
    scalar netIn = 0;
    forAll(fluidRegions, ri)
    {
        const fvMesh& fmesh = fluidRegions[ri];
        const PtrList<volScalarField>& Yf = thermoFluid[ri].composition().Y();
        const surfaceScalarField& phiF = phiFluid[ri];
        const volScalarField& rhof = rhoFluid[ri];

        forAll(Yf, j)
        {
            const volScalarField& Yi = Yf[j];
            bool isAdsorbate = false;
            forAll(fmesh.boundary(), patchi)
            {
                if (Yi.boundaryField()[patchi].type() == "speciesCoupledMixed")
                {
                    isAdsorbate = true;
                }
            }
            if (!isAdsorbate) continue;

            const scalar Dmol = DmolList[ri][j];
            forAll(fmesh.boundary(), patchi)
            {
                if (fmesh.boundary()[patchi].coupled()) continue;
                const fvPatchScalarField& Yp = Yi.boundaryField()[patchi];
                if (Yp.type() == "speciesCoupledMixed") continue;

                const fvPatch& pp = fmesh.boundary()[patchi];
                const scalarField& phip = phiF.boundaryField()[patchi];
                const scalarField& rhop = rhof.boundaryField()[patchi];
                netIn += sum(-phip*Yp);                                  // konveksi
                netIn += sum(rhop*Dmol*Yp.snGrad()*pp.magSf());          // difusi masuk
            }
        }
    }
    reduce(netIn, sumOp<scalar>());

    if (Pstream::master())
    {
        static autoPtr<OFstream> osPtr;
        if (!osPtr)
        {
            const fileName ppDir(runTime.globalPath()/"postProcessing");
            mkDir(ppDir);
            osPtr.reset(new OFstream(ppDir/"interfaceBalance.dat"));
            osPtr().precision(12);   // uji jumlah-ulang axialBins (B21) butuh > 6 digit
            *osPtr
                << "# time  Jf[kg/s]  Js[kg/s]  leak=Jf-Js[kg/s]  "
                << "channelGas[kg]  poreGas[kg]  adsorbed[kg]  feed[kg/s]  "
                << "netIn[kg/s]  globalErr[-]" << endl;
        }

        // Kumulatif sejak panggilan pertama run INI (restart = acuan baru).
        // Euler implisit: fluks pd t^{n+1} x dt konsisten dgn inventori t^{n+1}.
        const scalar inv = channelGas + poreGas + adsorbed;
        static bool globalInit = false;
        static scalar inv0 = 0, cumIn = 0, cumAbs = 0;
        scalar globalErr = 0;
        if (!globalInit)
        {
            inv0 = inv; globalInit = true;
        }
        else
        {
            const scalar dtG = runTime.deltaTValue();
            cumIn  += netIn*dtG;
            cumAbs += mag(netIn)*dtG;
            const scalar skalaG = max(cumAbs, mag(inv - inv0));
            globalErr = (skalaG > VSMALL) ? ((inv - inv0) - cumIn)/skalaG : 0;
        }

        OFstream& os = *osPtr;
        os  << runTime.timeOutputValue() << ' '
            << Jf << ' ' << Js << ' ' << (Jf - Js) << ' '
            << channelGas << ' ' << poreGas << ' ' << adsorbed << ' '
            << fluksUmpan << ' ' << netIn << ' ' << globalErr << endl;

        // --- GERBANG KONSERVASI ---------------------------------------------
        // `leak` mutlak tak dapat ditafsir tanpa skala: 1e-11 kg/s besar pada
        // kasus 12 sel dan tak berarti pada monolith penuh. Yang dapat dinilai
        // adalah leak RELATIF. Sampai 2026-09-03 kolom leak hanya ditulis ke
        // berkas dan tak pernah dinilai siapa pun, sehingga kampanye 5,6%
        // berjalan berbulan-bulan kehilangan 29-56% massanya tanpa satu pun
        // peringatan. Uji satuan case/uji_antarmuka membuktikan sebabnya BUKAN
        // rumus BC (ia konvergen ke 6e-9) melainkan kopling antar-region yang
        // belum konvergen: pada 5,6% `nOuterCorrectors` 8 memberi sisa 1e-2,
        // sedangkan 32 memberi 7e-7. Angka 8 itu dipilih untuk kestabilan
        // Courant, bukan untuk konservasi -- dua kriteria berbeda yang tak
        // pernah diikat satu sama lain. Peringatan ini yang mengikatnya.
        // Skala sama dengan interfaceCouplingControl.H: max(umpan, |Jf|, |Js|).
        // Relatif thd |Jf| saja memberi peringatan palsu saat dinding jenuh
        // (T30_u008 t 287-293 s, Jf 2e-12, bocor 2e-6 dari total fluks).
        const scalar skala = max(fluksUmpan, max(mag(Jf), mag(Js)));
        const scalar sisa =
            (skala > VSMALL) ? mag(Jf - Js)/skala : scalar(0);

        Info<< "interfaceBalance: Jf=" << Jf << " Js=" << Js
            << " leak=" << (Jf - Js) << " kg/s"
            << " sisa=" << sisa
            << " (adsorbed=" << adsorbed << " kg)"
            << " netIn=" << netIn << " globalErr=" << globalErr << endl;

        // Gerbang neraca global. Ambang 1e-2 (bukan 1e-3 spt antarmuka): kesalahan
        // kontinuitas kompresibel dan clamp Y ikut masuk di sini; kasus GR yg
        // rusak duduk di ~1. Nilai sehat dicatat dari regresi sebelum diketatkan.
        static label nPeringatanG = 0;
        if (mag(globalErr) > 1e-2 && cumAbs > 1e-3*max(inv, VSMALL) && nPeringatanG < 20)
        {
            ++nPeringatanG;
            WarningInFunction
                << "Neraca massa GLOBAL tidak tutup: (dInventori - INT netIn)/skala = "
                << globalErr << " pada t = " << runTime.timeOutputValue() << nl
                << "    Antarmuka bisa konservatif sementara massa hilang di dalam"
                << " region. Angka run ini TIDAK boleh dikutip apa adanya."
                << " (peringatan " << nPeringatanG << "/20)" << endl;
        }

        // Ambang 1e-3: uji satuan menunjukkan kasus SEHAT duduk di 1e-7..1e-9,
        // sedangkan yang bocor mulai dari 1e-2. Satu dekade di antaranya cukup
        // longgar untuk derau diskretisasi, cukup ketat untuk menangkap cacat.
        static label nPeringatan = 0;
        if (sisa > 1e-3 && nPeringatan < 20)
        {
            ++nPeringatan;
            WarningInFunction
                << "Antarmuka TIDAK konservatif: |Jf-Js|/max = " << sisa
                << " (ambang 1e-3) pada t = " << runTime.timeOutputValue() << nl
                << "    Massa yang diserahkan fluida tidak diterima padatan."
                << " Angka run ini TIDAK boleh dikutip apa adanya." << nl
                << "    Sebab yang sudah terbukti: kopling antar-region belum"
                << " konvergen dalam satu langkah waktu." << nl
                << "    PERBAIKAN: naikkan nOuterCorrectors, atau turunkan"
                << " maxCo/deltaT. Lihat case/uji_antarmuka." << nl
                << "    (peringatan " << nPeringatan << "/20)" << endl;
        }
    }
}
