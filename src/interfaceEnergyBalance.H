/*---------------------------------------------------------------------------*\
  interfaceEnergyBalance.H   (DIAGNOSTIK — tanpa fisika baru, hanya membaca field)

  Pasangan energi dari interfaceMassBalance.H. Species punya gerbang konservasi
  sejak 2026-09-03; energi TIDAK punya, dan GR (Grande) berjalan dgn plato suhu
  13-15% dari kenaikan adiabatik (~86% panas hilang) tanpa satu peringatan pun
  (REGISTER_RUN.md, "KOREKSI GR termal", 2026-09-20). File ini menutup celah itu.

  Dua neraca, sengaja berbeda basis:

  (A) FISIS, basis suhu, acuan T0 tetap (rata-rata volume T fluida pd panggilan
      pertama run INI; restart = acuan baru, sama spt neraca massa):
        E_f = INT rho Cp (T-T0) dV          fluida  [J]
        E_s = INT betav rho Cp (T-T0) dV    padatan [J]
      dE_f/dt = Hconv + condF + qIfF
      dE_s/dt = Sads  + condS + qIfS
        Hconv = SUM (-phi) Cp_b (T_b-T0)          batas luar fluida    [W, masuk +]
        condF = SUM kappaEff_b snGrad(T) |Sf|     batas luar fluida    [W, masuk +]
        condS = SUM betav kappa_b snGrad(T) |Sf|  batas luar padatan   [W, masuk +]
        qIfF  = konduksi MASUK fluida lewat antarmuka                  [W]
        qIfS  = konduksi MASUK padatan lewat antarmuka                 [W]
        Sads  = INT S_ads dV  (= deltaH rho_s dq/dt)                   [W]
      leakE = qIfF + qIfS  (antarmuka konservatif -> 0)

  (B) SOLVER, basis `he` (variabel yang benar-benar diselesaikan EEqn.H):
        d/dt INT rho he dV = SUM(-phi he_b) + SUM alphaEff snGrad(he)|Sf|
                             (luar + antarmuka) + INT dpdt dV
      Ini menguji diskretisasi EEqn terhadap dirinya sendiri.

  Cara membaca:
    errHe ~ 0 tapi errEf besar -> EEqn konsisten dgn dirinya, tetapi RUMUSnya
        tidak menyimpan energi fisis (mis. sink massa Smass di kontinuitas tanpa
        pasangan entalpi di EEqn, atau pertukaran komposisi pd he tetap).
    errHe besar                -> diskretisasi/kopling EEqn sendiri bocor.
    errEs besar                -> persamaan energi padatan / sumber Sads.
    leakE besar                -> BC kopling suhu antarmuka tak konservatif.

  Diabaikan (dicatat, bukan dilupakan): energi kinetik K (U ~ mm/s), fasa
  teradsorpsi sensibel (n_ads cp (T-T0), < 1% pd GR), radiasi (tak aktif).
  Padatan anisotropik: kappa isotropik dipakai -> kolom condS/qIfS tak berlaku.

  Aturan MPI (interfaceFeedFlux.H): jumlah LOKAL di semua loop, satu reduce.
  Output: postProcessing/energyBalance.dat (master).
\*---------------------------------------------------------------------------*/
{
    // Indeks akumulator lokal (satu reduce di akhir).
    enum
    {
        iFRhoCpT, iFRhoCp, iFTV, iFV,           // inventori fluida (basis T)
        iConvCpT, iConvCp,                      // konveksi batas luar fluida
        iCondF, iQIfF,                          // konduksi luar / antarmuka fluida
        iFRhoHe, iConvHe, iCondHe, iIfHe, iDpdt,// neraca he solver
        iSmass, iSmassHe,                       // sink massa & entalpinya (info)
        iSRhoCpT, iSRhoCp,                      // inventori padatan
        iSads, iCondS, iQIfS,                   // sumber / konduksi padatan
        iLoss,                                  // kehilangan panas dinding (wallHeatLoss)
        nE
    };
    List<scalar> e(nE, Zero);

    // Patch antarmuka suhu dikenali dari tipe BC-nya (keluarga *CoupledMixed).
    auto isIfT = [](const fvPatchScalarField& Tp)
    {
        return Tp.type().find("CoupledMixed") != std::string::npos;
    };

    bool heIsE = false;

    // ------------------------------------------------------------- FLUIDA
    forAll(fluidRegions, ri)
    {
        const fvMesh& fmesh = fluidRegions[ri];
        const rhoReactionThermo& th = thermoFluid[ri];
        const volScalarField& T = th.T();
        const volScalarField& he = th.he();
        const volScalarField& rhof = rhoFluid[ri];
        const surfaceScalarField& phiF = phiFluid[ri];
        const tmp<volScalarField> tCp(th.Cp());
        const volScalarField& Cp = tCp();
        const tmp<volScalarField> tKap(turbulenceFluid[ri].kappaEff());
        const tmp<volScalarField> tAlp(turbulenceFluid[ri].alphaEff());
        const volScalarField& kap = tKap();
        const volScalarField& alp = tAlp();
        const scalarField& V = fmesh.V();

        heIsE = heIsE || (he.name() == "e");

        e[iFRhoCpT] += sum(rhof.primitiveField()*Cp.primitiveField()*T.primitiveField()*V);
        e[iFRhoCp]  += sum(rhof.primitiveField()*Cp.primitiveField()*V);
        e[iFTV]     += sum(T.primitiveField()*V);
        e[iFV]      += sum(V);
        e[iFRhoHe]  += sum(rhof.primitiveField()*he.primitiveField()*V);
        e[iDpdt]    += sum(dpdtFluid[ri].primitiveField()*V);
        e[iSmass]   += sum(SmassFluid[ri].primitiveField()*V);
        e[iSmassHe] += sum(SmassFluid[ri].primitiveField()*he.primitiveField()*V);

        forAll(fmesh.boundary(), patchi)
        {
            const fvPatch& pp = fmesh.boundary()[patchi];
            if (pp.coupled()) continue;                 // processor/cyclic: internal

            const fvPatchScalarField& Tp = T.boundaryField()[patchi];
            const fvPatchScalarField& hp = he.boundaryField()[patchi];
            const scalarField magSf(pp.magSf());
            const scalarField kapp(kap.boundaryField()[patchi]);
            const scalarField alpp(alp.boundaryField()[patchi]);
            const scalarField qT(kapp*Tp.snGrad()*magSf);     // masuk domain +
            const scalarField qH(alpp*hp.snGrad()*magSf);

            if (isIfT(Tp))
            {
                e[iQIfF] += sum(qT);
                e[iIfHe] += sum(qH);
            }
            else
            {
                const scalarField& phip = phiF.boundaryField()[patchi];
                const scalarField& Cpp = Cp.boundaryField()[patchi];
                e[iConvCpT] += sum(-phip*Cpp*Tp);
                e[iConvCp]  += sum(-phip*Cpp);
                e[iCondF]   += sum(qT);
                e[iConvHe]  += sum(-phip*hp);
                e[iCondHe]  += sum(qH);
            }
        }
    }

    // ------------------------------------------------------------ PADATAN
    forAll(solidRegions, ri)
    {
        const fvMesh& smesh = solidRegions[ri];
        const solidThermo& th = thermos[ri];
        const volScalarField& T = th.T();
        const volScalarField& bv = betavSolid[ri];
        const tmp<volScalarField> tRho(th.rho());
        const tmp<volScalarField> tCp(th.Cp());
        const tmp<volScalarField> tKap(th.kappa());
        const scalarField& V = smesh.V();
        const scalarField w(bv.primitiveField()*tRho().primitiveField()*tCp().primitiveField()*V);

        e[iSRhoCpT] += sum(w*T.primitiveField());
        e[iSRhoCp]  += sum(w);
        e[iSads]    += sum(adsorption[ri].Sads().primitiveField()*V);
        if (hvLoss[ri] > 0)
        {
            e[iLoss] += sum(-hvLoss[ri]*(T.primitiveField() - TambLoss[ri])*V);
        }

        forAll(smesh.boundary(), patchi)
        {
            const fvPatch& pp = smesh.boundary()[patchi];
            if (pp.coupled()) continue;

            const fvPatchScalarField& Tp = T.boundaryField()[patchi];
            const scalarField q
            (
                bv.boundaryField()[patchi]*tKap().boundaryField()[patchi]
               *Tp.snGrad()*pp.magSf()
            );
            if (isIfT(Tp)) e[iQIfS] += sum(q);
            else           e[iCondS] += sum(q);
        }
    }

    Pstream::listCombineReduce(e, plusEqOp<scalar>());

    if (Pstream::master())
    {
        static autoPtr<OFstream> osE;
        static bool initE = false;
        static scalar T0 = 0, Ef0 = 0, Es0 = 0, He0 = 0;
        static scalar cumF = 0, cumS = 0, cumHe = 0, cumAbs = 0, cumSads = 0;
        static scalar cumConv = 0, cumCondF = 0, cumCondS = 0, cumLeak = 0, cumLoss = 0;

        if (!initE)
        {
            T0 = e[iFTV]/max(e[iFV], VSMALL);
        }

        // Inventori & konveksi relatif thd T0 (dirakit setelah reduce).
        const scalar Ef = e[iFRhoCpT] - T0*e[iFRhoCp];
        const scalar Es = e[iSRhoCpT] - T0*e[iSRhoCp];
        const scalar Hconv = e[iConvCpT] - T0*e[iConvCp];
        const scalar leakE = e[iQIfF] + e[iQIfS];

        scalar errEf = 0, errEs = 0, errHe = 0, errE = 0;
        if (!initE)
        {
            Ef0 = Ef; Es0 = Es; He0 = e[iFRhoHe];
            initE = true;

            const fileName ppDir(runTime.globalPath()/"postProcessing");
            mkDir(ppDir);
            osE.reset(new OFstream(ppDir/"energyBalance.dat"));
            *osE
                << "# T0 = " << T0 << " K (rata2 volume T fluida saat mulai)" << nl
                << "# satuan: E,cum,err* [J]; Sads..leakE [W]; errE [-]" << nl
                << "# time  Ef  Es  Sads  Hconv  condF  condS  qIfF  qIfS  leakE"
                << "  cumSads  cumConv  cumCondF  cumCondS  cumLeak"
                << "  errEf  errEs  errHe  errE  Smass[kg/s]  SmassHe[W]  Qloss[W]  cumLoss" << endl;
        }
        else
        {
            const scalar dt = runTime.deltaTValue();
            cumF    += (Hconv + e[iCondF] + e[iQIfF])*dt;
            cumS    += (e[iSads] + e[iCondS] + e[iQIfS] + e[iLoss])*dt;
            cumHe   += (e[iConvHe] + e[iCondHe] + e[iIfHe]
                      + (heIsE ? scalar(0) : e[iDpdt]))*dt;
            cumSads += e[iSads]*dt;
            cumConv += Hconv*dt;
            cumCondF += e[iCondF]*dt;
            cumCondS += e[iCondS]*dt;
            cumLeak += leakE*dt;
            cumLoss += e[iLoss]*dt;
            cumAbs  += (mag(e[iSads]) + mag(Hconv) + mag(e[iCondF])
                      + mag(e[iCondS]) + mag(e[iLoss]))*dt;

            errEf = (Ef - Ef0) - cumF;
            errEs = (Es - Es0) - cumS;
            errHe = (e[iFRhoHe] - He0) - cumHe;

            // Global fisis: antarmuka saling hapus -> hanya batas luar + Sads.
            const scalar dE = (Ef - Ef0) + (Es - Es0);
            const scalar cumIn = cumSads + cumConv + cumCondF + cumCondS + cumLoss;
            const scalar skala = max(cumAbs, mag(dE));
            errE = (skala > VSMALL) ? (dE - cumIn)/skala : scalar(0);
        }

        *osE
            << runTime.timeOutputValue() << ' '
            << Ef << ' ' << Es << ' ' << e[iSads] << ' ' << Hconv << ' '
            << e[iCondF] << ' ' << e[iCondS] << ' ' << e[iQIfF] << ' '
            << e[iQIfS] << ' ' << leakE << ' '
            << cumSads << ' ' << cumConv << ' ' << cumCondF << ' '
            << cumCondS << ' ' << cumLeak << ' '
            << errEf << ' ' << errEs << ' ' << errHe << ' ' << errE << ' '
            << e[iSmass] << ' ' << e[iSmassHe] << ' ' << e[iLoss] << ' ' << cumLoss << endl;

        Info<< "energyBalance: Sads=" << e[iSads] << " W  dE="
            << (Ef - Ef0) + (Es - Es0) << " J  errE=" << errE
            << "  (errEf=" << errEf << " errEs=" << errEs
            << " errHe=" << errHe << " J; leakE=" << leakE << " W)" << endl;

        // Gerbang. Ambang 1e-2 sama dgn neraca massa global; kasus GR yg rusak
        // diduga duduk di ~0,8. Aktif hanya setelah ada energi berarti mengalir.
        static label nPeringatanE = 0;
        if
        (
            mag(errE) > 1e-2
         && cumAbs > 1e-6*max(mag(Ef) + mag(Es), VSMALL)
         && nPeringatanE < 20
        )
        {
            ++nPeringatanE;
            WarningInFunction
                << "Neraca ENERGI global tidak tutup: (dE - INT masuk)/skala = "
                << errE << " pada t = " << runTime.timeOutputValue() << nl
                << "    errEf=" << errEf << " errEs=" << errEs << " errHe="
                << errHe << " J. Lihat kepala interfaceEnergyBalance.H untuk"
                << " cara membaca kolom mana yg bocor." << nl
                << "    Angka SUHU run ini TIDAK boleh dikutip apa adanya."
                << " (peringatan " << nPeringatanE << "/20)" << endl;
        }
    }
}
