/*---------------------------------------------------------------------------*\
  updatePorePressure.H  (Paper B, 2026-10-02)

  Tekanan total gas pori di region solid, dipakai oleh Bosanquet Ds(p,T)
  di regionAdsorption (medan "p" milik solidThermo).

  PIMPLE { porePressure frozen | mapped; }   (bawaan frozen)

  frozen : p solid tidak disentuh (nilai dari 0/<solid>/p). Semua kasus
           Paper A (1 atm, Bosanquet mati) memakai ini, jadi hasilnya tetap.
  mapped : p solid = solusi Laplace  laplacian(p) = 0  dengan
             - fixedValue di patch antarmuka = p gas region fluida
               (lewat mappedPatchBase patch tersebut),
             - patch constraint (cyclic, empty, processor, ...) tetap,
             - zeroGradient di patch lain (dinding luar, simetri).
           Ini medan tekanan Darcy kuasi-tunak tanpa kapasitas: tekanan di
           pori mengikuti tekanan kanal di sekitarnya. Pada desorpsi vakum
           Ds -> D_Knudsen, sehingga Ds tidak peka terhadap rincian p.
           Aliran viskos (fluks Darcy) sendiri BELUM dimodelkan.
  darcy  : (Paper B, 2026-10-04, tahap V1) tekanan pori TRANSIEN dari neraca mol total:
             d/dt(eps_s p/(R T)) = div[ (1/(R T)) (eps_p Ds + p B/mu) grad p ] - rho_s SUM dq_i/dt
           p = p gas di antarmuka (fixedValue lewat mappedPatchBase); eps_s = pori total (simpanan),
           eps_p = pori transport, Ds = difusivitas pori (Bosanquet bila aktif) species pertama,
           B = `porePermeability` [m2], mu = `poreViscosity` [Pa s] (PIMPLE region padat).
           Batas gas murni = dusty-gas (Knudsen + viskos). Suku dq/dt memakai q iterasi terakhir.
           `poreAdvection yes` (tahap V2): fluks Darcy phiD = -(B/mu) snGrad(p)|Sf| mengadveksi species
           di dinding (regionAdsorption: + div(phiD C_i)), dan di antarmuka gas pori masuk/keluar kanal:
           sisi fluida menambah sumber massa dan sumber species (fluidAdsorptionSource.H).
           Entalpi gas yang menyeberang masuk EEqn fluida (fluidPoreSource.H); kapasitas panas gas pori
           di sisi padatan diabaikan (dicatat).
\*---------------------------------------------------------------------------*/
{
    const word porePressureMode
    (
        mesh.solutionDict().subOrEmptyDict("PIMPLE")
           .getOrDefault<word>("porePressure", "frozen")
    );

    if (porePressureMode == "mapped")
    {
        volScalarField& pS = thermo.p();

        wordList pPTypes(mesh.boundary().size());
        bool anyFixed = false;

        forAll(pPTypes, patchi)
        {
            const fvPatch& fp = mesh.boundary()[patchi];

            if (polyPatch::constraintType(fp.patch().type()))
            {
                pPTypes[patchi] = fp.patch().type();
            }
            else if
            (
                isA<mappedPatchBase>(fp.patch())
             && refCast<const mappedPatchBase>(fp.patch()).sampleRegion()
             != mesh.name()
            )
            {
                pPTypes[patchi] = fixedValueFvPatchScalarField::typeName;
                anyFixed = true;
            }
            else
            {
                pPTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
            }
        }

        if (!anyFixed)
        {
            FatalErrorInFunction
                << "porePressure mapped: region " << mesh.name()
                << " has no mapped interface patch to a fluid region"
                << exit(FatalError);
        }

        volScalarField pPore
        (
            IOobject
            (
                "pPore",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            pS,
            pPTypes
        );

        forAll(pPore.boundaryField(), patchi)
        {
            if (pPTypes[patchi] != fixedValueFvPatchScalarField::typeName)
            {
                continue;
            }

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(mesh.boundary()[patchi].patch());
            const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
            const label samplePatchi = mpp.samplePolyPatch().index();

            scalarField pNbr
            (
                nbrMesh.lookupObject<volScalarField>("p")
                    .boundaryField()[samplePatchi]
            );
            mpp.distribute(pNbr);

            pPore.boundaryFieldRef()[patchi] == pNbr;
        }

        fvScalarMatrix pPoreEqn(fvm::laplacian(pPore));

        if (mesh.solution().solversDict().found("pPore"))
        {
            pPoreEqn.solve(mesh.solverDict("pPore"));
        }
        else
        {
            dictionary sd;
            sd.add("solver", "PCG");
            sd.add("preconditioner", "DIC");
            sd.add("tolerance", 1e-10);
            sd.add("relTol", 0.0);
            pPoreEqn.solve(sd);
        }

        pS.primitiveFieldRef() = pPore.primitiveField();
        pS.correctBoundaryConditions();

        Info<< "porePressure - " << mesh.name() << ": min/max = "
            << gMin(pS.primitiveField()) << " / "
            << gMax(pS.primitiveField()) << " Pa" << endl;
    }
    else if (porePressureMode == "darcy")
    {
        const dictionary& ppDict = mesh.solutionDict().subDict("PIMPLE");
        const dimensionedScalar Bperm("B", dimArea, ppDict.get<scalar>("porePermeability"));
        const dimensionedScalar muG
        (
            "mu", dimPressure*dimTime, ppDict.getOrDefault<scalar>("poreViscosity", 1.8e-5)
        );
        const dimensionedScalar Rgas
        (
            "R", dimEnergy/dimMoles/dimTemperature, 8.314462618
        );
        volScalarField& pS = thermo.p();

        if (!pPoreFields.set(i))
        {
            wordList pPTypes(mesh.boundary().size());
            forAll(pPTypes, patchi)
            {
                const fvPatch& fp = mesh.boundary()[patchi];
                if (polyPatch::constraintType(fp.patch().type()))
                {
                    pPTypes[patchi] = fp.patch().type();
                }
                else if
                (
                    isA<mappedPatchBase>(fp.patch())
                 && refCast<const mappedPatchBase>(fp.patch()).sampleRegion()
                 != mesh.name()
                )
                {
                    pPTypes[patchi] = fixedValueFvPatchScalarField::typeName;
                }
                else
                {
                    pPTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
                }
            }
            pPoreFields.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "pPore", runTime.timeName(), mesh,
                        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
                    ),
                    pS,
                    pPTypes
                )
            );
            pPoreFields[i].oldTime();
        }
        volScalarField& pP = pPoreFields[i];

        forAll(pP.boundaryField(), patchi)
        {
            if (!isA<fixedValueFvPatchScalarField>(pP.boundaryField()[patchi]))
            {
                continue;
            }
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(mesh.boundary()[patchi].patch());
            const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
            const label samplePatchi = mpp.samplePolyPatch().index();
            scalarField pNbr
            (
                nbrMesh.lookupObject<volScalarField>("p").boundaryField()[samplePatchi]
            );
            mpp.distribute(pNbr);
            pP.boundaryFieldRef()[patchi] == pNbr;
        }

        const volScalarField rRT("rRT", 1.0/(Rgas*thermo.T()));
        const dimensionedScalar epsS("epsS", dimless, adsorptionI.epsilonPtot());
        const dimensionedScalar epsP("epsP", dimless, adsorptionI.epsilonP());
        const dimensionedScalar rhoSd("rhoS", dimDensity, adsorptionI.rhoS());

        // laju serapan total [mol/m3/s] (q iterasi terakhir; desorpsi -> negatif -> sumber gas)
        volScalarField dqdt
        (
            IOobject("dqdtPore", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
            mesh,
            dimensionedScalar(dimMoles/dimVolume/dimTime, Zero)
        );
        forAll(adsorptionI.species(), si)
        {
            dqdt += rhoSd*(adsorptionI.q(si) - adsorptionI.q(si).oldTime())/runTime.deltaT();
        }

        dictionary sd;
        sd.add("solver", "PCG");
        sd.add("preconditioner", "DIC");
        sd.add("tolerance", 1e-12);
        sd.add("relTol", 0.0);
        const bool adaSolver = mesh.solution().solversDict().found("pPore");

        // Penstabil kopling tekanan-serapan (REGISTER PB2-STRIP-b/c; bawaan MATI: tak berpengaruh pada 1 kPa, diuji untuk 0,1 kPa). Suku -rho_s dq/dt memakai q
        // iterasi terakhir, padahal q = q*(p_i, T) bergantung kuat pada p pori: pada vakum waktu relaksasi
        // tekanan pori (~1e-4 s) << dt, sehingga iterasi luar tekanan <-> species berosilasi dan jumlah yang
        // terlepas bergantung pada jumlah iterasi. Linearisasi:  rho_s dq_i/dt ~ c_i dp/dt,
        //   c_i = (effStor_i - eps_s) x_i / (R T),  x_i = fraksi mol species i di pori.
        // Ditambahkan  + c (p - p_iterasiLalu)/dt  secara implisit: nol saat konvergen (solusi tak berubah),
        // tetapi meredam osilasi (kapasitas serapan bekerja sebagai simpanan tekanan).
        const volScalarField pPrevIter("pPorePrev", pP);
        volScalarField cStab
        (
            IOobject("cStabPore", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
            mesh,
            dimensionedScalar(dimMoles/dimVolume/dimPressure, Zero)
        );
        if (ppDict.getOrDefault<bool>("poreSorptionStabilisation", false))
        {
            forAll(adsorptionI.species(), si)
            {
                const dimensionedScalar Mi("Mi", dimMass/dimMoles, adsorptionI.M(si));
                cStab += max(adsorptionI.effStor(si) - epsS, scalar(0))
                        *min(max(adsorptionI.C(si)/Mi/(rRT*pP), scalar(0)), scalar(1))*rRT;
            }
        }
        const dimensionedScalar rDt("rDt", dimless/dimTime, 1.0/runTime.deltaTValue());
        // poreStorageOldT (REGISTER 2026-10-16, diagnosis 0,1 kPa): fvm::ddt(epsS*rRT, p) memakai rRT SEKARANG untuk
        // keadaan lama (rRT dibangun ulang tiap langkah, tanpa oldTime), sehingga suku pemuaian termal gas pori
        // -p d(1/RT)/dt hilang: saat dinding memanas jumlah mol total pori tak turun, sedangkan C.CO2 dihitung benar ->
        // fraksi mol CO2 pori naik sebesar rasio suhu (x_CO2 > 1 pada 0,1 kPa). Bentuk benar (Euler):
        //   epsS (rRT p - rRT_lama p_lama)/dt.   Bawaan MATI (run lama bit-identik).
        const bool poreOldT = ppDict.getOrDefault<bool>("poreStorageOldT", false);
        // poreThermalConsistent (REGISTER 2026-10-16/17): di bawah gradien suhu, fluks mol total persamaan tekanan
        // (difusi (1/RT)_muka grad p, viskos (p/RT)_muka) tak sama dengan jumlah fluks species (difusi grad C,
        // adveksi phiDarcy p_muka/p_hulu x C_hulu tanpa 1/RT). Selisihnya = fluks inert semu yang mendorong
        // fraksi mol CO2 pori > 1 di sisi dingin (0,1 kPa). Bentuk konsisten (persis secara diskret pada mesh seragam):
        //   difusi total  = Dk_f snGrad(p/RT) = Dk_f [(1/RT)_f snGrad p + p_f snGrad(1/RT)]
        //   adveksi species phiC = phiDarcy (p/RT)_lin/(p/RT)_hulu;  buang di antarmuka phiDarcy_b (p/RT)_b/(p/RT)_sel.
        // Bawaan MATI (run lama bit-identik).
        const bool poreThermCons = ppDict.getOrDefault<bool>("poreThermalConsistent", false);
        // poreInterfaceSeries (REGISTER 2026-10-17; butuh poreThermalConsistent): difusi mol TOTAL di muka antarmuka memakai
        // konduktansi pori saja (Dk deltaCoeffs_padat), sedangkan difusi CO2 lewat kopling species memakai konduktansi SERI
        // pori + kanal G = Kf Ks/(Kf + Ks). Selisihnya muncul di kanal sebagai sumber inert mInert = W (Nt - SUM F_j/W_j)
        // yang dapat negatif saat gas hampir CO2 murni (N2 ditarik dari sel tanpa N2 -> Y_CO2 > 1; PB2-PORI-K (g)).
        // Di sini Dk muka antarmuka dikali w = Kf/(Kf + Ks) = valueFraction BC kopling C (padatan), sehingga difusi total
        // = G (c_tot,kanal - c_tot,pori) dan difusi inert tersirat = G (c_inert,kanal - c_inert,pori). Bawaan MATI.
        const bool poreIfSeries = poreThermCons && ppDict.getOrDefault<bool>("poreInterfaceSeries", false);
        // poreAdvectionBounded (REGISTER 2026-10-17, massSourceDiag D1/D2; butuh poreThermalConsistent): fluks adveksi species
        // = phiDarcy x c_tot,muka x x_j,hulu dengan x_j = C_j/(W_j c_tot). Tekanan pori (c_tot) dan C diselesaikan terpisah,
        // jadi SUM x_j dapat melewati 1 sedikit; adveksi lalu membawa CO2 > fluks mol total, sisa inert yang diserahkan ke kanal
        // (NtPore - SUM species) menjadi negatif, Y_CO2 kanal > 1, dan BC kopling mendorong C pori naik lagi (umpan balik;
        // D1: Adot - Mdot = +2,5 % fluks, Y_CO2 1,027). Di sini pembagi c_tot,hulu diganti max(c_tot, SUM C_j/W_j)_hulu:
        // fraksi mol yang terbawa adveksi <= 1, inert terbawa >= 0. Tak berpengaruh selama SUM x_j <= 1. Bawaan MATI.
        const bool poreAdvBound = poreThermCons && ppDict.getOrDefault<bool>("poreAdvectionBounded", false);
        auto cBagi = [&](const volScalarField& cTot)
        {
            volScalarField c("cBagiPore", cTot);
            if (poreAdvBound)
            {
                volScalarField cSp("cSpPore", 0.0*cTot);
                forAll(adsorptionI.species(), sj)
                {
                    cSp += adsorptionI.C(sj)/dimensionedScalar(dimMass/dimMoles, adsorptionI.M(sj));
                }
                c = max(cTot, cSp);
            }
            return c;
        };
        auto DkMuka = [&]()
        {
            surfaceScalarField Dkf(fvc::interpolate(epsP*adsorptionI.Ds(0)));
            if (poreIfSeries)
            {
                forAll(pP.boundaryField(), patchi)
                {
                    if (!isA<fixedValueFvPatchScalarField>(pP.boundaryField()[patchi])) continue;
                    const fvPatchScalarField& Cb = adsorptionI.C(0).boundaryField()[patchi];
                    if (!isA<mixedFvPatchScalarField>(Cb)) continue;
                    Dkf.boundaryFieldRef()[patchi] *= refCast<const mixedFvPatchScalarField>(Cb).valueFraction();
                }
            }
            return Dkf;
        };
        const volScalarField rRT0("rRT0", 1.0/(Rgas*thermo.T().oldTime()));

        for (label picard = 0; picard < 3; ++picard)
        {
            const volScalarField Kp
            (
                "KpPore",
                rRT*(epsP*adsorptionI.Ds(0) + pP*Bperm/muG)
            );
            tmp<fvScalarMatrix> tSimpan
            (
                poreOldT
              ? fvm::Sp(epsS*rRT*rDt, pP) - epsS*rRT0*rDt*pP.oldTime()
              : fvm::ddt(epsS*rRT, pP)
            );
            tmp<fvScalarMatrix> tLap;
            if (poreThermCons)
            {
                const surfaceScalarField Dkf(DkMuka());
                tLap = fvm::laplacian(fvc::interpolate(rRT*pP*Bperm/muG) + Dkf*fvc::interpolate(rRT), pP, "laplacian(KpPore,pPore)")
                     + fvc::div(Dkf*fvc::interpolate(pP)*fvc::snGrad(rRT)*mesh.magSf());
            }
            else
            {
                tLap = fvm::laplacian(Kp, pP, "laplacian(KpPore,pPore)");
            }
            fvScalarMatrix pPEqn
            (
                tSimpan
              + fvm::Sp(cStab*rDt, pP)
              - tLap
             ==
              - dqdt
              + cStab*rDt*pPrevIter
            );
            if (adaSolver) { pPEqn.solve(mesh.solverDict("pPore")); }
            else { pPEqn.solve(sd); }
        }

        pS.primitiveFieldRef() = pP.primitiveField();
        pS.correctBoundaryConditions();

        // Tahap V2: fluks Darcy membawa species (dan gas pori masuk/keluar kanal lewat antarmuka).
        if (ppDict.getOrDefault<bool>("poreAdvection", false))
        {
            if (speciesCouplerPtr)
            {
                // Kopling species implisit: Y fluida diselesaikan bersama C padatan SESUDAH pEqn, dengan sumber
                // adveksi dari pembacaan fluks yang sama dengan pEqn iterasi ini -> buku utang species lag 1.
                fadvLag = 1;
            }
            const surfaceScalarField snGp(fvc::snGrad(pP)*mesh.magSf());
            const volScalarField KpNow("KpPore", rRT*(epsP*adsorptionI.Ds(0) + pP*Bperm/muG));
            if (!phiDarcyFields.set(i))
            {
                phiDarcyFields.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject("phiDarcy", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                        -(Bperm/muG)*snGp
                    )
                );
                NtPoreFields.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject("NtPoreRaw", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                        poreThermCons
                      ? surfaceScalarField
                        (
                            -fvc::interpolate(rRT*pP*Bperm/muG)*snGp
                            - DkMuka()*(fvc::interpolate(rRT)*snGp + fvc::interpolate(pP)*fvc::snGrad(rRT)*mesh.magSf())
                        )
                      : surfaceScalarField(-fvc::interpolate(KpNow)*snGp)
                    )
                );
                phiCFields.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject("phiDarcyC", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                        phiDarcyFields[i]
                    )
                );
                phiCOutFields.set(i, new surfaceScalarField("phiDarcyCOut", 0.0*phiDarcyFields[i]));
                poreSinkFields.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject("poreSink", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                        mesh,
                        dimensionedScalar(dimless/dimTime, Zero)
                    )
                );
                adsorptionI.setPoreFlux(&phiCFields[i], &poreSinkFields[i]);
            }
            else
            {
                phiDarcyFields[i] = -(Bperm/muG)*snGp;
                if (poreThermCons)
                {
                    NtPoreFields[i] =
                        -fvc::interpolate(rRT*pP*Bperm/muG)*snGp
                        - DkMuka()*(fvc::interpolate(rRT)*snGp + fvc::interpolate(pP)*fvc::snGrad(rRT)*mesh.magSf());
                }
                else
                {
                    NtPoreFields[i] = -fvc::interpolate(KpNow)*snGp;
                }
            }
            if (!GporeFields.set(i))
            {
                GporeFields.set(i, new surfaceScalarField
                (
                    IOobject("GPore", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                    fvc::interpolate(KpNow)*mesh.magSf()*mesh.deltaCoeffs()
                ));
            }
            else
            {
                GporeFields[i] = fvc::interpolate(KpNow)*mesh.magSf()*mesh.deltaCoeffs();
            }

            // Fluks adveksi species konsisten-mol. Persamaan tekanan memakai p muka LINEAR (interpolate(Kp)),
            // species memakai C hulu; supaya SUM_j C_j/W_j x fluks = fluks mol viskos persamaan tekanan
            // (komposisi seragam tetap seragam):  phiC = phiDarcy p_linear/p_hulu.
            // Di antarmuka, gas yang KELUAR membawa komposisi SEL pori (hulu), bukan nilai batas BC kopling
            // (yang = komposisi kanal dan membuat selisih komposisi tumbuh sendiri, REGISTER PB0-V2-3e):
            // muka keluar di-nol-kan di phiC dan diganti pembuangan implisit  Sp = phiDarcy (p_b/p_sel)/V.
            {
                surfaceScalarField& phiC = phiCFields[i];
                const volScalarField cBagiSel(cBagi(rRT*pP));     // = c_tot bila opsi mati
                const surfaceScalarField& phiD = phiDarcyFields[i];
                const surfaceScalarField pLin(fvc::interpolate(pP));
                const surfaceScalarField pUp(upwind<scalar>(mesh, phiD).interpolate(pP));
                if (poreThermCons)
                {
                    const volScalarField cTot("cTotPore", rRT*pP);
                    phiC = phiD*fvc::interpolate(cTot)/upwind<scalar>(mesh, phiD).interpolate(cBagi(cTot));
                }
                else
                {
                    phiC = phiD*pLin/pUp;
                }
                phiCOutFields[i] = dimensionedScalar(phiD.dimensions(), Zero);
                scalarField& sink = poreSinkFields[i].primitiveFieldRef();
                sink = 0;
                const scalarField& Vc = mesh.V();
                forAll(pP.boundaryField(), patchi)
                {
                    if (!isA<fixedValueFvPatchScalarField>(pP.boundaryField()[patchi]))
                    {
                        continue;
                    }
                    const labelUList& fc = mesh.boundary()[patchi].faceCells();
                    const scalarField& pb = pP.boundaryField()[patchi];
                    const scalarField& phib = phiD.boundaryField()[patchi];
                    scalarField& phiCb = phiC.boundaryFieldRef()[patchi];
                    scalarField& outb = phiCOutFields[i].boundaryFieldRef()[patchi];
                    forAll(fc, f)
                    {
                        if (phib[f] > 0)
                        {
                            outb[f] = poreThermCons
                              ? phib[f]*pb[f]*rRT.boundaryField()[patchi][f]/cBagiSel[fc[f]]
                              : phib[f]*pb[f]/pP[fc[f]];
                            sink[fc[f]] += outb[f]/Vc[fc[f]];
                            phiCb[f] = 0;
                        }
                        else
                        {
                            phiCb[f] = 0;      // arah masuk ditangani sebagai sumber eksplisit di bawah
                        }
                    }
                }

                // Arah MASUK (kanal -> pori): hulu = sel kanal. Massa species j yang masuk per muka
                //   F_in = |phiDarcy_b| (rho Y_j)_sel kanal ; padatan: + F_in / V_sel ; kanal: - F_in (lewat FadvPore).
                const label nSp = adsorptionI.species().size();
                if (poreInSrcFields[i].size() != nSp)
                {
                    poreInSrcFields[i].setSize(nSp); poreInFluxFields[i].setSize(nSp);
                    forAll(adsorptionI.species(), si)
                    {
                        poreInSrcFields[i].set
                        (
                            si,
                            new volScalarField
                            (
                                IOobject("poreInSrc." + adsorptionI.species()[si], runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh, dimensionedScalar(dimDensity/dimTime, Zero)
                            )
                        );
                        poreInFluxFields[i].set
                        (
                            si,
                            new surfaceScalarField
                            (
                                IOobject("poreInFlux." + adsorptionI.species()[si], runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh, dimensionedScalar(dimMass/dimTime, Zero)
                            )
                        );
                    }
                    adsorptionI.setPoreSource(&poreInSrcFields[i]);
                }
                forAll(adsorptionI.species(), si)
                {
                    scalarField& src = poreInSrcFields[i][si].primitiveFieldRef();
                    src = 0;
                    poreInFluxFields[i][si] = dimensionedScalar(dimMass/dimTime, Zero);
                    forAll(pP.boundaryField(), patchi)
                    {
                        if (!isA<fixedValueFvPatchScalarField>(pP.boundaryField()[patchi]))
                        {
                            continue;
                        }
                        const mappedPatchBase& mpp =
                            refCast<const mappedPatchBase>(mesh.boundary()[patchi].patch());
                        const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
                        const label sP = mpp.samplePolyPatch().index();
                        if (!nbrMesh.foundObject<volScalarField>(adsorptionI.species()[si]))
                        {
                            continue;
                        }
                        scalarField cNbr
                        (
                            nbrMesh.lookupObject<volScalarField>("rho").boundaryField()[sP].patchInternalField()
                           *nbrMesh.lookupObject<volScalarField>(adsorptionI.species()[si]).boundaryField()[sP].patchInternalField()
                        );
                        mpp.distribute(cNbr);
                        const labelUList& fc = mesh.boundary()[patchi].faceCells();
                        const scalarField& phib = phiD.boundaryField()[patchi];
                        scalarField& fin = poreInFluxFields[i][si].boundaryFieldRef()[patchi];
                        forAll(fc, f)
                        {
                            if (phib[f] < 0)
                            {
                                const scalar F = -phib[f]*cNbr[f];        // kg/s masuk pori
                                src[fc[f]] += F/Vc[fc[f]];
                                fin[f] = -F;                               // tanda: keluar padatan positif
                            }
                        }
                    }
                }
            }
        }

        Info<< "porePressure(darcy) - " << mesh.name() << ": min/max = "
            << gMin(pS.primitiveField()) << " / "
            << gMax(pS.primitiveField()) << " Pa" << endl;
    }
}
