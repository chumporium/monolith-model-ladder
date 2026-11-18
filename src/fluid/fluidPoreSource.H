// fluidPoreSource.H (Paper B, tahap V2) -- di-include di solveFluid.H SESUDAH EEqn, sebelum pEqn.
// Urutan ini disengaja (REGISTER PB0-V2-3d): YEqn/EEqn iterasi luar k memakai (rho, phi) hasil pEqn k-1,
// jadi sumber species/entalpi-nya harus berasal dari pembacaan fluks yang SAMA dengan sumber massa pEqn k-1.
// Sumber massa dipakai pEqn pada iterasi yang sama; sumber species dan entalpi pada iterasi berikutnya.
{
    // ------------------------------------------------------------------
    // Tahap V2 (Paper B, 2026-10-04): gas pori yang menyeberang antarmuka.
    // Region padat dgn `porePressure darcy` + `poreAdvection yes` mendaftarkan
    //   NtPore        [mol/s] fluks mol total (Knudsen/difusif + viskos), keluar padatan = masuk fluida
    //   FadvPore.<sp> [kg/s]  adveksi Darcy species terlacak (C muka x phiDarcy)
    //   (keduanya diterbitkan updatePoreInterfaceFlux.H lewat poreFluxCarry: konservatif thd iterasi luar)
    // Masuk fluida per muka:
    //   F_adv,j  = FadvPore.j                            [kg/s] species terlacak, adveksi
    //   F_diff,j = rho D (Y_muka - Y_sel) deltaCoeffs A [kg/s] difusif (sudah ada di Smass di atas)
    //   inert    = W_inert [ NtPore - SUM_j (F_adv,j + F_diff,j)/W_j ]
    // Tambahan massa (di luar bagian difusif) = SUM_j F_adv,j + inert  -> Spore [kg/m3/s] (Smass -= Spore)
    // Sumber species j di YEqn: + F_adv,j / V.
    // Entalpi: SentPore = SUM F h_k(T asal)/V masuk EEqn (`poreEnthalpy`, bawaan yes); T asal = muka
    // dinding bila gas masuk fluida, sel bila keluar. Tanpa suku ini gas datang dgn he = 0 (bergantung acuan).
    // ------------------------------------------------------------------
    {
        if (SadvFluid[i].size() != Y.size())
        {
            SadvFluid[i].setSize(Y.size());
        }
        forAll(Y, j)
        {
            if (SadvFluid[i].set(j))
            {
                SadvFluid[i][j] == dimensionedScalar(dimDensity/dimTime, Zero);
            }
        }
        const scalar Winert =
            (inertIndex >= 0 ? thermo.composition().W(inertIndex)*1e-3 : 0.0280134);

        if (!SporeFluid.set(i))
        {
            SporeFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject("Spore", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                    mesh,
                    dimensionedScalar(dimDensity/dimTime, Zero)
                )
            );
            SentPoreFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject("SentPore", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                    mesh,
                    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
                )
            );
        }
        SporeFluid[i] == dimensionedScalar(dimDensity/dimTime, Zero);
        poreDiagFluid[i] = Zero;
        // poreSourceImplicit (RANCANGAN_PEQN_IMPLISIT.md): dS/dp = W_campuran x G_muka / V; p* = p_rgh saat ini.
        const bool poreImpl =
            mesh.solutionDict().subDict("PIMPLE").getOrDefault<bool>("poreSourceImplicit", false);
        if (poreImpl && !GporeFluid.set(i))
        {
            GporeFluid.set(i, new volScalarField
            (
                IOobject("Gpore", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedScalar(dimDensity/dimTime/dimPressure, Zero)
            ));
            pStarPoreFluid.set(i, new volScalarField("p_rghStarPore", p_rgh));
        }
        if (poreImpl)
        {
            GporeFluid[i] == dimensionedScalar(dimDensity/dimTime/dimPressure, Zero);
            pStarPoreFluid[i] = p_rgh;
        }
        SentPoreFluid[i] == dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero);
        scalarField& sCells = SporeFluid[i].primitiveFieldRef();
        scalarField& hCells = SentPoreFluid[i].primitiveFieldRef();
        const scalarField& V = mesh.V();
        const volScalarField& Tg = thermo.T();
        const volScalarField& pg = thermo.p();
        const bool poreEnthalpy =
            mesh.solutionDict().subDict("PIMPLE").getOrDefault<bool>("poreEnthalpy", true);
        // entalpi sensibel species k pada suhu asal gas: muka dinding bila masuk fluida, sel bila keluar
        auto hAsal = [&](const label k, const scalar F, const label patchi, const label f, const label c)
        {
            return F > 0
              ? thermo.composition().HE(k, pg.boundaryField()[patchi][f], Tg.boundaryField()[patchi][f])
              : thermo.composition().HE(k, pg[c], Tg[c]);
        };

        const bool inertLimit =
            mesh.solutionDict().subDict("PIMPLE").getOrDefault<bool>("inertAvailabilityLimit", false);
        // poreMassFromSpecies (REGISTER 2026-10-18, PB2-N2; bawaan MATI): bila SEMUA komponen gas nyata diangkut sebagai species
        // terkopel (N2 ikut sebagai species padatan berkapasitas nol, inertSpecie = penutup Y = 0), sumber massa kanal = jumlah fluks
        // species saja; sisa W_inert (NtPore - SUM species) TIDAK diserahkan ke kanal. Ketidakcocokan tekanan pori vs jumlah species
        // tetap di dalam pori (dilaporkan massSourceDiag kolom mInertMentah), tidak lagi menggeser komposisi kanal.
        const bool massFromSpecies =
            mesh.solutionDict().subDict("PIMPLE").getOrDefault<bool>("poreMassFromSpecies", false);
        scalarField nMukaIf(mesh.nCells(), Zero);      // jumlah muka antarmuka pori per sel (pembagi pembatas)
        if (inertLimit)
        {
            forAll(mesh.boundary(), patchi)
            {
                if (!isA<mappedPatchBase>(mesh.boundary()[patchi].patch())) continue;
                for (const label c : mesh.boundary()[patchi].faceCells()) nMukaIf[c] += 1;
            }
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& pp = mesh.boundary()[patchi];
            if (!isA<mappedPatchBase>(pp.patch()))
            {
                continue;
            }
            const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp.patch());
            const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
            if (!nbrMesh.foundObject<surfaceScalarField>("NtPore"))
            {
                continue;
            }
            const label sP = mpp.samplePolyPatch().index();
            const surfaceScalarField& Nt =
                nbrMesh.lookupObject<surfaceScalarField>("NtPore");

            scalarField Ft(Nt.boundaryField()[sP]);
            mpp.distribute(Ft);                               // mol/s masuk fluida

            const labelUList& fc = pp.faceCells();
            if (poreImpl && nbrMesh.foundObject<surfaceScalarField>("GPore"))
            {
                scalarField Gf(nbrMesh.lookupObject<surfaceScalarField>("GPore").boundaryField()[sP]);
                mpp.distribute(Gf);                          // mol/s/Pa per muka
                const volScalarField Wmix(thermo.W());       // kg/kmol
                scalarField& Gc = GporeFluid[i].primitiveFieldRef();
                forAll(fc, f)
                {
                    Gc[fc[f]] += 1e-3*Wmix[fc[f]]*Gf[f]/V[fc[f]];
                }
            }
            scalarField molSpecies(pp.size(), Zero);
            scalarField massIn(pp.size(), Zero);

            forAll(Y, j)
            {
                if (j == inertIndex || !composition.active(j))
                {
                    continue;
                }
                const word Fname("FadvPore." + Y[j].name());
                if (!nbrMesh.foundObject<surfaceScalarField>(Fname))
                {
                    continue;
                }
                const scalar Wj = thermo.composition().W(j)*1e-3;

                scalarField Fadv
                (
                    nbrMesh.lookupObject<surfaceScalarField>(Fname).boundaryField()[sP]
                );
                mpp.distribute(Fadv);                         // kg/s masuk fluida

                const fvPatchScalarField& Yp = Y[j].boundaryField()[patchi];
                scalarField Fdiff(pp.size(), Zero);
                if (Yp.type() == "speciesCoupledMixed")
                {
                    Fdiff =
                        rho.boundaryField()[patchi]*DmolList[i][j]
                       *(Yp - Yp.patchInternalField())*pp.deltaCoeffs()*pp.magSf();
                }
                molSpecies += (Fadv + Fdiff)/Wj;
                massIn += Fadv;
                poreDiagFluid[i][1] += sum(Fadv);
                poreDiagFluid[i][2] += sum(Fdiff);

                if (!SadvFluid[i].set(j))
                {
                    SadvFluid[i].set
                    (
                        j,
                        new volScalarField
                        (
                            IOobject
                            (
                                "Sadv." + Y[j].name(), runTime.timeName(), mesh,
                                IOobject::NO_READ, IOobject::NO_WRITE
                            ),
                            mesh,
                            dimensionedScalar(dimDensity/dimTime, Zero)
                        )
                    );
                }
                scalarField& sj = SadvFluid[i][j].primitiveFieldRef();
                forAll(fc, f)
                {
                    sj[fc[f]] += Fadv[f]/V[fc[f]];
                    if (poreEnthalpy)
                    {
                        hCells[fc[f]] += Fadv[f]*hAsal(j, Fadv[f], patchi, f, fc[f])/V[fc[f]];
                    }
                }
            }

            scalarField mInert(Winert*(Ft - molSpecies));
            poreDiagFluid[i][5] += sum(mInert);
            // inertAvailabilityLimit (REGISTER 2026-10-17, massSourceDiag D3/D4; bawaan MATI): sisa inert negatif = N2 mengalir dari kanal
            // ke pori. Itu fisis selama sel kanal punya N2, tetapi rumus sisa tetap "menarik" N2 setelah N2 kanal habis (N2 dipotong di 0),
            // sehingga massa total yang masuk < CO2 yang masuk dan Y_CO2 > 1 (D4: -3,3 % fluks pada 50 s). Pembatas positivitas: N2 yang
            // ditarik per langkah <= isi N2 sel kanal (dibagi rata antar muka antarmuka sel itu).
            if (massFromSpecies)
            {
                mInert = 0;
            }
            else if (inertLimit && inertIndex >= 0)
            {
                const scalarField& Yin = Y[inertIndex].primitiveField();
                const scalar dtNow = runTime.deltaTValue();
                forAll(fc, f)
                {
                    const label c = fc[f];
                    const scalar batas = -rho[c]*max(Yin[c], scalar(0))*V[c]/(dtNow*max(nMukaIf[c], scalar(1)));
                    mInert[f] = max(mInert[f], batas);
                }
            }
            massIn += mInert;
            poreDiagFluid[i][0] += sum(Ft);
            poreDiagFluid[i][3] += sum(molSpecies);
            poreDiagFluid[i][4] += sum(mInert);
            forAll(fc, f)
            {
                sCells[fc[f]] += massIn[f]/V[fc[f]];
                if (poreEnthalpy && inertIndex >= 0)
                {
                    hCells[fc[f]] += mInert[f]*hAsal(inertIndex, mInert[f], patchi, f, fc[f])/V[fc[f]];
                }
            }
        }
    }
    if (SporeFluid.set(i))
    {
        SmassFluid[i] -= SporeFluid[i];
    }
}
