    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimple.getOrDefault<int>("nOuterCorrectors", 1);

    // Kriteria konvergensi kopling species antarmuka; <= 0 (bawaan) = mati.
    // Lihat interfaceCouplingControl.H.
    const scalar interfaceTol =
        pimple.getOrDefault<scalar>("interfaceTolerance", -1);

    // Relaksasi kopling species antarmuka: none (bawaan, bit-identik) | aitken.
    // Lihat interfaceAitken.H.
    const word interfaceRelaxation =
        pimple.getOrDefault<word>("interfaceRelaxation", "none");
    if (interfaceRelaxation != "none" && interfaceRelaxation != "aitken")
    {
        FatalIOErrorInFunction(pimple)
            << "interfaceRelaxation '" << interfaceRelaxation
            << "' tidak dikenal. Pilihan: none, aitken." << exit(FatalIOError);
    }
    const bool interfaceAitken = (interfaceRelaxation == "aitken");

    // Urutan iterasi luar untuk kasus yang seluruh alirannya berasal dari gas terlepas (desorpsi vakum, V2)
    // dengan energi implisit: energi -> pEqn -> tekanan pori -> species terkopel. Bawaan no (urutan lama:
    // species padatan sebelum energi implisit dan pEqn). Lihat REGISTER PB2-STRIP-b.
    const bool speciesAfterPressure =
        energyImplicit && pimple.getOrDefault<bool>("speciesAfterPressure", false);
