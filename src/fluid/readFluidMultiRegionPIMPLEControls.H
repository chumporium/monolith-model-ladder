    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    const int nCorr =
        pimple.getOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimple.getOrDefault("momentumPredictor", true);

    pimple.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    // Kopling massa dua-arah (Bab 5.1/8): sink massa adsorpsi di interface
    // masuk ke kontinuitas. ~0 otomatis pada regime trace (DAC 400 ppm).
    const bool twoWayCoupling =
        pimple.getOrDefault("twoWayCoupling", true);
