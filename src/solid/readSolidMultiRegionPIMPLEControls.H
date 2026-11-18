const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

int nNonOrthCorr =
    pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

// Sub-iterasi kopling q*<->T pada region solid (Bab 10.1).
// Regime encer (DAC trace): 1 cukup. Regime pekat/TVSA: 2-5.
const int nSolidCouple =
    pimple.getOrDefault<int>("nSolidCouple", 1);
