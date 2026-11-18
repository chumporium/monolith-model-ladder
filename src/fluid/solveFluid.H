if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

if (frozenFlow)
{
    // Aliran beku (Paper B, 2026-10-03): U, p, rho dan phi tidak diselesaikan
    // (medan tunak dibaca dari waktu awal); species DAN energi tetap diselesaikan.
    // Bawaan chtMultiRegionFoam hanya energi -- untuk adsorber itu salah, karena
    // transport adsorbat justru yang dicari. Tak ada kasus lama yang memakai opsi ini.
    #include "YEqn.H"
    #include "EEqn.H"
    if (twoWayCoupling)
    {
        #include "fluidPoreSource.H"   // sumber gas pori untuk pEqn iterasi ini dan Y/E iterasi berikut
    }
}
else
{
    if (oCorr == 0)
    {
        #include "rhoEqn.H"
    }

    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"
    if (twoWayCoupling)
    {
        #include "fluidPoreSource.H"   // sumber gas pori untuk pEqn iterasi ini dan Y/E iterasi berikut
    }

    if (!coupled && !energyImplicit)
    {
        Info<< "\nSolving for fluid region " << fluidRegions[i].name() << endl;

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "pEqn.H"
        }

        turbulence.correct();

        rho = thermo.rho();

        #include "fluidFinalCorrectors.H"   // korektor species final + skalar modus-peluruhan
    }
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
