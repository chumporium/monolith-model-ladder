/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    adsorptionMultiRegionFoam

Description
    Transient solver for buoyant, turbulent fluid flow with conjugate heat
    transfer and adsorption between fluid and solid regions.

    Extends chtMultiRegionFoam with:
      - Intraparticle diffusion (Ci) in solid region
      - Local equilibrium Toth isotherm (qi*)
      - Adsorption heat source in solid energy equation
      - Concentration coupling at fluid-solid interface

    Ref: Stampi-Bombelli et al., Ind. Eng. Chem. Res. 2024

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedPatchBase.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "regionAdsorption.H"
#include "poreFluxCarry.H"
#include "mixedFvPatchFields.H"   // poreInterfaceSeries: valueFraction BC kopling C
#include "upwind.H"   // library adsorptionModels (file-scope, defines namespace)
#include "multivariateScheme.H"
#include "speciesImplicitCouplerFluid.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for fluid flow with conjugate heat transfer"
        " and adsorption between fluid and solid regions."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "adsorptionCouplingNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "createCoupledRegions.H"
    #include "createSpeciesCoupling.H"
    #include "createEnergyCoupling.H"
    #include "createAxialBins.H"

    // Faktor relaksasi Aitken, dibawa antar langkah waktu (interfaceAitken.H)
    scalar aitkenOmega = 1;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "adsorptionCouplingNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        // Store loading old-time BEFORE PIMPLE loop (q and closure state)
        forAll(solidRegions, i)
        {
            adsorption[i].storeOldTime();
        }

        #include "interfaceAitkenInit.H"

        // --- PIMPLE loop
        // kopKonvergen diisi interfaceCouplingControl.H; selama kontrol itu
        // mati ia tetap false dan loop identik dengan chtMultiRegionFoam.
        bool kopKonvergen = false;
        bool modeFinal = false;     // sekali FINAL, tetap FINAL sampai keluar
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1) || modeFinal;

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"   // U, p, Y_i (adsorbat), energi
            }

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                if (!speciesAfterPressure)
                {
                    #include "updatePorePressure.H"   // p pori (Bosanquet), bawaan frozen
                    #include "solveSolidSpecies.H"   // Ci equation (adsorption)
                    #include "updatePoreInterfaceFlux.H"
                }
                #include "solveSolid.H"           // energy equation
            }

            if (energyImplicit)
            {
                #include "solveEnergyImplicit.H"
            }

            if (speciesAfterPressure)
            {
                // species terkopel SESUDAH pEqn iterasi ini: memakai (rho, phi) dan sumber gas pori
                // dari pembacaan yang sama dengan pEqn (konsisten dengan kontinuitas)
                forAll(solidRegions, i)
                {
                    fvMesh& mesh = solidRegions[i];

                    #include "readSolidMultiRegionPIMPLEControls.H"
                    #include "setRegionSolidFields.H"
                    #include "updatePorePressure.H"
                    #include "solveSolidSpecies.H"
                    #include "updatePoreInterfaceFlux.H"
                }
            }


            if (coupled)
            {
                Info<< "\nSolving energy coupled regions " << endl;
                fvMatrixAssemblyPtr->solve();
                #include "correctThermos.H"

                forAll(fluidRegions, i)
                {
                    fvMesh& mesh = fluidRegions[i];

                    #include "readFluidMultiRegionPIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    if (!frozenFlow)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        for (int corr=0; corr<nCorr; corr++)
                        {
                            #include "pEqn.H"
                        }
                        turbulence.correct();
                    }

                    rho = thermo.rho();
                    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                        << max(thermo.T()).value() << endl;
                }

                fvMatrixAssemblyPtr->clear();
            }

            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        fvMesh& mesh = fluidRegions[i];

                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        #include "readFluidMultiRegionPIMPLEControls.H"
                        #include "setRegionFluidFields.H"
                        frozenFlow = true;
                        #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        fvMesh& mesh = solidRegions[i];

                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "setRegionSolidFields.H"
                        #include "updatePorePressure.H"
                        #include "solveSolidSpecies.H"
                        #include "updatePoreInterfaceFlux.H"
                        #include "solveSolid.H"
                    }

                    if (coupled)
                    {
                        Info<< "\nSolving energy coupled regions " << endl;
                        fvMatrixAssemblyPtr->solve();
                        #include "correctThermos.H"

                        forAll(fluidRegions, i)
                        {
                            #include "setRegionFluidFields.H"
                            rho = thermo.rho();
                        }

                        fvMatrixAssemblyPtr->clear();
                    }
                }
            }

            // Relaksasi Aitken nilai C antarmuka yang dilihat fluida (bila aktif)
            #include "interfaceAitken.H"

            // Diagnostik: sisa konservasi antarmuka PER ITERASI LUAR.
            // Mati kecuali `probeInterfaceLeak true` di controlDict.
            #include "interfaceLeakPerOuter.H"

            // Antarmuka dinilai SESUDAH tiap iterasi, termasuk iterasi FINAL:
            // yang dinilai harus sama dengan yang disimpan (bug 2026-09-13).
            #include "interfaceCouplingControl.H"

            if (interfaceTol <= 0)
            {
                if (finalIter)
                {
                    break;      // perilaku lama: FINAL hanya di iterasi terakhir
                }
            }
            else if (kopKonvergen)
            {
                if (finalIter)
                {
                    break;      // konservatif SESUDAH iterasi final
                }
                modeFinal = true;
            }
        }

        // Diagnostik P1: konservasi massa species di interface (tanpa fisika baru)
        #include "interfaceMassBalance.H"

        // Diagnostik energi: pasangan neraca massa (gerbang konservasi energi)
        #include "interfaceEnergyBalance.H"

        // Diagnostik neraca massa total gas vs adsorbat (Y_CO2 > 1); mati kecuali diminta
        #include "massSourceDiag.H"

        // Tahap 2.1: rata-rata per potongan aksial (aktif bila system/axialBins ada)
        #include "axialBins.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
