/*---------------------------------------------------------------------------*\
  solveEnergyImplicit.H  -- energyInterfaceCoupling implicit, tiap iterasi luar

  Dipanggil sesudah semua region fluida (EEqn.H) dan padatan (solveSolid.H)
  menyimpan persamaan energinya ke energyEqns (lihat createEnergyCoupling.H).
  Rakit -> selesaikan -> salin hCpl ke he -> perbarui thermo -> baru pEqn
  (urutan sama dgn jalur `coupled` bawaan: pEqn memakai T yang baru).
\*---------------------------------------------------------------------------*/
{
    fvMatrix<scalar> asmb(hCplFluid[0], dimEnergy/dimTime);
    forAll(energyEqns, k)
    {
        asmb.addFvMatrix(energyEqns[k]);
    }

    Info<< "\nSolving energy implicit (fluid + solid)" << endl;
    asmb.solve
    (
        fluidRegions[0].solverDict(finalIter ? word("hFinal") : word("h"))
    );
    energyEqns.clear();

    forAll(fluidRegions, i)
    {
        volScalarField& he = thermoFluid[i].he();
        he.primitiveFieldRef() = hCplFluid[i].primitiveField();
        he.correctBoundaryConditions();
    }
    forAll(solidRegions, i)
    {
        volScalarField& h = thermos[i].he();
        h.primitiveFieldRef() = hCplSolid[i].primitiveField();
        h.correctBoundaryConditions();
    }

    #include "correctThermos.H"

    forAll(fluidRegions, i)
    {
        fvMesh& mesh = fluidRegions[i];

        #include "readFluidMultiRegionPIMPLEControls.H"
        #include "setRegionFluidFields.H"
        if (!frozenFlow)
        {
            Info<< "\nSolving for fluid region " << fluidRegions[i].name() << endl;
            for (int corr=0; corr<nCorr; corr++)
            {
                #include "pEqn.H"
            }
            turbulence.correct();
        }

        rho = thermo.rho();
        if (!frozenFlow)
        {
            #include "fluidFinalCorrectors.H"
        }
        Info<< "Min/max T (fluid):" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
    forAll(solidRegions, i)
    {
        Info<< "Min/max T (solid):" << min(thermos[i].T()).value() << ' '
            << max(thermos[i].T()).value() << endl;
    }
}
