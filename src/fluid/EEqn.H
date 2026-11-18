{
    volScalarField& he = thermo.he();

    #include "energySpeciesEnthalpy.H"

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence.alphaEff(), he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
      + Qdot
      + SentF
      + fvOptions(rho, he)
    );

    if (SentPoreFluid.set(i))
    {
        EEqn -= SentPoreFluid[i];      // + entalpi gas pori yang menyeberang (tahap V2)
    }

    // Koreksi kontinuitas (energyContinuityCorrection, createEnergyCoupling.H): EEqn -= he x residu kontinuitas.
    // Smass total yang dipakai pEqn = bagian difusif (SmassFluid) - gas pori (SporeFluid, iterasi sebelumnya).
    if (energyContinuityCorrection)
    {
        volScalarField contErr("contErrEnergy", fvc::ddt(rho) + fvc::div(phi));
        if (SmassFluid.set(i))
        {
            contErr += SmassFluid[i];
        }
        if (SporeFluid.set(i))
        {
            contErr -= SporeFluid[i];
        }
        EEqn -= fvm::Sp(contErr, he);
    }

    // Iterasi final tanpa relaksasi (pasangan solveSolid.H, REGISTER 2026-09-22 cacat 1).
    // Pola kamus spt "h.*" ikut cocok dgn hFinal, jadi pengecualian dipasang di kode.
    if (!finalIter)
    {
        EEqn.relax();
    }

    fvOptions.constrain(EEqn);

    if (energyImplicit)
    {
        // Kopling energi implisit: simpan, dirakit di solveEnergyImplicit.H
        kappaCplFluid[i] = turbulence.kappaEff();
        CpCplFluid[i] = thermo.Cp();
        if (energyLinearisationConsistent)
        {
            hCplFluid[i] == thermo.he(thermo.p(), thermo.T());   // h* = h(T*, Y kini): linearisasi antarmuka eksak
        }
        else
        {
            hCplFluid[i] == he;
        }
        simpanEnergi(EEqn, hCplFluid[i]);
    }
    else if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(EEqn);
    }
    else
    {
        EEqn.solve(he.select(finalIter));
        fvOptions.correct(he);

        thermo.correct();
        rad.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}
