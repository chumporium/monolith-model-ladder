/*---------------------------------------------------------------------------*\
  solveSolid.H  (pluggable adsorption version)

  Persamaan energi solid dengan source panas adsorpsi (Bab 5.2):

    d(betav rho h)/dt = div(betav alpha grad h) + S_ads

  S_ads = sum_i deltaH_i rho_s dq_i/dt   [W/m3]  disediakan regionAdsorption
  (konvensi: deltaH_i > 0 eksotermik -> adsorpsi memanaskan solid).

  S_ads sudah dihitung di solveSolidSpecies.H (adsorptionI.correct) yang
  dipanggil sebelum file ini.
\*---------------------------------------------------------------------------*/
{
    const volScalarField& Sads = adsorptionI.Sads();

    Info<< "Sads - " << solidRegions[i].name()
        << ": min/max = " << min(Sads).value()
        << " / " << max(Sads).value() << " W/m3" << endl;

    fvScalarMatrix hEqn
    (
        fvm::ddt(betav*rho, h)
        - (
            thermo.isotropic()
            ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
            : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
          )
       ==
        Sads
      + fvOptions(rho, h)
    );

    if (hvLoss[i] > 0)
    {
        const dimensionedScalar hvD("hv", dimPower/dimVolume/dimTemperature, hvLoss[i]);
        const dimensionedScalar TaD("Tamb", dimTemperature, TambLoss[i]);
        const volScalarField& Ts = thermo.T();
        // Sumber ruas kanan R = -hv (T - Tamb) = -hv/cp h + hv (h*/cp - (T* - Tamb)).
        // hEqn menyimpan (kiri - kanan), jadi  hEqn -= R :
        //   -(-fvm::Sp(hv/cp, h)) -> +=  ;  -(+eksplisit) -> -=
        // (konvensi sama dgn radiationModel::Sh: sink ditulis "- fvm::Sp" di ruas kanan).
        // Versi pertama (0fe92dd) tandanya terbalik -> memanaskan ~hv (T - Tamb); tertangkap
        // gerbang neraca energi dlm < 500 s simulasi (errE +0,05...+0,10).
        hEqn += fvm::Sp(hvD/cp, h);
        hEqn -= hvD*(h/cp - (Ts - TaD));
    }

    // DIAG-EN (2026-09-22): lacak ke mana Sads pergi. H = INT betav rho h dV.
    //   Hawal  : sebelum solve (sama dgn Hakhir iterasi sebelumnya bila tak ada yg mereset h)
    //   Hold   : dari h.oldTime() -- acuan langkah waktu
    //   dHsrc  : INT Sads dV * dt   (yang seharusnya ditambahkan relatif Hold)
    const scalar diagHold = gSum(betav.primitiveField()*rho.primitiveField()
        *h.oldTime().primitiveField()*mesh.V());
    const scalar diagHawal = gSum(betav.primitiveField()*rho.primitiveField()
        *h.primitiveField()*mesh.V());
    const scalar diagSrc = gSum(Sads.primitiveField()*mesh.V())*runTime.deltaTValue();

    // Relaksasi energi HANYA di iterasi luar non-final (REGISTER 2026-09-22, cacat 1).
    // relax() memaksa D = max(D0, sum|off|)/alpha dan menambah (D - D0) psi ke sumber:
    // bila dt >> waktu konduksi sel (dx^2/alpha ~ 1e-3 s), Laplacian mendominasi D sehingga
    // hanya ~tau_sel/dt dari Sads diterapkan per solve. Dulu relax() dipanggil SEBELUM
    // setFinalIteration, jadi iterasi final pun direlaksasi dan sisa sumber dibuang saat
    // langkah waktu maju (padatan tertutup: 2,3e-4 Sads diterapkan). Iterasi final kini
    // selalu tanpa relaksasi, apa pun isi relaxationFactors -> neraca padatan eksak.
    if (!finalIter)
    {
        hEqn.relax();
    }
    fvOptions.constrain(hEqn);

    if (energyImplicit)
    {
        kappaCplSolid[i] = betav*thermo.kappa();
        CpCplSolid[i] = cp;
        hCplSolid[i] == h;
        simpanEnergi(hEqn, hCplSolid[i]);
    }
    else if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(hEqn);
    }
    else
    {
        Info<< "\nSolving for solid region "<< solidRegions[i].name() << endl;

        if (finalIter) mesh.data().setFinalIteration(true);

        hEqn.solve(h.select(finalIter));
        fvOptions.correct(h);
        thermo.correct();

        const scalar diagHakhir = gSum(betav.primitiveField()*rho.primitiveField()
            *h.primitiveField()*mesh.V());
        Info<< "DIAG-EN solid: Hold=" << diagHold << " Hawal-Hold=" << diagHawal - diagHold
            << " Hakhir-Hold=" << diagHakhir - diagHold << " SadsdtV=" << diagSrc
            << " rasio=" << (diagHakhir - diagHold)/max(diagSrc, VSMALL)
            << " dt=" << runTime.deltaTValue() << endl;

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;

        if (finalIter) mesh.data().setFinalIteration(false);
    }
}
