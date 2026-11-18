/*---------------------------------------------------------------------------*\
  massSourceDiag.H  (DIAGNOSTIK, 2026-10-17; mati kecuali `massSourceDiag true` di controlDict)

  Pertanyaan: mengapa Y_CO2 kanal > 1 di semua sel (REGISTER PB2-PORI, pemeriksaan medan C7_A)?
  Per langkah waktu, per region fluida, dari medan AKHIR langkah:
      Mdot   = (INT rho - INT rho0)/dt + SUM_batas phi            [kg/s] sumber massa total yang TEREALISASI
      Adot   = (INT rho Ya - INT rho0 Ya0)/dt + SUM_batas phi Ya  [kg/s] sumber adsorbat (non-inert) yang terealisasi
      Sused  = - INT Smass (+ bagian implisit) pada pEqn terakhir  [kg/s] sumber massa yang DIPAKAI kontinuitas
      Minert = INT rho Y_inert,   Mneg = INT rho (1 - SUM Ya)      [kg]   (Mneg < 0: adsorbat melebihi massa gas)
  Adot > Mdot terus-menerus -> Y_adsorbat naik melewati 1. Mdot != Sused -> kontinuitas tak terealisasi seperti
  yang dipakai pEqn (mis. rho berubah sesudah pEqn). Keluaran postProcessing/massSourceDiag.dat (master).
\*---------------------------------------------------------------------------*/
if (runTime.controlDict().getOrDefault<bool>("massSourceDiag", false))
{
    forAll(fluidRegions, ri)
    {
        const fvMesh& fmesh = fluidRegions[ri];
        const rhoReactionThermo& th = thermoFluid[ri];
        const PtrList<volScalarField>& Yf = th.composition().Y();
        const label iIn = th.composition().species().find(th.get<word>("inertSpecie"));
        const volScalarField& rhof = rhoFluid[ri];
        const surfaceScalarField& phiF = phiFluid[ri];
        const scalar dt = runTime.deltaTValue();

        scalar M = gSum(rhof.primitiveField()*fmesh.V());
        scalar M0 = gSum(rhof.oldTime().primitiveField()*fmesh.V());
        scalar A = 0, A0 = 0, Mneg = 0;
        scalarField sumYa(fmesh.nCells(), Zero);
        scalar out = 0, outA = 0;
        forAll(fmesh.boundary(), patchi)
        {
            if (fmesh.boundary()[patchi].coupled()) continue;
            out += sum(phiF.boundaryField()[patchi]);
        }
        forAll(Yf, j)
        {
            if (j == iIn) continue;
            const volScalarField& Yj = Yf[j];
            sumYa += Yj.primitiveField();
            A += gSum(rhof.primitiveField()*Yj.primitiveField()*fmesh.V());
            A0 += gSum(rhof.oldTime().primitiveField()*Yj.oldTime().primitiveField()*fmesh.V());
            forAll(fmesh.boundary(), patchi)
            {
                if (fmesh.boundary()[patchi].coupled()) continue;
                outA += sum(phiF.boundaryField()[patchi]*Yj.boundaryField()[patchi]);
            }
        }
        reduce(out, sumOp<scalar>());
        reduce(outA, sumOp<scalar>());
        Mneg = gSum(rhof.primitiveField()*(1 - sumYa)*fmesh.V());
        const scalar Min = gSum(rhof.primitiveField()*Yf[iIn].primitiveField()*fmesh.V());
        const scalar Mdot = (M - M0)/dt + out;
        const scalar Adot = (A - A0)/dt + outA;
        const scalar Sused = -SmassUsedTot[ri];   // domainIntegrate sudah global
        const scalar maxSumYa = gMax(sumYa);      // kolektif: dipanggil semua proses
        FixedList<scalar, 6> pd(poreDiagFluid[ri]);
        for (scalar& x : pd) { reduce(x, sumOp<scalar>()); }
        if (Pstream::master())
        {
            static autoPtr<OFstream> osPtr;
            if (!osPtr)
            {
                mkDir(runTime.globalPath()/"postProcessing");
                osPtr.reset(new OFstream(runTime.globalPath()/"postProcessing"/"massSourceDiag.dat"));
                osPtr() << "# time  region  Mdot[kg/s]  Adot[kg/s]  Sused[kg/s]  M[kg]  A[kg]  Minert[kg]  Mneg[kg]  maxSumYa  NtPore[mol/s]  Fadv[kg/s]  Fdiff[kg/s]  molAds[mol/s]  mInert[kg/s]  mInertMentah[kg/s]" << endl;
            }
            osPtr() << runTime.timeOutputValue() << ' ' << ri << ' ' << Mdot << ' ' << Adot << ' ' << Sused << ' '
                    << M << ' ' << A << ' ' << Min << ' ' << Mneg << ' ' << maxSumYa << ' ' << pd[0] << ' ' << pd[1] << ' ' << pd[2] << ' ' << pd[3] << ' ' << pd[4] << ' ' << pd[5] << endl;
        }
    }
}
