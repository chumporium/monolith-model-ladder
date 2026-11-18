/*---------------------------------------------------------------------------*\
  createFluidAdsorptionFields.H

  Sisi fluida memakai FRAKSI MASSA Y_i native thermo multikomponen kompresibel
  (Bab 5.1 / 10.3). Konsentrasi C_i = rho Y_i terikat EOS gas -> prasyarat TVSA
  emergent (Bab 5.4). TIDAK ada passive scalar konsentrasi yang bebas-tekanan.

  Disiapkan di sini:
    - DmolList[i][j] : difusivitas molekular species j di region i [m2/s],
      dibaca dari thermophysicalProperties subdict speciesDiffusivity
      (default DmolDefault). Dipakai YEqn.H.
    - SmassFluid[i]  : sink massa adsorpsi interface [kg/m3/s] untuk kopling
      dua-arah (Bab 5.1). Dihitung di fluidAdsorptionSource.H, dipakai
      rhoEqn.H & pEqn.H. ~0 otomatis pada regime trace.
\*---------------------------------------------------------------------------*/

// Difusivitas molekular default [m2/s] bila species tak terdaftar
const scalar DmolDefault = 1.6e-05;

// Turbulent Schmidt number untuk bagian turbulen difusi species.
// Dibaca dari speciesDiffusivity { Sct ...; } region fluida pertama; bawaan 0,7 (nilai lama).
// Tak berpengaruh pada kasus laminar (mut = 0). Untuk analogi panas-massa pakai Sct = Prt.
const scalar SctFluid =
(
    fluidRegions.size() && thermoFluid[0].found("speciesDiffusivity")
  ? thermoFluid[0].subDict("speciesDiffusivity").getOrDefault<scalar>("Sct", 0.7)
  : 0.7
);
Info<< nl << "*** Schmidt turbulen Sct = " << SctFluid << endl;

// Difusivitas molekular per region per species (urutan composition.species())
List<scalarList> DmolList(fluidRegions.size());

// Sink massa adsorpsi [kg/m3/s] (kopling dua-arah)
PtrList<volScalarField> SmassFluid(fluidRegions.size());

// Tahap V2 (poreAdvection): sumber species [kg/m3/s] dari gas pori yang masuk kanal oleh fluks Darcy.
// Kosong (tak diset) bila tak dipakai -> YEqn tak berubah.
List<PtrList<volScalarField>> SadvFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    const rhoReactionThermo& thermo = thermoFluid[i];
    const auto& composition = thermo.composition();
    const auto& specieNames = composition.species();

    dictionary diffDict;
    if (thermo.found("speciesDiffusivity"))
    {
        diffDict = thermo.subDict("speciesDiffusivity");
    }

    DmolList[i].setSize(specieNames.size(), DmolDefault);
    forAll(specieNames, j)
    {
        DmolList[i][j] =
            diffDict.getOrDefault<scalar>(specieNames[j], DmolDefault);
    }

    Info<< "\n*** Fluid region " << fluidRegions[i].name()
        << " species diffusivities [m2/s]:" << nl;
    forAll(specieNames, j)
    {
        Info<< "      " << specieNames[j] << " = " << DmolList[i][j] << nl;
    }

    SmassFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Smass",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("0", dimDensity/dimTime, 0)
        )
    );
}
// Tahap V2: sumber massa [kg/m3/s] dan entalpi [W/m3] dari gas pori yang menyeberang (fluidPoreSource.H).
PtrList<volScalarField> SporeFluid(fluidRegions.size());
PtrList<volScalarField> SentPoreFluid(fluidRegions.size());
// poreSourceImplicit: G = dS_pori/dp kanal [kg/m3/s/Pa] per sel fluida dan p_rgh saat Spore dibaca (p*), lihat pEqn.H.
PtrList<volScalarField> GporeFluid(fluidRegions.size());
PtrList<volScalarField> pStarPoreFluid(fluidRegions.size());
// Skalar modus-peluruhan (decayScalar.H): satu daftar per region fluida; kosong bila tak dipakai.
List<PtrList<volScalarField>> decFields(fluidRegions.size());
// Sumber massa iterasi luar sebelumnya (massSourceRelaxation < 1); kosong bila tak dipakai.
PtrList<volScalarField> SmassPrevFluid(fluidRegions.size());
// massSourceDiag: integral sumber massa [kg/s] yang DIPAKAI pEqn terakhir (Smass + bagian implisit G(p - p*)); hanya diagnostik.
List<scalar> SmassUsedTot(fluidRegions.size(), Zero);
// massSourceDiag: rincian fluks pori pembacaan terakhir fluidPoreSource.H (lokal per proses): NtPore [mol/s], adveksi adsorbat [kg/s],
// difusi adsorbat [kg/s] dan adsorbat [mol/s] (adveksi + difusi), sisa inert [kg/s].
List<FixedList<scalar, 6>> poreDiagFluid(fluidRegions.size(), FixedList<scalar, 6>(Zero));   // [5] = sisa inert SEBELUM pembatas
