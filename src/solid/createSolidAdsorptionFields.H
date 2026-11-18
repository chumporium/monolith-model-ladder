/*---------------------------------------------------------------------------*\
  createSolidAdsorptionFields.H

  Membuat koordinator adsorpsi pluggable (regionAdsorption) untuk tiap solid
  region. Koordinator memiliki & menyelesaikan:
    - C.<sp>  konsentrasi pori [kg/m3]   (MUST_READ dari 0/<region>/)
    - q.<sp>  loading [mol/kg]           (READ_IF_PRESENT, AUTO_WRITE)
    - qStar, Ds, Sads

  Komponen dipilih saat runtime dari constant/<region>/adsorptionProperties:
    isothermModel (Langmuir/Toth/Sips/DSL), mixingRule (ideal/extendedLangmuir/
    IAST/coAdsorption), kineticClosure (localEquilibrium/LDF/dualSite).

  Lihat library adsorptionModels.
  (regionAdsorption.H di-include di file-scope adsorptionMultiRegionFoam.C,
   bukan di sini, karena file ini di-include di dalam main().)
\*---------------------------------------------------------------------------*/

PtrList<regionAdsorption> adsorption(solidRegions.size());

forAll(solidRegions, i)
{
    Info<< "\n*** Constructing adsorption model for solid region "
        << solidRegions[i].name() << nl << endl;

    adsorption.set
    (
        i,
        new regionAdsorption(solidRegions[i], thermos[i])
    );
}

// Tekanan pori transien (updatePorePressure.H, mode darcy): satu medan per region padat,
// dibuat saat pertama dipakai. Tak dipakai (mode frozen/mapped) -> tetap kosong.
PtrList<volScalarField> pPoreFields(solidRegions.size());
// Tahap V2 (poreAdvection yes): fluks volumetrik Darcy [m3/s] dan fluks mol total [mol/s] per muka,
// "phiDarcy" terdaftar di mesh padat; fluks mol mentah "NtPoreRaw".
PtrList<surfaceScalarField> phiDarcyFields(solidRegions.size());
PtrList<surfaceScalarField> NtPoreFields(solidRegions.size());
// Konduktansi mol pori per muka [mol/s/Pa] = Kp |S_f| deltaCoeffs, terdaftar "GPore" di mesh padat; dipakai sisi fluida
// untuk sumber massa pori IMPLISIT di pEqn (PIMPLE fluida: poreSourceImplicit; RANCANGAN_PEQN_IMPLISIT.md).
PtrList<surfaceScalarField> GporeFields(solidRegions.size());
// Penerbit fluks antarmuka konservatif (poreFluxCarry.H): "NtPore" dan "FadvPore.<sp>" terdaftar di mesh padat.
// NtPoreFields di atas kini medan MENTAH ("NtPoreRaw").
List<poreFluxCarry> NtCarry(solidRegions.size());
List<PtrList<poreFluxCarry>> FadvCarry(solidRegions.size());
// Fluks adveksi species yang konsisten-mol (updatePorePressure.H): phiC = phiDarcy p_muka/p_hulu, muka antarmuka
// keluar di-nol-kan dan digantikan pembuangan implisit di sel (poreSink) dengan koefisien muka phiCOut.
PtrList<surfaceScalarField> phiCFields(solidRegions.size());
PtrList<surfaceScalarField> phiCOutFields(solidRegions.size());
PtrList<volScalarField> poreSinkFields(solidRegions.size());
// Lag pembukuan sumber species gas pori: 2 = dipakai YEqn iterasi luar berikutnya (bawaan);
// 1 = korektor species final aktif (solveFluid.H), sumber dipakai pada iterasi pembacaannya.
label fadvLag = 2;
// Gas kanal yang MASUK pori lewat antarmuka (phiDarcy < 0): komposisi hulu = sel kanal. Per region padat, per species:
// sumber volumetrik di sel antarmuka [kg/m3/s] dan fluks massa per muka [kg/s] (negatif = meninggalkan kanal).
List<PtrList<volScalarField>> poreInSrcFields(solidRegions.size());
List<PtrList<surfaceScalarField>> poreInFluxFields(solidRegions.size());
