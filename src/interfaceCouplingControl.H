/*---------------------------------------------------------------------------*\
  interfaceCouplingControl.H

  ITERASI LUAR SAMPAI ANTARMUKA KONSERVATIF -- menggantikan hitungan tetap
  `nOuterCorrectors` dengan kriteria yang memang dituntut fisikanya.

  KENAPA
  ------
  Kopling species fluida-padatan (speciesCoupledMixed) EKSPLISIT: tiap region
  memakai nilai tetangga dari iterasi luar sebelumnya. Iterasi titik-tetap itu
  konvergen ke solusi monolitik, tetapi hanya kalau dijalankan cukup lama.
  Sampai sekarang "cukup lama" diwakili dua angka hasil KALIBRASI:

      maxLambda 0,2          ambang terukur (lulus 0,272 / gagal 0,327, nOuter 8)
      nOuter 32 + dt 0,0107  pasangan yang lolos uji NO32_dt0107

  Keduanya hanya benar pada kasus tempat ia diukur. Kriteria yang tidak
  bergantung kasus adalah neraca itu sendiri:

      sisa = |Jf - Js| / max(|Jf|, |Js|)   <   interfaceTolerance

  dihitung SETELAH tiap iterasi luar (fluida lalu padatan). Begitu terpenuhi,
  iterasi berikutnya ditandai FINAL (solver *Final, relTol 0). Loop baru keluar
  bila sisa < tol juga SESUDAH iterasi final; kalau tidak, iterasi FINAL
  diulang. `nOuterCorrectors` menjadi BATAS PENGAMAN, bukan jumlah iterasi.

  BUG 2026-09-13 (versi c58dbb5): loop keluar TANPA menilai iterasi final.
  Iterasi final masih menggeser fluks, dan interfaceMassBalance.H mengukur
  sesudahnya -> 16 peringatan terisolasi (sisa s/d 7,2e-3) pada sapuan tau
  5,6% walau kebocoran kumulatifnya hanya <= 3,6e-5. Yang dinilai kini sama
  dengan yang disimpan.

  Tol tidak dipilih agar hasil cocok: ia anggaran galat massa (1e-4 = 0,01%
  fluks antarmuka per langkah) dan kecukupannya DIBUKTIKAN dengan uji tol/10
  serta uji dt/2 (doc/tahap1/REGISTER_RUN.md, uji KOP).

  MPI
  ---
  gSum dipanggil di SEMUA rank setiap iterasi (tanpa guard master). Keputusan
  keluar-loop menentukan apakah kolektif berikutnya dijalankan, jadi ia
  di-reduce(andOp) supaya identik di semua rank -- pelajaran buntu
  interfaceLeakPerOuter.H (2026-09-05).

  BAWAAN: interfaceTolerance tidak ditulis (<= 0) -> blok ini tidak berbuat
  apa-apa dan finalIter = (oCorr == nOuterCorr-1) seperti sebelumnya, sehingga
  seluruh case lama bit-identik.

  Hanya species. Kopling ENERGI tetap memakai loop yang sama tetapi tidak
  ikut dinilai di sini.
\*---------------------------------------------------------------------------*/
if (interfaceTol > 0)
{
    scalar Jf = 0, Js = 0;

    forAll(fluidRegions, ri)
    {
        const fvMesh& fmesh = fluidRegions[ri];
        const auto& comp = thermoFluid[ri].composition();
        const PtrList<volScalarField>& Yf = comp.Y();
        const volScalarField& rhof = rhoFluid[ri];

        forAll(Yf, j)
        {
            const volScalarField& Yi = Yf[j];
            const scalar Dmol = DmolList[ri][j];

            forAll(fmesh.boundary(), patchi)
            {
                const fvPatchScalarField& Yp = Yi.boundaryField()[patchi];
                if (Yp.type() != "speciesCoupledMixed") continue;
                const fvPatch& pp = fmesh.boundary()[patchi];
                const scalarField rhop(rhof.boundaryField()[patchi]);
                Jf += gSum(-rhop*Dmol*Yp.snGrad()*pp.magSf());
            }
        }
    }

    forAll(solidRegions, ri)
    {
        const fvMesh& smesh = solidRegions[ri];
        const regionAdsorption& ads = adsorption[ri];
        const scalar epsP = ads.epsilonP();

        forAll(ads.species(), j)
        {
            const volScalarField& Cj = ads.C(j);
            const volScalarField& Dsj = ads.Ds(j);

            forAll(smesh.boundary(), patchi)
            {
                const fvPatchScalarField& Cp = Cj.boundaryField()[patchi];
                if (Cp.type() != "speciesCoupledMixed") continue;
                const fvPatch& pp = smesh.boundary()[patchi];
                const scalarField Dsp(Dsj.boundaryField()[patchi]);
                Js += gSum(epsP*Dsp*Cp.snGrad()*pp.magSf());
            }
        }
    }

    // Skala = max(umpan, |Jf|, |Js|): fraksi umpan yang hilang; tak lenyap saat
    // dinding jenuh. Lihat interfaceFeedFlux.H.
    #include "interfaceFeedFlux.H"
    const scalar skala = max(fluksUmpan, max(mag(Jf), mag(Js)));
    const scalar sisaKop = (skala > VSMALL) ? mag(Jf - Js)/skala : scalar(0);

    kopKonvergen = (sisaKop < interfaceTol);
    reduce(kopKonvergen, andOp<bool>());

    if (kopKonvergen && finalIter)
    {
        Info<< "interfaceCoupling: SELESAI oCorr=" << (oCorr + 1)
            << " sisa=" << sisaKop << " (sesudah FINAL)" << endl;
    }
    else if (kopKonvergen)
    {
        Info<< "interfaceCoupling: konvergen oCorr=" << (oCorr + 1)
            << " sisa=" << sisaKop << " -> iterasi berikut FINAL" << endl;
    }
    else if (finalIter && oCorr < nOuterCorr-1)
    {
        // Kejadian yang dulu lolos tanpa dinilai.
        Info<< "interfaceCoupling: FINAL belum konservatif oCorr=" << (oCorr + 1)
            << " sisa=" << sisaKop << " -> ulang FINAL" << endl;
    }
    else if (oCorr == nOuterCorr-1)
    {
        // Pencacah maju di semua rank (tak ada kolektif yang dijaganya, tapi
        // aturannya dipegang konsisten).
        static label nBatas = 0;
        ++nBatas;
        Info<< "interfaceCoupling: BATAS nOuterCorrectors=" << nOuterCorr
            << " tercapai, sisa=" << sisaKop << " > tol " << interfaceTol
            << " (kejadian ke-" << nBatas << ")" << endl;
        if (nBatas <= 20)
        {
            WarningInFunction
                << "Kopling antarmuka BELUM konvergen dalam "
                << nOuterCorr << " iterasi luar pada t = "
                << runTime.timeOutputValue() << " (sisa " << sisaKop
                << "). Naikkan nOuterCorrectors atau turunkan dt." << endl;
        }
    }
}
