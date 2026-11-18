/*---------------------------------------------------------------------------*\
  interfaceLeakPerOuter.H   (DIAGNOSTIK — tanpa fisika baru)

  Mencetak sisa konservasi antarmuka PER ITERASI LUAR di dalam satu langkah
  waktu, bukan sekali di akhir seperti interfaceMassBalance.H.

  KENAPA INI ADA
  --------------
  Kampanye 5,6% kehilangan 40-50% massanya karena kopling antar-region tidak
  konvergen. Diukur di geometri sebenarnya, menaikkan nOuterCorrectors hampir
  tak menolong: 8/16/32 memberi sisa 0,608/0,426/0,323 -- faktor ~0,96 per
  iterasi, dan MEMBURUK. Itu menyisakan satu pertanyaan yang menentukan besar
  kecilnya pekerjaan perbaikan, dan ia tak dapat dijawab dari hasil akhir:

    sisa TURUN monoton dari iterasi ke iterasi
        -> iterasi titik-tetapnya kontraktif. Cukup ganti hitungan tetap
           `nOuterCorrectors` dengan kriteria residual. Perubahan kecil.

    sisa MENDATAR atau berosilasi
        -> iterasinya TIDAK kontraktif; berapa pun nOuter tak akan menolong,
           dan koplingnya harus dibuat IMPLISIT. Perubahan solver besar.

  Satu langkah waktu cukup untuk membedakannya. Karena itu berkas ini ada.

  BATAS
  -----
  Hanya fluks; inventori tidak dihitung (itu tugas interfaceMassBalance.H di
  akhir langkah). Tidak menulis berkas apa pun, tidak menyentuh satu pun medan
  -- jadi ia tak dapat menggeser hasil run mana pun.

  MENYALAKAN: di system/controlDict
      probeInterfaceLeak   true;      // bawaan false
      probeLeakMaxSteps    3;         // berapa langkah waktu pertama dicetak

  BUNTU YANG PERNAH TERJADI (2026-09-05) -- kenapa pencacahnya di luar master
  ---------------------------------------------------------------------------
  Blok ini memanggil gSum(), yaitu KOLEKTIF MPI: setiap rank harus memasukinya
  atau tak satu pun boleh. Versi pertama menaikkan `langkahTercetak` di dalam
  `if (Pstream::master())`, sehingga hanya rank 0 yang pernah melewati ambang
  `probeMaks`. Sesudah itu rank 0 melompati seluruh blok sementara rank 1..n
  tetap memasukinya dan menunggu gSum yang tak akan pernah datang.

  Gejalanya menyesatkan: run `probe_C56/runs/OUT` berhenti PERSIS di oCorr=1
  langkah ke-4 (probeLeakMaxSteps 3), tanpa error, tanpa exit, empat core di
  99,9% CPU selama 3 jam 19 menit -- terlihat seperti run yang lambat, bukan
  run yang mati. Pembandingnya yang membuktikan: `OUT2` memakai batas 500 dan
  hanya menjalani 44 langkah, jadi ia tak pernah menyentuh ambang dan selesai
  normal dengan "End". Jadi yang membunuh adalah AMBANGNYA, bukan fisikanya.

  Aturan yang ditinggalkan: pencacah apa pun yang ikut menentukan apakah suatu
  kolektif dieksekusi TIDAK BOLEH hidup di dalam guard master-only.
  (`interfaceMassBalance.H` aman: seluruh gSum di sana berada di LUAR guard
  master-nya, dan `nPeringatan` tidak menjaga satu kolektif pun.)
\*---------------------------------------------------------------------------*/
{
    static const bool probeAktif =
        runTime.controlDict().getOrDefault<bool>("probeInterfaceLeak", false);
    static const label probeMaks =
        runTime.controlDict().getOrDefault<label>("probeLeakMaxSteps", 3);
    static label langkahTercetak = 0;
    static scalar waktuTerakhir = -GREAT;

    if (probeAktif && langkahTercetak <= probeMaks)
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

        const scalar skala = max(mag(Jf), mag(Js));
        const scalar sisa = (skala > VSMALL) ? mag(Jf - Js)/skala : scalar(0);

        // Pencacah WAJIB maju di SEMUA rank -- lihat catatan BUNTU di kepala
        // berkas. `Info` sendiri sudah hanya mencetak di master, jadi tidak
        // ada guard Pstream::master() di sini dan memang tidak boleh ada.
        if (runTime.timeOutputValue() != waktuTerakhir)
        {
            waktuTerakhir = runTime.timeOutputValue();
            ++langkahTercetak;

            if (langkahTercetak <= probeMaks)
            {
                Info<< "probeLeak: --- t = " << waktuTerakhir
                    << " (langkah " << langkahTercetak << "/" << probeMaks
                    << ") ---" << endl;
            }
        }

        Info<< "probeLeak: oCorr=" << (oCorr + 1) << "/" << nOuterCorr
            << "  Jf=" << Jf << "  Js=" << Js
            << "  sisa=" << sisa << endl;
    }
}
