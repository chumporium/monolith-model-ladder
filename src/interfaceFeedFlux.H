/*---------------------------------------------------------------------------*\
  interfaceFeedFlux.H

  Mendefinisikan  scalar fluksUmpan  [kg/s] = laju massa adsorbat yang MASUK
  domain fluida lewat muka berfluks-masuk (phi < 0) di patch non-antarmuka:

      fluksUmpan = sum ( -min(phi, 0) * Y )   untuk species ber-speciesCoupledMixed

  KENAPA
  ------
  `sisa` antarmuka semula dinilai relatif thd max(|Jf|,|Js|). Fluks itu lenyap
  saat dinding jenuh (T30_u008: 3,0e-9 -> 2,3e-12 kg/s), sehingga tol yang sama
  menjadi 1300x lebih ketat secara mutlak tepat ketika hampir tak ada CO2 yang
  berpindah -> iterasi 26 -> 40-48, BATAS, dan peringatan palsu (2026-09-14).
  Skala baru max(fluksUmpan, |Jf|, |Js|) = "fraksi umpan yang hilang": sama
  dengan skala lama selama front lewat (Jf ~ umpan), tak lenyap saat jenuh, dan
  tetap |Jf|/|Js| bila fluks antarmuka melebihi umpan (mis. desorpsi).

  BUNTU MPI 2026-09-14 (versi ab10a79) -- kenapa jumlahnya LOKAL lalu reduce
  -----------------------------------------------------------------------
  Versi pertama memanggil gSum di dalam loop atas SEMUA patch. Patch prosesor
  (processorPolyPatch) jumlah dan urutannya BERBEDA per rank, sehingga tiap rank
  memanggil kolektif sebanyak yang berbeda -> semua rank saling menunggu. Tiga
  run u = 0,0009 macet di langkah pertama sesudah dilanjutkan, dua siklus 8 jam
  habis oleh detektor HANG. Aturan: jangan panggil kolektif di dalam loop yang
  panjangnya bergantung rank; jumlahkan lokal, reduce SEKALI.
\*---------------------------------------------------------------------------*/
scalar fluksUmpan = 0;

forAll(fluidRegions, ri)
{
    const fvMesh& fmesh = fluidRegions[ri];
    const PtrList<volScalarField>& Yf = thermoFluid[ri].composition().Y();
    const surfaceScalarField& phiF = phiFluid[ri];

    forAll(Yf, j)
    {
        const volScalarField& Yi = Yf[j];

        bool isAdsorbate = false;
        forAll(fmesh.boundary(), patchi)
        {
            if (Yi.boundaryField()[patchi].type() == "speciesCoupledMixed")
            {
                isAdsorbate = true;
            }
        }
        if (!isAdsorbate) continue;

        forAll(fmesh.boundary(), patchi)
        {
            // Patch prosesor/cyclic bukan batas domain: tak ada umpan di sana.
            if (fmesh.boundary()[patchi].coupled()) continue;

            const fvPatchScalarField& Yp = Yi.boundaryField()[patchi];
            if (Yp.type() == "speciesCoupledMixed") continue;
            const scalarField& phip = phiF.boundaryField()[patchi];
            fluksUmpan += sum(-min(phip, scalar(0))*Yp);   // LOKAL, tanpa kolektif
        }
    }
}

// Satu-satunya kolektif di berkas ini: dipanggil tepat sekali di semua rank.
reduce(fluksUmpan, sumOp<scalar>());
