/*---------------------------------------------------------------------------*\
  interfaceAitken.H   -- relaksasi dinamis Aitken kopling species antarmuka

  AKTIF hanya bila  PIMPLE { interfaceRelaxation aitken; }  (bawaan none ->
  jalur lama bit-identik).

  MASALAH (Grande 2006, 2026-09-17)
  ---------------------------------
  Kopling antar-region adalah iterasi titik-tetap Gauss-Seidel: fluida
  diselesaikan dgn C padatan dari iterasi sebelumnya, lalu padatan dgn Y fluida
  yang baru. Iterasi itu hanya kontraktif untuk Lambda kecil (adsorptionCouplingNo.H;
  C56 5,6%: gagal di atas ~0,3). Zeolit 4A + propilena: Lambda 2,4 dan
  sebanding dt -> sisa antarmuka membesar 3e-4 -> 7e-2 tiap langkah, BATAS
  100 iterasi. Membatasi dt lewat maxLambda butuh dt ~1,5e-3 s untuk 60 000 s
  (~4e7 langkah): tidak dapat dijalankan.

  SKEMA
  -----
  Variabel titik-tetap x = nilai C padatan di sel-sel antarmuka yang DILIHAT
  fluida (field "<C>Aitken", dibaca speciesCoupledMixed sisi fluida).
  Satu iterasi luar memetakan x_k -> F(x_k) = C antarmuka hasil solve padatan.

      r_k     = F(x_k) - x_k
      omega_k = -omega_{k-1} r_{k-1}.(r_k - r_{k-1}) / |r_k - r_{k-1}|^2   (k >= 1)
      x_{k+1} = x_k + omega_k r_k

  (Irons & Tuck 1969; Kuttler & Wall 2008, Comput. Mech. 43:61.) omega_0 tiap
  langkah waktu = |omega| terakhir langkah sebelumnya, dibatasi <= 1.
  Tidak ada parameter bebas: omega dihitung dari dua residu terakhir.

  Konservasi TIDAK diserahkan ke skema ini: iterasi tetap berhenti hanya bila
  |Jf - Js|/skala < interfaceTolerance (interfaceCouplingControl.H), dan neraca
  global diperiksa di interfaceMassBalance.H. Pada titik tetap r = 0 -> x = C
  antarmuka padatan -> BC identik dgn skema tanpa relaksasi.

  MPI: jumlah lokal per patch non-prosesor, reduce SEKALI (aturan di
  interfaceFeedFlux.H).
\*---------------------------------------------------------------------------*/

if (interfaceAitken)
{
    scalar sumR2 = 0;       // |r_k|^2
    scalar sumF2 = 0;       // |F|^2 (skala laporan)
    scalar sumRpDr = 0;     // r_{k-1}.(r_k - r_{k-1})
    scalar sumDr2 = 0;      // |r_k - r_{k-1}|^2

    // --- lintasan 1: residu dan hasil kali dalam ---
    forAll(solidRegions, ri)
    {
        const fvMesh& smesh = solidRegions[ri];
        const regionAdsorption& ads = adsorption[ri];

        forAll(ads.species(), j)
        {
            const volScalarField& Cj = ads.C(j);
            volScalarField& xF =
                smesh.lookupObjectRef<volScalarField>(Cj.name() + "Aitken");
            volScalarField& rF =
                smesh.lookupObjectRef<volScalarField>(Cj.name() + "AitkenR");

            forAll(smesh.boundary(), patchi)
            {
                if (smesh.boundary()[patchi].coupled()) continue;
                if (Cj.boundaryField()[patchi].type() != "speciesCoupledMixed")
                {
                    continue;
                }

                const scalarField F(Cj.boundaryField()[patchi].patchInternalField());
                const scalarField& x = xF.boundaryField()[patchi];
                const scalarField r(F - x);

                sumR2 += sum(r*r);
                sumF2 += sum(F*F);

                if (oCorr > 0)
                {
                    const scalarField& rPrev = rF.boundaryField()[patchi];
                    const scalarField dr(r - rPrev);
                    sumRpDr += sum(rPrev*dr);
                    sumDr2 += sum(dr*dr);
                }
            }
        }
    }

    reduce(sumR2, sumOp<scalar>());
    reduce(sumF2, sumOp<scalar>());
    reduce(sumRpDr, sumOp<scalar>());
    reduce(sumDr2, sumOp<scalar>());

    if (oCorr == 0)
    {
        aitkenOmega = min(mag(aitkenOmega), scalar(1));
    }
    else if (sumDr2 > VSMALL)
    {
        aitkenOmega = -aitkenOmega*sumRpDr/sumDr2;
    }
    // sumDr2 == 0: residu tak berubah (sudah di titik tetap) -> omega tetap.

    // --- lintasan 2: x <- x + omega r, simpan r ---
    forAll(solidRegions, ri)
    {
        const fvMesh& smesh = solidRegions[ri];
        const regionAdsorption& ads = adsorption[ri];

        forAll(ads.species(), j)
        {
            const volScalarField& Cj = ads.C(j);
            volScalarField& xF =
                smesh.lookupObjectRef<volScalarField>(Cj.name() + "Aitken");
            volScalarField& rF =
                smesh.lookupObjectRef<volScalarField>(Cj.name() + "AitkenR");

            forAll(smesh.boundary(), patchi)
            {
                if (smesh.boundary()[patchi].coupled()) continue;
                if (Cj.boundaryField()[patchi].type() != "speciesCoupledMixed")
                {
                    continue;
                }

                const scalarField F(Cj.boundaryField()[patchi].patchInternalField());
                scalarField& x = xF.boundaryFieldRef()[patchi];
                const scalarField r(F - x);

                x += aitkenOmega*r;
                rF.boundaryFieldRef()[patchi] = r;
            }
        }
    }

    Info<< "interfaceAitken: oCorr=" << oCorr << " omega=" << aitkenOmega
        << " |r|/|C|=" << Foam::sqrt(sumR2/max(sumF2, VSMALL)) << endl;
}
