/*---------------------------------------------------------------------------*\
  fluidAdsorptionSource.H

  Kopling massa dua-arah (Bab 5.1 / 8). Menghitung sink massa volumetrik
  SmassFluid[i] [kg/m3/s] di sel-sel yang bersinggungan dengan interface dari
  fluks difusif species adsorbat yang menyeberang ke solid:

      mdot_face = -rho D_i (snGrad Y_i) |Sf|        [kg/s]   (>0 = teradsorpsi)
      Smass_cell = sum_face mdot_face / V_cell

  Tanda bertanda penuh: adsorpsi (Y turun ke dinding) -> Smass>0 (massa hilang
  dari gas); desorpsi/TVSA (Y naik ke dinding) -> Smass<0 (massa ditambahkan).
  Pada regime trace (DAC 400 ppm) nilainya ~0 secara otomatis.

  Patch interface dikenali dari tipe BC species "speciesCoupledMixed".
\*---------------------------------------------------------------------------*/
{
    volScalarField& Smass = SmassFluid[i];
    Smass == dimensionedScalar("0", dimDensity/dimTime, 0);

    scalarField& sCells = Smass.primitiveFieldRef();
    const scalarField& V = mesh.V();

    forAll(Y, j)
    {
        if (j == inertIndex || !composition.active(j))
        {
            continue;
        }

        const volScalarField& Yi = Y[j];
        const scalar Dmol = DmolList[i][j];

        forAll(mesh.boundary(), patchi)
        {
            const fvPatchScalarField& Yp = Yi.boundaryField()[patchi];

            if (Yp.type() != "speciesCoupledMixed")
            {
                continue;
            }

            const fvPatch& pp = mesh.boundary()[patchi];
            const labelUList& fc = pp.faceCells();
            const scalarField magSf(pp.magSf());
            const scalarField rhop(rho.boundaryField()[patchi]);

            const scalarField Yint(Yp.patchInternalField());
            const scalarField snG((Yp - Yint)*pp.deltaCoeffs());  // outward
            const scalarField mdot(-rhop*Dmol*snG*magSf);          // kg/s

            forAll(fc, f)
            {
                sCells[fc[f]] += mdot[f]/V[fc[f]];
            }
        }
    }

    // Relaksasi sumber massa antar-iterasi luar (REGISTER PB2-STRIP-b, 2026-10-10). Pada desorpsi vakum
    // sumber massa <-> tekanan <-> komposisi membentuk daur dua-iterasi: Sads berpindah-pindah antara dua
    // keadaan (-3e5 dan +-6e7 W/m3) tiap iterasi luar, sehingga jumlah iterasi genap berapa pun memberi hasil
    // sama tetapi neraca energi bocor. Rata-rata dengan iterasi sebelumnya (faktor 0,5) mematikan daur itu;
    // saat konvergen nilainya tak berubah. Bawaan 1 (tanpa relaksasi; kasus lama bit-identik).
    {
        const scalar aS = pimple.getOrDefault<scalar>("massSourceRelaxation", 1.0);
        if (aS < 1.0)
        {
            if (!SmassPrevFluid.set(i))
            {
                SmassPrevFluid.set(i, new volScalarField("SmassPrev", Smass));
            }
            else
            {
                sCells = aS*sCells + (1.0 - aS)*SmassPrevFluid[i].primitiveField();
                SmassPrevFluid[i].primitiveFieldRef() = sCells;
            }
        }
    }

    // Tahap V2: bagian gas pori (Spore) dikurangkan di fluidPoreSource.H, sesudah EEqn.
}
