/*---------------------------------------------------------------------------*\
  axialBins.H   (Tahap 2.1 -- DIAGNOSTIK, hanya membaca field)

  Rata-rata per potongan aksial untuk mesin ekstraksi parameter (Tahap 2):
  sel padatan & fluida dikelompokkan menurut z pusat sel ke nBins potongan
  seragam [zmin, zmax] (batas global mesh fluida).

  Padatan, rata-rata VOLUME per potongan:
      qbar_i    = INT q_i dV / INT dV          [mol/kg]
      qsbar_i   = INT qStar_i dV / INT dV      [mol/kg]  (pembanding saja: rata2 q*(C lokal)
                                                         != q*(c bulk) -- Tahap 2.2 WAJIB
                                                         menghitung q* dari c bulk)
      Cbar_i    = INT C_i dV / INT dV          [kg/m3]   gas pori
      dqdt_i    = INT (q_i - q_i^old)/dt dV / INT dV     laju SAMA dgn yg dipakai Sads
      Ts        = INT T dV / INT dV
  Fluida:
      Yv_i, Tv  = rata-rata volume (~ rata-rata luas, mesh aksial seragam)
      Yc_i, Tc  = rata-rata berbobot FLUKS (cup-mixing), w = rho U_z:
                  INT w Y dV / INT w dV   -- cup-mixing dirata-rata sepanjang potongan.
                  Profil laminar parabolik: Yc != Yv; memakai Yv menggeser k secara
                  sistematis (RENCANA 2.1).
      p, rho    = rata-rata volume
  Uji konsistensi (Tahap 2.1): SUM_bins rho_s M_i qbar_i Vs_bin == 'adsorbed' di
  interfaceBalance.dat -- kolom Vs ditulis supaya bisa dijumlah ulang.

  Aktif bila system/axialBins ada:
      nBins     40;        // jumlah potongan
      interval  10;        // [s] waktu simulasi antar keluaran (terpisah dari writeTime:
                           //  tak menambah direktori/inode, tak memicu cacat tulis MPI)
  Keluaran: postProcessing/axialBins.dat (master), satu baris per potongan per keluaran.
  Aturan MPI: jumlah LOKAL, satu reduce.
\*---------------------------------------------------------------------------*/
if (axialBinsN > 0)
{
    const scalar tNow = runTime.value();
    const label slot = label(std::floor(tNow/axialBinsInterval + 1e-9));
    static label slotLast = -1;

    if (slot != slotLast)
    {
        slotLast = slot;

        const regionAdsorption& adsB = adsorption[0];
        const wordList& spB = adsB.species();
        const label nSp = spB.size();
        const label nBase = 8;          // Vs TsV Vf TfV wV wTV pV rhoV
        const label nPer = nBase + 6*nSp; // per sp: qV qsV CV dqV YV wYV
        const label nb = axialBinsN;
        List<scalar> a(nb*nPer, Zero);

        auto binOf = [&](const scalar z)
        {
            label k = label((z - axialBinsZmin)/(axialBinsZmax - axialBinsZmin)*nb);
            return max(label(0), min(nb - 1, k));
        };

        // ---------------- padatan ----------------
        {
            const fvMesh& sm = solidRegions[0];
            const scalarField& V = sm.V();
            const vectorField& Cc = sm.C();
            const scalarField& Ts = thermos[0].T().primitiveField();
            const scalar dt = runTime.deltaTValue();

            forAll(V, c)
            {
                const label o = binOf(Cc[c].z())*nPer;
                a[o + 0] += V[c];
                a[o + 1] += Ts[c]*V[c];
            }
            forAll(spB, i)
            {
                const scalarField& q = adsB.q(i).primitiveField();
                const scalarField& qo = adsB.q(i).oldTime().primitiveField();
                const scalarField& C = adsB.C(i).primitiveField();
                const scalarField& qs =
                    sm.lookupObject<volScalarField>("qStar." + spB[i]).primitiveField();
                forAll(V, c)
                {
                    const label o = binOf(Cc[c].z())*nPer + nBase + 6*i;
                    a[o + 0] += q[c]*V[c];
                    a[o + 1] += qs[c]*V[c];
                    a[o + 2] += C[c]*V[c];
                    a[o + 3] += (q[c] - qo[c])/dt*V[c];
                }
            }
        }

        // ---------------- fluida ----------------
        {
            const fvMesh& fm = fluidRegions[0];
            const scalarField& V = fm.V();
            const vectorField& Cc = fm.C();
            const rhoReactionThermo& th = thermoFluid[0];
            const scalarField& Tf = th.T().primitiveField();
            const scalarField& p = th.p().primitiveField();
            const scalarField& rho = rhoFluid[0].primitiveField();
            const vectorField& U = UFluid[0].primitiveField();
            const auto& Yf = th.composition().Y();

            forAll(V, c)
            {
                const label o = binOf(Cc[c].z())*nPer;
                const scalar w = rho[c]*U[c].z();
                a[o + 2] += V[c];
                a[o + 3] += Tf[c]*V[c];
                a[o + 4] += w*V[c];
                a[o + 5] += w*Tf[c]*V[c];
                a[o + 6] += p[c]*V[c];
                a[o + 7] += rho[c]*V[c];
            }
            forAll(spB, i)
            {
                const label j = th.composition().species().find(spB[i]);
                if (j < 0) continue;
                const scalarField& Y = Yf[j].primitiveField();
                forAll(V, c)
                {
                    const label o = binOf(Cc[c].z())*nPer + nBase + 6*i;
                    a[o + 4] += Y[c]*V[c];
                    a[o + 5] += rho[c]*U[c].z()*Y[c]*V[c];
                }
            }
        }

        Pstream::listCombineReduce(a, plusEqOp<scalar>());

        if (Pstream::master())
        {
            static autoPtr<OFstream> osB;
            if (!osB)
            {
                const fileName ppDir(runTime.globalPath()/"postProcessing");
                mkDir(ppDir);
                osB.reset(new OFstream(ppDir/"axialBins.dat"));
                osB().precision(12);   // dq/dt & jumlah ulang butuh > 6 digit (uji B21)
                *osB<< "# axialBins.H (Tahap 2.1): nBins " << nb << ", z " << axialBinsZmin
                    << " .. " << axialBinsZmax << " m, rho_s " << adsB.rhoS();
                forAll(spB, i) { *osB<< ", M_" << spB[i] << " " << adsB.M(i); }
                *osB<< nl << "# time bin zc Vs Ts Vf Tv Tc p rho";
                forAll(spB, i)
                {
                    const word& s = spB[i];
                    *osB<< " qbar_" << s << " qsbar_" << s << " Cbar_" << s
                        << " dqdt_" << s << " Yv_" << s << " Yc_" << s;
                }
                *osB<< endl;
            }
            OFstream& os = *osB;
            const scalar dz = (axialBinsZmax - axialBinsZmin)/nb;
            for (label b = 0; b < nb; ++b)
            {
                const label o = b*nPer;
                const scalar Vs = max(a[o + 0], VSMALL);
                const scalar Vf = max(a[o + 2], VSMALL);
                const scalar W = a[o + 4];
                const scalar Wd = (mag(W) > VSMALL) ? W : VSMALL;
                os  << runTime.timeOutputValue() << ' ' << b << ' '
                    << axialBinsZmin + (b + 0.5)*dz << ' ' << a[o + 0] << ' '
                    << a[o + 1]/Vs << ' ' << a[o + 2] << ' '
                    << a[o + 3]/Vf << ' ' << a[o + 5]/Wd << ' '
                    << a[o + 6]/Vf << ' ' << a[o + 7]/Vf;
                for (label i = 0; i < nSp; ++i)
                {
                    const label s = o + nBase + 6*i;
                    os  << ' ' << a[s + 0]/Vs << ' ' << a[s + 1]/Vs << ' ' << a[s + 2]/Vs
                        << ' ' << a[s + 3]/Vs << ' ' << a[s + 4]/Vf << ' ' << a[s + 5]/Wd;
                }
                os  << nl;
            }
            os.flush();
        }
    }
}
