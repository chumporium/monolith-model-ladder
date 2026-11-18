/*---------------------------------------------------------------------------*\
  createEnergyCoupling.H  -- sekali saat start

  PIMPLE { energyInterfaceCoupling partitioned | implicit; }   (bawaan partitioned)

  partitioned : jalur lama (BC suhu mixed + iterasi luar), bit-identik.
  implicit    : EEqn fluida + hEqn padatan dirakit SATU matriks tiap iterasi
                luar (solveEnergyImplicit.H). Koefisien antarmuka ditulis
                energyCoupledImplicit::manipulateMatrix pada field bantu hCpl.

  Kenapa (REGISTER 2026-09-22, cacat 3): waktu konduksi antarmuka ~0,03 s << dt
  ~2 s; kopling partisi konvergen sangat lambat (30 iterasi luar: 60% panas
  adsorpsi masih tak terkirim) dan sisanya dibuang tiap langkah.

  Field per region (terdaftar, NO_WRITE):
    hCpl      salinan he; BC antarmuka energyCoupledImplicit (useImplicit),
              patch terkopel (processor) = tipe he, lainnya calculated
    kappaCpl  gas: turbulence.kappaEff(); padatan: betav*kappa     [W/m/K]
    CpCpl     Cp                                                   [J/kg/K]

  Syarat: satu region fluida + satu padatan; kedua sisi muka antarmuka di
  prosesor sama (dekomposisi aksial); tidak digabung dgn coupledEnergyField.
  Disarankan Tref = T umpan di hConst: nilai h absolut besar (T - Tstd) ikut
  menentukan normalisasi residual solver rakitan (ENmono: 0 iterasi).
\*---------------------------------------------------------------------------*/

word energyCouplingMode("partitioned");
{
    fvSolution solDict(runTime);
    energyCouplingMode =
        solDict.subDict("PIMPLE").getOrDefault<word>("energyInterfaceCoupling", "partitioned");

    if (energyCouplingMode != "partitioned" && energyCouplingMode != "implicit")
    {
        FatalIOErrorInFunction(solDict.subDict("PIMPLE"))
            << "energyInterfaceCoupling '" << energyCouplingMode
            << "' tidak dikenal. Pilihan: partitioned, implicit."
            << exit(FatalIOError);
    }
}
const bool energyImplicit = (energyCouplingMode == "implicit");

// Fluks entalpi species di EEqn (REGISTER 2026-09-22, cacat 2) -- fluid/energySpeciesEnthalpy.H
word energySpeciesMode("legacy");
{
    fvSolution solDict(runTime);
    energySpeciesMode =
        solDict.subDict("PIMPLE").getOrDefault<word>("energySpeciesEnthalpy", "legacy");
    if (energySpeciesMode != "legacy" && energySpeciesMode != "consistent")
    {
        FatalIOErrorInFunction(solDict.subDict("PIMPLE"))
            << "energySpeciesEnthalpy '" << energySpeciesMode
            << "' tidak dikenal. Pilihan: legacy, consistent." << exit(FatalIOError);
    }
    Info<< "energySpeciesEnthalpy " << energySpeciesMode << nl << endl;
}
const bool energySpeciesConsistent = (energySpeciesMode == "consistent");

// Koreksi kontinuitas di EEqn (Paper A_v2, REGISTER_V2 V2-SF1b, 2026-10-25). Bentuk konservatif
// ddt(rho,he) + div(phi,he) menyimpan he x residu kontinuitas lokal (ddt(rho) + div(phi) + Smass); pada acuan
// entalpi jauh dari T, he berbeda besar antar-komposisi sehingga residu itu tidak saling menghapus di front.
// Opsi ini mengurangkan he x residu (bentuk rho Dhe/Dt; nol bila kontinuitas terpenuhi). Bawaan mati.
// Titik linearisasi kopling energi implisit (V2-SF1c, 2026-10-25). BC energyCoupledImplicit memakai
// T = T* + (h - h*)/Cp, eksak hanya bila h* = h(T*, Y) dgn komposisi KINI. h* = he lama dihitung dgn komposisi
// SEBELUM YEqn iterasi ini, sehingga suhu antarmuka menyimpang h*(1/Cp_lama - 1/Cp_baru) ~ offset acuan entalpi.
// Opsi ini memakai h* = thermo.he(p, T*) dgn komposisi kini. Bawaan mati (kasus lama bit-identik).
bool energyLinearisationConsistent = false;
{
    fvSolution solDict(runTime);
    energyLinearisationConsistent =
        solDict.subDict("PIMPLE").getOrDefault<bool>("energyLinearisationConsistent", false);
    Info<< "energyLinearisationConsistent " << energyLinearisationConsistent << nl << endl;
}
bool energyContinuityCorrection = false;
{
    fvSolution solDict(runTime);
    energyContinuityCorrection =
        solDict.subDict("PIMPLE").getOrDefault<bool>("energyContinuityCorrection", false);
    Info<< "energyContinuityCorrection " << energyContinuityCorrection << nl << endl;
}

if (energyImplicit)
{
    if (coupled)
    {
        FatalErrorInFunction
            << "energyInterfaceCoupling implicit tidak boleh digabung dgn jalur"
            << " coupledEnergyField / useImplicit pada BC suhu (bawaan)."
            << exit(FatalError);
    }
    if (fluidRegions.size() != 1 || solidRegions.size() != 1)
    {
        FatalErrorInFunction
            << "energyInterfaceCoupling implicit saat ini hanya untuk SATU region fluida"
            << " dan SATU region padatan (ada " << fluidRegions.size() << " dan "
            << solidRegions.size() << ")." << exit(FatalError);
    }
}

PtrList<volScalarField> hCplFluid(fluidRegions.size());
PtrList<volScalarField> kappaCplFluid(fluidRegions.size());
PtrList<volScalarField> CpCplFluid(fluidRegions.size());
PtrList<volScalarField> hCplSolid(solidRegions.size());
PtrList<volScalarField> kappaCplSolid(solidRegions.size());
PtrList<volScalarField> CpCplSolid(solidRegions.size());

// Persamaan energi region yg menunggu dirakit dlm iterasi luar ini
PtrList<fvScalarMatrix> energyEqns;

if (energyImplicit)
{
    // Tipe patch hCpl: antarmuka -> energyCoupledImplicit; terkopel -> tipe he
    auto tipeCpl = [](const volScalarField& T, const volScalarField& he)
    {
        wordList t(he.boundaryField().size());
        forAll(t, pi)
        {
            if (he.boundaryField()[pi].patch().coupled())
            {
                t[pi] = he.boundaryField()[pi].type();
            }
            else if
            (
                T.boundaryField()[pi].type().find("CoupledMixed")
             != std::string::npos
            )
            {
                t[pi] = "energyCoupledImplicit";
            }
            else
            {
                t[pi] = calculatedFvPatchScalarField::typeName;
            }
        }
        return t;
    };

    auto buatCpl = [&runTime, &tipeCpl]
    (
        const fvMesh& mesh,
        const volScalarField& T,
        const volScalarField& he,
        const tmp<volScalarField>& kappa,
        const tmp<volScalarField>& Cp,
        PtrList<volScalarField>& hL,
        PtrList<volScalarField>& kL,
        PtrList<volScalarField>& cL,
        const label i
    )
    {
        const auto io = [&](const word& n)
        {
            return IOobject
            (
                n, runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, true
            );
        };

        hL.set(i, new volScalarField(io("hCpl"), he, tipeCpl(T, he)));
        kL.set(i, new volScalarField(io("kappaCpl"), kappa()));
        cL.set(i, new volScalarField(io("CpCpl"), Cp()));

        label n = 0;
        forAll(hL[i].boundaryField(), pi)
        {
            auto& pf = hL[i].boundaryFieldRef()[pi];
            if (pf.type() == "energyCoupledImplicit")
            {
                pf.useImplicit(true);
                ++n;
            }
        }
        if (returnReduce(n, sumOp<label>()) == 0)
        {
            FatalErrorInFunction
                << "region " << mesh.name() << ": tidak ada patch antarmuka"
                << " suhu (*CoupledMixed) untuk energyInterfaceCoupling implicit."
                << exit(FatalError);
        }
        return n;
    };

    label nImpl = 0;
    forAll(fluidRegions, i)
    {
        rhoReactionThermo& th = thermoFluid[i];
        nImpl += buatCpl
        (
            fluidRegions[i], th.T(), th.he(),
            turbulenceFluid[i].kappaEff(), th.Cp(),
            hCplFluid, kappaCplFluid, CpCplFluid, i
        );
    }
    forAll(solidRegions, i)
    {
        solidThermo& th = thermos[i];
        nImpl += buatCpl
        (
            solidRegions[i], th.T(), th.he(),
            betavSolid[i]*th.kappa(), th.Cp(),
            hCplSolid, kappaCplSolid, CpCplSolid, i
        );
    }

    Info<< "energyInterfaceCoupling implicit: " << nImpl
        << " patch antarmuka energi dirakit implisit" << nl << endl;
}

// Salin koefisien persamaan energi region ke matriks di atas hCpl, simpan.
// Patch antarmuka hCpl implisit -> perakit mengabaikan internal/boundaryCoeffs
// patch itu (kontribusi partisi BC suhu), lalu energyCoupledImplicit menulis
// kopling penuh. Patch lain (inlet, outlet, processor) terbawa apa adanya.
auto simpanEnergi = [&energyEqns]
(
    const fvScalarMatrix& eqn,
    volScalarField& hCpl
)
{
    fvScalarMatrix M(hCpl, eqn.dimensions());
    M.lduMatrix::operator=(eqn);
    M.source() = eqn.source();
    M.internalCoeffs() = eqn.internalCoeffs();
    M.boundaryCoeffs() = eqn.boundaryCoeffs();
    energyEqns.append(new fvScalarMatrix(M));
};

// Kehilangan panas dinding kolom (REGISTER 2026-09-24, GRua) -- opsional per padatan.
// constant/<padatan>/wallHeatLoss { hv <W/m3/K>; Tambient <K>; }
// Sumber volumetrik padatan  S = -hv (T - Tambient)  (unit-cell mewakili kanal rata2
// kolom yg melepas panas lewat dinding luar ke oven). Tanpa berkas: hv = 0, bit-identik.
List<scalar> hvLoss(solidRegions.size(), Zero);
List<scalar> TambLoss(solidRegions.size(), Zero);
forAll(solidRegions, i)
{
    IOdictionary whl
    (
        IOobject
        (
            "wallHeatLoss", runTime.constant(), solidRegions[i],
            IOobject::READ_IF_PRESENT, IOobject::NO_WRITE, false
        )
    );
    if (whl.headerOk())
    {
        hvLoss[i] = whl.get<scalar>("hv");
        TambLoss[i] = whl.get<scalar>("Tambient");
        if (hvLoss[i] < 0)
        {
            FatalIOErrorInFunction(whl) << "hv harus >= 0" << exit(FatalIOError);
        }
        Info<< "wallHeatLoss " << solidRegions[i].name() << ": hv = " << hvLoss[i]
            << " W/m3/K, Tambient = " << TambLoss[i] << " K" << nl << endl;
    }
}
