/*---------------------------------------------------------------------------*\
  energySpeciesEnthalpy.H  -- dipanggil dari EEqn.H (REGISTER 2026-09-22, cacat 2)

  PIMPLE { energySpeciesEnthalpy legacy | consistent; }   (bawaan legacy)

  EEqn bawaan mendifusikan ENTALPI CAMPURAN dgn alphaEff:  -div(alpha grad h),
  padahal  alpha grad h = kappa grad T + alpha SUM h_i grad Y_i  (hConst).
  Suku kedua mengangkut entalpi species dgn alpha, bukan dgn fluks massa
  species J_i = -DEff_i grad Y_i. Pada Le = alpha/(rho D) ~ 2,3 dan h_i
  relatif Tstd 298 K (h_He - h_C3H6 ~ 400 kJ/kg pd 423 K) suku semu itu
  sebesar panas adsorpsi (EN/ENr1: padatan MENDINGIN 1,3 K). Di dinding,
  BC energi membawa selisih he akibat Y dinding != Y sel yg tak dibayar
  padatan. Dan massa yg diserap keluar dari kontinuitas (Smass) tanpa
  membawa entalpinya dari EEqn. Ketiganya bergantung pada acuan entalpi.

  consistent:
    fluks difusi energi fisis (koreksi tertunda thd laplacian implisit):
        q = -kappaEff grad T + SUM_i h_i J_i ,   J_inert = -SUM_aktif J_i
        dq_f = [alphaEff snGrad(he) - kappaEff snGrad(T)]|Sf| + SUM_i h_i,f J_i,f
        EEqn += -div(dq)       (saat konvergen fvm/fvc laplacian saling hapus)
    patch antarmuka: tanpa suku species (massa adsorbat ditangani sink);
        partisi: dq = alphaEff snGrad(he) - kappaEff snGrad(T) (hapus selisih
        komposisi BC dinding); implisit energi: dq = 0 (energyCoupledImplicit).
    sink dinding: massa adsorbat membawa entalpinya sendiri,
        EEqn += -mdot_i h_i(T_P)/V_P   (twoWayCoupling; tanpa itu Y adsorbat
        ditukar ke inert: -mdot_i (h_i - h_inert)).
  Hasilnya: suhu tak bergantung pada acuan entalpi (Tref/Tstd).

  Diabaikan (dicatat): kapasitas panas fasa teradsorpsi (~2,6% kapasitas
  padatan pd GR); entalpi adsorbat yg keluar dari gas tak masuk inventori
  mana pun di neraca berbasis T -> sistematik ~ -0,4% errE pd GR.
\*---------------------------------------------------------------------------*/

volScalarField SentF
(
    IOobject("SentF", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

if (energySpeciesConsistent && Y.size())
{
    const volScalarField& T = thermo.T();
    const volScalarField& pF = thermo.p();
    const volScalarField& heF = thermo.he();
    const tmp<volScalarField> tAlpE(turbulence.alphaEff());
    const tmp<volScalarField> tKapE(turbulence.kappaEff());
    const basicSpecieMixture& compE = thermo.composition();

    // Entalpi sensibel species k pd (p, T) sel & batas
    auto hSpecies = [&](const label k)
    {
        tmp<volScalarField> th
        (
            new volScalarField
            (
                IOobject("hs." + Y[k].name(), runTime.timeName(), mesh),
                mesh,
                dimensionedScalar(dimEnergy/dimMass, Zero)
            )
        );
        volScalarField& hk = th.ref();
        forAll(hk, c)
        {
            hk[c] = compE.HE(k, pF[c], T[c]);
        }
        forAll(hk.boundaryField(), pi)
        {
            if (hk.boundaryField()[pi].coupled()) continue;
            fvPatchScalarField& hp = hk.boundaryFieldRef()[pi];
            const fvPatchScalarField& Tp = T.boundaryField()[pi];
            const fvPatchScalarField& pp = pF.boundaryField()[pi];
            forAll(hp, f)
            {
                hp[f] = compE.HE(k, pp[f], Tp[f]);
            }
        }
        hk.correctBoundaryConditions();     // patch processor
        return th;
    };

    surfaceScalarField dq
    (
        "dqEnergySpecies",
        fvc::interpolate(tAlpE())*fvc::snGrad(heF)*mesh.magSf()
      - fvc::interpolate(tKapE())*fvc::snGrad(T)*mesh.magSf()
    );

    surfaceScalarField Jaktif
    (
        IOobject("JaktifEnergy", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimMass/dimTime, Zero)
    );

    forAll(Y, j)
    {
        if (j == inertIndex || !composition.active(j)) continue;

        const dimensionedScalar DmolJ("Dmol", dimArea/dimTime, DmolList[i][j]);
        const volScalarField DEffJ(rho*DmolJ + turbulence.mut()/SctFluid);
        const surfaceScalarField Jj
        (
            -fvc::interpolate(DEffJ)*fvc::snGrad(Y[j])*mesh.magSf()
        );
        dq += fvc::interpolate(hSpecies(j)())*Jj;
        Jaktif += Jj;
    }
    if (inertIndex >= 0)
    {
        dq -= fvc::interpolate(hSpecies(inertIndex)())*Jaktif;
    }

    // Patch antarmuka suhu: tanpa suku species; implisit -> 0
    forAll(dq.boundaryField(), pi)
    {
        if (T.boundaryField()[pi].type().find("CoupledMixed") == std::string::npos)
        {
            continue;
        }
        fvsPatchScalarField& dqp = dq.boundaryFieldRef()[pi];
        if (energyImplicit)
        {
            dqp = 0;
        }
        else
        {
            const fvPatch& pp = mesh.boundary()[pi];
            dqp =
                tAlpE().boundaryField()[pi]*heF.boundaryField()[pi].snGrad()*pp.magSf()
              - tKapE().boundaryField()[pi]*T.boundaryField()[pi].snGrad()*pp.magSf();
        }
    }

    SentF = -fvc::div(dq);

    // Sink dinding: adsorbat keluar dari gas membawa entalpinya sendiri
    scalarField& S = SentF.primitiveFieldRef();
    const scalarField& V = mesh.V();
    forAll(Y, j)
    {
        if (j == inertIndex || !composition.active(j)) continue;

        forAll(Y[j].boundaryField(), pi)
        {
            const fvPatchScalarField& Yp = Y[j].boundaryField()[pi];
            if (Yp.type() != "speciesCoupledMixed") continue;

            const fvPatch& pp = mesh.boundary()[pi];
            const labelUList& fc = pp.faceCells();
            const scalarField rhop(rho.boundaryField()[pi]);
            const scalarField snG((Yp - Yp.patchInternalField())*pp.deltaCoeffs());
            const scalarField mdot(-rhop*DmolList[i][j]*snG*pp.magSf());  // >0 terserap

            forAll(fc, f)
            {
                const label c = fc[f];
                scalar hPer = compE.HE(j, pF[c], T[c]);
                if (!twoWayCoupling && inertIndex >= 0)
                {
                    hPer -= compE.HE(inertIndex, pF[c], T[c]);
                }
                S[c] -= mdot[f]*hPer/V[c];
            }
        }
    }
}
