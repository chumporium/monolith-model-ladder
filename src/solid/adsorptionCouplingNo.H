/*---------------------------------------------------------------------------*\
  adsorptionCouplingNo.H

  BILANGAN KOPLING  Lambda  -- pengendali dt ketiga, sejajar Co dan Di.

      Lambda = (eps_p Ds / effStor) * dt * deltaCoeffs^2
             = Ds dt / (R dx^2),      R = effStor/eps_p = retardasi

  KENAPA INI PERLU
  ----------------
  Persamaan species di padatan berbentuk

      effStor dC/dt = div(eps_p Ds grad C)

  jadi difusivitas EFEKTIF untuk C adalah eps_p Ds / effStor, BUKAN Ds
  telanjang. Kopling antar-region diselesaikan EKSPLISIT (tiap region memakai
  nilai tetangga yang tertinggal satu iterasi), dan iterasi titik-tetap itu
  hanya kontraktif selama Lambda cukup kecil.

  Diukur pada kampanye C56 5,6% (doc/tahap1 4h): faktor kontraksi per iterasi
  luar memburuk sistematis dengan dt --

      dt 0,0012 -> 0,62 | 0,0074 -> 0,61 | 0,0089 -> 0,73 | 0,0185 -> 0,92 gagal

  dan transisinya terjepit rapat antara Lambda 0,272 (lulus) dan 0,327 (gagal).
  Kampanye lapangan berjalan pada Lambda = 9,3, yaitu 31x di atas ambang, dan
  kehilangan 40-50% massanya tanpa satu pun peringatan.

  KENAPA Co DAN Di TIDAK MENANGKAPNYA
  -----------------------------------
  maxCo  : Courant konveksi FLUIDA -- tak menyentuh padatan sama sekali.
  maxDi  : solidRegionDiffNo memakai kappa/(rho Cp), yaitu difusi PANAS.
           Difusi SPECIES di padatan tak pernah punya kendali dt apa pun,
           dan retardasi R tak pernah masuk hitungan mana pun.

  R itulah yang membedakan rezim: 3649 pada 400 ppm vs 49 pada 5,6% -- 74x.
  Pada 400 ppm isoterm curam sehingga padatan menyerap gangguan tiap iterasi
  (Lambda = 0,013, 23x DI BAWAH ambang); pada 5,6% isoterm sudah datar, R
  runtuh, dan nilai antarmuka berayun.

  Bawaan `maxLambda` = GREAT -> tak berpengaruh, seluruh case lama bit-identik.
\*---------------------------------------------------------------------------*/

scalar LambdaNum = -GREAT;

forAll(solidRegions, i)
{
    const fvMesh& smesh = solidRegions[i];
    const regionAdsorption& ads = adsorption[i];
    const scalar epsP = ads.epsilonP();

    forAll(ads.species(), j)
    {
        const volScalarField& Dsj = ads.Ds(j);
        const volScalarField& eSj = ads.effStor(j);

        // Difusivitas efektif untuk C, di MUKA sel (konsisten dengan
        // laplacian yang benar-benar didiskretkan).
        surfaceScalarField DeffByDelta
        (
            sqr(smesh.surfaceInterpolation::deltaCoeffs())
           *fvc::interpolate(epsP*Dsj)
           /fvc::interpolate(eSj)
        );

        LambdaNum = max
        (
            gMax(DeffByDelta.primitiveField())*runTime.deltaTValue(),
            LambdaNum
        );
    }
}

if (LambdaNum > -GREAT)
{
    Info<< "Region: adsorpsi  bilangan kopling Lambda max: "
        << LambdaNum << endl;
}
