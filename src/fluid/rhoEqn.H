/*---------------------------------------------------------------------------*\
  rhoEqn.H  (lokal — menggantikan rhoEqn.H standar OpenFOAM via -I./fluid)

  Kontinuitas dengan sink massa adsorpsi (kopling dua-arah, Bab 5.1):

      d(rho)/dt + div(phi) = fvOptions(rho) - Smass

  Smass > 0 (adsorpsi) mengurangi massa gas; Smass < 0 (desorpsi/TVSA)
  menambah. SmassFluid[i] dihitung di fluidAdsorptionSource.H. Bila
  twoWayCoupling nonaktif, SmassFluid tetap 0.
\*---------------------------------------------------------------------------*/

{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
     ==
        fvOptions(rho)
      - SmassFluid[i]
    );

    rhoEqn.relax();

    rhoEqn.solve();

    fvOptions.correct(rho);
}
