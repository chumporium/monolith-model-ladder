/*---------------------------------------------------------------------------*\
  solveSolidSpecies.H  (pluggable adsorption version)

  Menyelesaikan persamaan species pori + isoterm + closure kinetik melalui
  koordinator regionAdsorption, dengan sub-iterasi q*<->T (nSolidCouple).

  Semua fisika (isoterm, aturan campuran, closure, mesin konversi, Bosanquet
  Ds(P), source energi S_ads) ditangani di library adsorptionModels.
\*---------------------------------------------------------------------------*/
{
    if (speciesCouplerPtr)
    {
        speciesCouplerPtr->setFinalIter(finalIter);
    }
    adsorptionI.correct(nSolidCouple, speciesCouplerPtr.get());
}
