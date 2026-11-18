/*---------------------------------------------------------------------------*\
  createSpeciesCoupling.H  -- sekali saat start

  PIMPLE { speciesCoupling partitioned | implicit; }   (bawaan partitioned)

  partitioned : jalur lama (BC mixed, iterasi luar), bit-identik.
  implicit    : Y adsorbat fluida + C padatan dirakit satu matriks tiap
                iterasi Newton padatan (speciesImplicitCouplerFluid.H,
                speciesCoupledMixed::manipulateMatrix). Konvergensi tidak lagi
                dibatasi bilangan kopling Lambda (adsorptionCouplingNo.H).
\*---------------------------------------------------------------------------*/

word speciesCouplingMode("partitioned");
{
    fvSolution solDict(runTime);
    speciesCouplingMode =
        solDict.subDict("PIMPLE").getOrDefault<word>("speciesCoupling", "partitioned");

    if (speciesCouplingMode != "partitioned" && speciesCouplingMode != "implicit")
    {
        FatalIOErrorInFunction(solDict.subDict("PIMPLE"))
            << "speciesCoupling '" << speciesCouplingMode
            << "' tidak dikenal. Pilihan: partitioned, implicit."
            << exit(FatalIOError);
    }

    if (speciesCouplingMode == "implicit")
    {
        if
        (
            solDict.subDict("PIMPLE").getOrDefault<word>
            (
                "interfaceRelaxation", "none"
            ) != "none"
        )
        {
            FatalIOErrorInFunction(solDict.subDict("PIMPLE"))
                << "speciesCoupling implicit tidak memerlukan (dan tidak boleh"
                << " digabung dengan) interfaceRelaxation." << exit(FatalIOError);
        }

        if (fluidRegions.size() != 1 || solidRegions.size() != 1)
        {
            FatalErrorInFunction
                << "speciesCoupling implicit saat ini hanya untuk SATU region"
                << " fluida dan SATU region padatan (ada "
                << fluidRegions.size() << " dan " << solidRegions.size() << ")."
                << exit(FatalError);
        }
    }
}

autoPtr<speciesImplicitCouplerFluid> speciesCouplerPtr;

if (speciesCouplingMode == "implicit")
{
    // useImplicit pada patch speciesCoupledMixed kedua sisi
    label nImpl = 0;
    forAll(fluidRegions, ri)
    {
        PtrList<volScalarField>& Yf = thermoFluid[ri].composition().Y();
        forAll(Yf, j)
        {
            forAll(Yf[j].boundaryField(), patchi)
            {
                auto& pf = Yf[j].boundaryFieldRef()[patchi];
                if (pf.type() == "speciesCoupledMixed")
                {
                    pf.useImplicit(true);
                    ++nImpl;
                }
            }
        }
    }
    forAll(solidRegions, ri)
    {
        const regionAdsorption& ads = adsorption[ri];
        forAll(ads.species(), j)
        {
            volScalarField& Cj = const_cast<volScalarField&>(ads.C(j));
            forAll(Cj.boundaryField(), patchi)
            {
                auto& pf = Cj.boundaryFieldRef()[patchi];
                if (pf.type() == "speciesCoupledMixed")
                {
                    pf.useImplicit(true);
                    ++nImpl;
                }
            }
        }
    }

    speciesCouplerPtr.reset
    (
        new speciesImplicitCouplerFluid
        (
            reactionFluid, rhoFluid, phiFluid, turbulenceFluid, fieldsFluid,
            fluidFvOptions, DmolList, SctFluid, 0, adsorption[0].species()
        )
    );

    speciesCouplerPtr->setAdvSource(&SadvFluid);   // V2: sumber adveksi gas pori ikut dirakit

    Info<< "speciesCoupling implicit: " << nImpl
        << " patch speciesCoupledMixed dirakit implisit" << nl << endl;
}
