/*---------------------------------------------------------------------------*\
  interfaceAitkenInit.H  -- awal tiap langkah waktu (sebelum loop PIMPLE)

  Membuat (sekali) dan mengisi field x = "<C>Aitken" dan residu "<C>AitkenR"
  di region padatan. x_0 = C sel antarmuka saat ini, yaitu persis nilai yang
  dilihat fluida tanpa relaksasi -> iterasi luar pertama tidak berubah.
  Lihat interfaceAitken.H.
\*---------------------------------------------------------------------------*/

if (interfaceAitken)
{
    forAll(solidRegions, ri)
    {
        fvMesh& smesh = solidRegions[ri];
        const regionAdsorption& ads = adsorption[ri];

        forAll(ads.species(), j)
        {
            const volScalarField& Cj = ads.C(j);

            for (const word& suffix : {word("Aitken"), word("AitkenR")})
            {
                const word nm(Cj.name() + suffix);
                if (!smesh.foundObject<volScalarField>(nm))
                {
                    auto* fPtr = new volScalarField
                    (
                        IOobject
                        (
                            nm,
                            runTime.timeName(),
                            smesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            IOobject::REGISTER
                        ),
                        smesh,
                        dimensionedScalar(Cj.dimensions(), Zero),
                        calculatedFvPatchScalarField::typeName
                    );
                    fPtr->store();
                }
            }

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
                xF.boundaryFieldRef()[patchi] =
                    Cj.boundaryField()[patchi].patchInternalField();
                rF.boundaryFieldRef()[patchi] = Zero;
            }
        }
    }
}
