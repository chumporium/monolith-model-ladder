/*---------------------------------------------------------------------------*\
  danckwertsInletTemperatureFvPatchScalarField.C   -- lihat .H
\*---------------------------------------------------------------------------*/

#include "danckwertsInletTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basicThermo.H"

namespace Foam
{

danckwertsInletTemperatureFvPatchScalarField::danckwertsInletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    feedValue_(0),
    phiName_("phi"),
    explicitBlend_(false)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 1;
}


danckwertsInletTemperatureFvPatchScalarField::danckwertsInletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    feedValue_(dict.get<scalar>("feedValue")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    explicitBlend_(dict.getOrDefault<bool>("explicitBlend", false))
{
    this->refValue() = feedValue_;
    this->refGrad() = 0;
    this->valueFraction() = 1;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(feedValue_);
    }
}


danckwertsInletTemperatureFvPatchScalarField::danckwertsInletTemperatureFvPatchScalarField
(
    const danckwertsInletTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    feedValue_(ptf.feedValue_),
    phiName_(ptf.phiName_),
    explicitBlend_(ptf.explicitBlend_)
{}


danckwertsInletTemperatureFvPatchScalarField::danckwertsInletTemperatureFvPatchScalarField
(
    const danckwertsInletTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    feedValue_(ptf.feedValue_),
    phiName_(ptf.phiName_),
    explicitBlend_(ptf.explicitBlend_)
{}


danckwertsInletTemperatureFvPatchScalarField::danckwertsInletTemperatureFvPatchScalarField
(
    const danckwertsInletTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    feedValue_(ptf.feedValue_),
    phiName_(ptf.phiName_),
    explicitBlend_(ptf.explicitBlend_)
{}


void danckwertsInletTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto* phiPtr = this->db().findObject<surfaceScalarField>(phiName_);
    if (!phiPtr)
    {
        FatalErrorInFunction
            << "patch " << patch().name() << ": field '" << phiName_
            << "' tidak ditemukan" << exit(FatalError);
    }

    const label pi = patch().index();
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const scalarField& phip = phiPtr->boundaryField()[pi];
    const scalarField& pp = thermo.p().boundaryField()[pi];

    // Kapasitas konveksi [W/K] dan konduktansi muka [W/K]
    const scalarField Cpb(thermo.Cp(pp, *this, pi));
    const scalarField inflow(max(-phip, scalar(0))*Cpb);
    const scalarField G(thermo.kappa(pi)*patch().deltaCoeffs()*patch().magSf());

    if (explicitBlend_)
    {
        const scalarField f(inflow/(inflow + G + VSMALL));
        this->refValue() = f*feedValue_ + (1.0 - f)*this->patchInternalField();
        this->refGrad() = 0;
        this->valueFraction() = 1;
    }
    else
    {
        this->refValue() = feedValue_;
        this->refGrad() = 0;
        this->valueFraction() = inflow/(inflow + G + VSMALL);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void danckwertsInletTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("feedValue", feedValue_);
    os.writeEntry("phi", phiName_);
    if (explicitBlend_) os.writeEntry("explicitBlend", explicitBlend_);
    this->writeEntry("value", os);
}


makePatchTypeField
(
    fvPatchScalarField,
    danckwertsInletTemperatureFvPatchScalarField
);

} // End namespace Foam
