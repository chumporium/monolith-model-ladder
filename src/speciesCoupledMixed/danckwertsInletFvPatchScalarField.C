/*---------------------------------------------------------------------------*\
  danckwertsInletFvPatchScalarField.C   -- lihat .H
\*---------------------------------------------------------------------------*/

#include "danckwertsInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

danckwertsInletFvPatchScalarField::danckwertsInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    feedValue_(0),
    D_(0),
    rhoName_("rho"),
    phiName_("phi")
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 1;
}


danckwertsInletFvPatchScalarField::danckwertsInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    feedValue_(dict.get<scalar>("feedValue")),
    D_(dict.get<scalar>("D")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
{
    if (D_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "D harus > 0 (difusivitas molekular species) pada patch "
            << p.name() << exit(FatalIOError);
    }

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


danckwertsInletFvPatchScalarField::danckwertsInletFvPatchScalarField
(
    const danckwertsInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    feedValue_(ptf.feedValue_),
    D_(ptf.D_),
    rhoName_(ptf.rhoName_),
    phiName_(ptf.phiName_)
{}


danckwertsInletFvPatchScalarField::danckwertsInletFvPatchScalarField
(
    const danckwertsInletFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    feedValue_(ptf.feedValue_),
    D_(ptf.D_),
    rhoName_(ptf.rhoName_),
    phiName_(ptf.phiName_)
{}


danckwertsInletFvPatchScalarField::danckwertsInletFvPatchScalarField
(
    const danckwertsInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    feedValue_(ptf.feedValue_),
    D_(ptf.D_),
    rhoName_(ptf.rhoName_),
    phiName_(ptf.phiName_)
{}


void danckwertsInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto* phiPtr = this->db().findObject<surfaceScalarField>(phiName_);
    const auto* rhoPtr = this->db().findObject<volScalarField>(rhoName_);

    if (!phiPtr || !rhoPtr)
    {
        FatalErrorInFunction
            << "patch " << patch().name() << ": field '" << phiName_ << "' atau '"
            << rhoName_ << "' tidak ditemukan" << exit(FatalError);
    }

    const scalarField& phip = phiPtr->boundaryField()[patch().index()];
    const scalarField& rhop = rhoPtr->boundaryField()[patch().index()];

    // Konduktansi difusif muka [kg/s]: rho D delta |Sf|
    const scalarField G(rhop*D_*patch().deltaCoeffs()*patch().magSf());
    const scalarField inflow(max(-phip, scalar(0)));

    this->refValue() = feedValue_;
    this->refGrad() = 0;
    this->valueFraction() = inflow/(inflow + G + VSMALL);

    mixedFvPatchScalarField::updateCoeffs();
}


void danckwertsInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("feedValue", feedValue_);
    os.writeEntry("D", D_);
    os.writeEntry("rho", rhoName_);
    os.writeEntry("phi", phiName_);
    this->writeEntry("value", os);
}


makePatchTypeField
(
    fvPatchScalarField,
    danckwertsInletFvPatchScalarField
);

} // End namespace Foam
