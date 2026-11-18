#include "energyCoupledImplicitFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "mappedPatchBase.H"
#include "lduPrimitiveMeshAssembly.H"

namespace Foam
{

energyCoupledImplicitFvPatchScalarField::energyCoupledImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF)
{}


energyCoupledImplicitFvPatchScalarField::energyCoupledImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << p.name() << " bukan mapped (tipe "
            << p.patch().type() << "); energyCoupledImplicit butuh pasangan"
            << " antarmuka antar-region." << exit(FatalIOError);
    }
}


energyCoupledImplicitFvPatchScalarField::energyCoupledImplicitFvPatchScalarField
(
    const energyCoupledImplicitFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


energyCoupledImplicitFvPatchScalarField::energyCoupledImplicitFvPatchScalarField
(
    const energyCoupledImplicitFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf)
{}


energyCoupledImplicitFvPatchScalarField::energyCoupledImplicitFvPatchScalarField
(
    const energyCoupledImplicitFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF)
{}


tmp<scalarField> energyCoupledImplicitFvPatchScalarField::conductance() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const label nbrPatchId = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch = nbrMesh.boundary()[nbrPatchId];

    const scalarField K
    (
        db().lookupObject<volScalarField>("kappaCpl")
           .boundaryField()[patch().index()]
       *patch().deltaCoeffs()
    );
    const scalarField Knbr
    (
        nbrMesh.lookupObject<volScalarField>("kappaCpl")
           .boundaryField()[nbrPatchId]
       *nbrPatch.deltaCoeffs()
    );

    if (Knbr.size() != K.size())
    {
        FatalErrorInFunction
            << "patch " << patch().name() << ": jumlah muka tak cocok (lokal "
            << K.size() << ", tetangga " << Knbr.size() << "). Kopling energi"
            << " implisit menuntut kedua sisi antarmuka di prosesor yang sama."
            << exit(FatalError);
    }

    return K*Knbr/(K + Knbr + VSMALL)*patch().magSf();
}


void energyCoupledImplicitFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const label mat,
    const direction
)
{
    const label index = patch().index();
    const label nbrPatchId = patch().patch().neighbPolyPatchID();

    const lduPrimitiveMeshAssembly& asmb = matrix.lduMeshAssembly();
    const label meshNbrId = asmb.findNbrMeshId(patch().patch(), mat);

    // Sama persis dgn speciesCoupledMixed::manipulateMatrix (terbukti):
    // myCells = sel region INI, nbrCells = sel region tetangga (indeks global).
    const labelList& faceMap = asmb.faceBoundMap()[mat][index];
    const labelList& myCells = asmb.cellBoundMap()[meshNbrId][nbrPatchId];
    const labelList& nbrCells = asmb.cellBoundMap()[mat][index];

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const fvPatch& nbrPatch = nbrMesh.boundary()[nbrPatchId];

    if (faceMap.size() != patch().size() || nbrPatch.size() != patch().size())
    {
        FatalErrorInFunction
            << "patch " << patch().name() << ": peta rakitan " << faceMap.size()
            << ", lokal " << patch().size() << ", tetangga " << nbrPatch.size()
            << ". Kopling energi implisit menuntut kedua sisi di prosesor sama."
            << exit(FatalError);
    }

    const scalarField G(conductance());

    // Keadaan linearisasi sisi ini (P) dan tetangga (Q), nilai sel di muka.
    const scalarField TP
    (
        db().lookupObject<volScalarField>("T").boundaryField()[index]
           .patchInternalField()
    );
    const scalarField CpP
    (
        db().lookupObject<volScalarField>("CpCpl").boundaryField()[index]
           .patchInternalField()
    );
    const scalarField hP(this->patchInternalField());

    const volScalarField& hNbrFld =
        nbrMesh.lookupObject<volScalarField>(this->internalField().name());
    const scalarField TQ
    (
        nbrMesh.lookupObject<volScalarField>("T").boundaryField()[nbrPatchId]
           .patchInternalField()
    );
    const scalarField CpQ
    (
        nbrMesh.lookupObject<volScalarField>("CpCpl").boundaryField()[nbrPatchId]
           .patchInternalField()
    );
    const scalarField hQ(hNbrFld.boundaryField()[nbrPatchId].patchInternalField());

    // Alokasi lower SEBELUM menulis upper (lihat speciesCoupledMixed).
    matrix.lower();

    forAll(faceMap, j)
    {
        const label facei = faceMap[j];
        const label myCell = myCells[j];
        const label nbrCell = nbrCells[j];

        const scalar cP = max(CpP[j], SMALL);
        const scalar cQ = max(CpQ[j], SMALL);

        matrix.diag()[myCell] += G[j]/cP;

        if (myCell < nbrCell)
        {
            matrix.upper()[facei] -= G[j]/cQ;
        }
        else
        {
            matrix.lower()[facei] -= G[j]/cQ;
        }

        matrix.source()[myCell] +=
            G[j]*((TQ[j] - hQ[j]/cQ) - (TP[j] - hP[j]/cP));
    }
}


makePatchTypeField
(
    fvPatchScalarField,
    energyCoupledImplicitFvPatchScalarField
);

} // End namespace Foam
