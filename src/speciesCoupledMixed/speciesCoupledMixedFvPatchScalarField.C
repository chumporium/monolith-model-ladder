/*---------------------------------------------------------------------------*\
  speciesCoupledMixedFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#include "speciesCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "IOdictionary.H"
#include "fvMatrix.H"
#include "lduPrimitiveMeshAssembly.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

// Sisi solid: eps_p diambil dari SATU sumber -- constant/<region>/
// adsorptionProperties. Lihat CATATAN 23.80[1]/23.82: sebelum ini eps_p hidup
// di DUA berkas (dict + entri `D` BC ini) tanpa pemeriksaan silang, sehingga
// mengubah material di satu tempat saja memberi persamaan pori dan fluks
// antarmuka yang memakai eps_p BERBEDA -- ketidakcocokan yang hanya terlihat
// setelah satu run, lewat neraca massa antarmuka.
Foam::scalar speciesCoupledMixedFvPatchScalarField::readEpsilonP
(
    const fvPatch& p,
    const dictionary& dict
)
{
    const fvMesh& mesh = p.boundaryMesh().mesh();

    IOdictionary adsDict
    (
        IOobject
        (
            "adsorptionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const scalar epsP = adsDict.get<scalar>("epsilonP");

    // Entri `D` pada sisi solid sekarang USANG. Ia dibiarkan lolos hanya bila
    // nilainya masih sama -- supaya case lama tetap berjalan identik (provenans
    // 23.38) -- dan MENGGAGALKAN run bila sudah menyimpang.
    if (dict.found("D"))
    {
        const scalar Dold = dict.get<scalar>("D");

        if (mag(Dold - epsP) > 1e-6*max(mag(epsP), SMALL))
        {
            // Berkas & barisnya dicetak sendiri oleh FatalIOError(dict).
            FatalIOErrorInFunction(dict)
                << "eps_p TIDAK COCOK antara dua berkas pada patch "
                << p.name() << ':' << nl
                << "    entri `D` (USANG) di BC ini      : " << Dold << nl
                << "    `epsilonP` di adsorptionProperties: " << epsP << nl
                << nl
                << "Sejak eps_p bersumber tunggal, entri `D` sisi solid tidak"
                << " lagi dipakai." << nl
                << "PERBAIKAN: HAPUS entri `D` dari patch ini. Nilai yang"
                << " berlaku adalah `epsilonP`." << nl
                << "Kalau dibiarkan, persamaan pori dan fluks antarmuka akan"
                << " memakai eps_p berbeda." << nl
                << exit(FatalIOError);
        }

        WarningInFunction
            << "patch " << p.name() << ": entri `D` pada sisi solid USANG dan"
            << " diabaikan (eps_p = " << epsP << " diambil dari"
            << " adsorptionProperties). Hapus entri itu dari dict." << endl;
    }

    return epsP;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

speciesCoupledMixedFvPatchScalarField::speciesCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phase_("fluid"),
    nbrFieldName_("undefined"),
    D_(0),
    DFieldName_("none"),
    rhoName_("rho")
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 1;
}


speciesCoupledMixedFvPatchScalarField::speciesCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phase_(dict.getOrDefault<word>("phase", "fluid")),
    nbrFieldName_(dict.get<word>("nbrField")),
    // fluid : D = D molekular species, tetap dari dict ini.
    // solid : D = eps_p, BERSUMBER TUNGGAL dari adsorptionProperties.
    //         `phase` dibaca ulang dari dict (bukan phase_) supaya nilainya
    //         tidak bergantung urutan inisialisasi anggota.
    D_
    (
        dict.getOrDefault<word>("phase", "fluid") == "solid"
      ? readEpsilonP(p, dict)
      : dict.get<scalar>("D")
    ),
    DFieldName_(dict.getOrDefault<word>("DField", "none")),
    rhoName_(dict.getOrDefault<word>("rho", "rho"))
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "patch " << p.name()
            << " is not of mappedPatchBase type" << nl
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        this->refValue() = scalarField("refValue", dict, p.size());
        this->refGrad() = scalarField("refGradient", dict, p.size());
        this->valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        this->refValue() = *this;
        this->refGrad() = 0;
        this->valueFraction() = 1;
    }
}


speciesCoupledMixedFvPatchScalarField::speciesCoupledMixedFvPatchScalarField
(
    const speciesCoupledMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phase_(ptf.phase_),
    nbrFieldName_(ptf.nbrFieldName_),
    D_(ptf.D_),
    DFieldName_(ptf.DFieldName_),
    rhoName_(ptf.rhoName_)
{}


speciesCoupledMixedFvPatchScalarField::speciesCoupledMixedFvPatchScalarField
(
    const speciesCoupledMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phase_(ptf.phase_),
    nbrFieldName_(ptf.nbrFieldName_),
    D_(ptf.D_),
    DFieldName_(ptf.DFieldName_),
    rhoName_(ptf.rhoName_)
{}


speciesCoupledMixedFvPatchScalarField::speciesCoupledMixedFvPatchScalarField
(
    const speciesCoupledMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phase_(ptf.phase_),
    nbrFieldName_(ptf.nbrFieldName_),
    D_(ptf.D_),
    DFieldName_(ptf.DFieldName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<scalarField> speciesCoupledMixedFvPatchScalarField::effectiveD() const
{
    auto tD = tmp<scalarField>::New(patch().size(), D_);

    if (DFieldName_ != "none")
    {
        const auto* DfPtr =
            this->db().findObject<volScalarField>(DFieldName_);
        if (DfPtr)
        {
            tD.ref() *= DfPtr->boundaryField()[patch().index()];
        }
    }

    return tD;
}


void speciesCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);
    const fvPatch& nbrPatch = nbrFvMesh.boundary()[samplePatchi];

    const volScalarField& nbrFld =
        nbrFvMesh.lookupObject<volScalarField>(nbrFieldName_);

    const speciesCoupledMixedFvPatchScalarField& nbrField =
        refCast<const speciesCoupledMixedFvPatchScalarField>
        (
            nbrFld.boundaryField()[samplePatchi]
        );

    // Nilai internal tetangga (dalam variabel tetangga: Y atau C).
    // Sisi fluida + relaksasi Aitken aktif: pakai nilai C antarmuka yang
    // DIRELAKSASI ("<C>Aitken", dikelola solver -- interfaceAitken.H).
    // Tanpa field itu jalur lama tak berubah.
    const auto* aitkenPtr =
        (phase_ == "fluid")
      ? nbrFvMesh.findObject<volScalarField>(nbrFieldName_ + "Aitken")
      : nullptr;

    scalarField nbrIntFld
    (
        aitkenPtr
      ? scalarField(aitkenPtr->boundaryField()[samplePatchi])
      : nbrField.patchInternalField()()
    );
    mpp.distribute(nbrIntFld);

    // Konduktansi ruang konsentrasi
    scalarField nbrKDelta(nbrField.effectiveD()*nbrPatch.deltaCoeffs());
    mpp.distribute(nbrKDelta);

    scalarField KDelta(effectiveD()*patch().deltaCoeffs());

    // --- Konversi nilai tetangga ke variabel sisi ini via rho ---
    scalarField refConverted(nbrIntFld);

    if (phase_ == "fluid")
    {
        // field ini = Y; tetangga = C  ->  Y = C / rho_lokal
        const auto* rhoPtr = this->db().findObject<volScalarField>(rhoName_);
        if (rhoPtr)
        {
            const scalarField rhoLoc(rhoPtr->boundaryField()[patch().index()]);
            forAll(refConverted, f)
            {
                refConverted[f] = nbrIntFld[f]/max(rhoLoc[f], SMALL);
            }
        }
    }
    else // phase_ == "solid"
    {
        // field ini = C; tetangga = Y  ->  C = Y * rho_fluid(tetangga)
        const auto* rhoPtr =
            nbrFvMesh.findObject<volScalarField>(rhoName_);
        if (rhoPtr)
        {
            scalarField rhoNbr(rhoPtr->boundaryField()[samplePatchi]);
            mpp.distribute(rhoNbr);
            forAll(refConverted, f)
            {
                refConverted[f] = nbrIntFld[f]*rhoNbr[f];
            }
        }
    }

    this->valueFraction() = nbrKDelta/(nbrKDelta + KDelta + VSMALL);
    this->refValue() = refConverted;
    this->refGrad() = 0;

    mixedFvPatchScalarField::updateCoeffs();

    UPstream::msgType() = oldTag;
}


// KOPLING IMPLISIT (2026-09-17)
// ----------------------------
// Kopling partisi (updateCoeffs di atas) hanya konvergen untuk bilangan kopling
// Lambda kecil; zeolit dgn difusi kristal lambat + pori cepat memberi Lambda
// 40-100 (Grande 2006) dan iterasi luar tak konvergen, juga dgn Aitken. Dgn
// `useImplicit true` di kedua sisi, Y fluida dan C padatan dirakit ke SATU
// matriks (fvMatrix::addFvMatrix, mekanisme yg sama dgn kopling energi
// chtMultiRegionFoam) sehingga fluks antarmuka implisit penuh.
//
// Muka antarmuka dgn konduktansi seri (ruang konsentrasi, per luas):
//     Kf = D_f deltaCoeffs_f,  Ks = eps_p Ds deltaCoeffs_s,  G = Kf Ks/(Kf+Ks)
//     C_muka = rho_b Y_muka  ->  J (padatan -> fluida) = G |Sf| (C_Q - rho_b Y_P)
// Identik dgn nilai muka mixed-BC partisi saat konvergen (valueFraction =
// Knbr/(Knbr+K)), jadi hasil konvergen kedua skema sama; bedanya hanya cara
// sampai ke sana. Baris fluida P: diag += G|Sf| rho_b, off(P,Q) = -G|Sf|.
// Baris padatan Q: diag += G|Sf|, off(Q,P) = -G|Sf| rho_b. Asimetris (rho_b),
// jadi lower() selalu dialokasikan lebih dulu.
//
// Perakit mengabaikan internalCoeffs/boundaryCoeffs patch implisit
// (fvMatrix::addBoundaryDiag/Source melewati globalPatchID == -1), maka
// seluruh kopling HARUS ditulis di sini. Pasangan muka dianggap kolokasi
// berindeks sama (mappedWallPolyPatch::mapCollocatedFaces), dan perakit
// menuntut kedua sisi di prosesor yang sama.
void speciesCoupledMixedFvPatchScalarField::manipulateMatrix
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

    const labelList& faceMap = asmb.faceBoundMap()[mat][index];
    const labelList& myCells = asmb.cellBoundMap()[meshNbrId][nbrPatchId];
    const labelList& nbrCells = asmb.cellBoundMap()[mat][index];

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const fvPatch& nbrPatch = nbrMesh.boundary()[nbrPatchId];

    const speciesCoupledMixedFvPatchScalarField& nbrField =
        refCast<const speciesCoupledMixedFvPatchScalarField>
        (
            nbrMesh.lookupObject<volScalarField>(nbrFieldName_)
           .boundaryField()[nbrPatchId]
        );

    if (nbrPatch.size() != patch().size() || faceMap.size() != patch().size())
    {
        FatalErrorInFunction
            << "patch " << patch().name() << ": jumlah muka tak cocok (lokal "
            << patch().size() << ", tetangga " << nbrPatch.size()
            << ", peta rakitan " << faceMap.size() << ")." << nl
            << "Kopling implisit menuntut kedua sisi antarmuka di prosesor"
            << " yang sama." << exit(FatalError);
    }

    const scalarField K(effectiveD()*patch().deltaCoeffs());
    const scalarField Knbr(nbrField.effectiveD()*nbrPatch.deltaCoeffs());
    const scalarField& magSf = patch().magSf();

    const scalarField rhoB
    (
        phase_ == "fluid"
      ? db().lookupObject<volScalarField>(rhoName_).boundaryField()[index]
      : nbrMesh.lookupObject<volScalarField>(rhoName_)
           .boundaryField()[nbrPatchId]
    );

    // Alokasi lower SEBELUM menulis upper: bila matriks masih simetris,
    // lower() menyalin upper saat pertama dipanggil.
    matrix.lower();

    forAll(faceMap, j)
    {
        const label facei = faceMap[j];
        const label myCell = myCells[j];
        const label nbrCell = nbrCells[j];

        const scalar G = K[j]*Knbr[j]/(K[j] + Knbr[j] + VSMALL)*magSf[j];

        const scalar dCoef = (phase_ == "fluid") ? G*rhoB[j] : G;
        const scalar oCoef = (phase_ == "fluid") ? -G : -G*rhoB[j];

        matrix.diag()[myCell] += dCoef;

        if (myCell < nbrCell)
        {
            matrix.upper()[facei] += oCoef;
        }
        else
        {
            matrix.lower()[facei] += oCoef;
        }
    }
}


void speciesCoupledMixedFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("phase", phase_);
    os.writeEntry("nbrField", nbrFieldName_);

    // Sisi solid: `D` (= eps_p) SENGAJA TIDAK ditulis. Kalau ia ditulis, tiap
    // direktori waktu melahirkan kembali duplikat yang baru saja dihapus, dan
    // restart dari t>0 akan membawa eps_p lama meski dict sudah diubah.
    if (phase_ != "solid")
    {
        os.writeEntry("D", D_);
    }

    os.writeEntry("DField", DFieldName_);
    os.writeEntry("rho", rhoName_);
}


makePatchTypeField
(
    fvPatchScalarField,
    speciesCoupledMixedFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
