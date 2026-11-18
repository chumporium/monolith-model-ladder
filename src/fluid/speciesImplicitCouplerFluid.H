/*---------------------------------------------------------------------------*\
  speciesImplicitCouplerFluid.H

  Implementasi kait regionAdsorption::speciesImplicitCoupler di sisi solver:
  merakit persamaan Y_i fluida (bentuk IDENTIK dengan YEqn.H) bersama
  persamaan C_i padatan hasil linearisasi Newton, lalu menyelesaikannya dalam
  satu matriks. Koefisien antarmuka ditulis speciesCoupledMixed::manipulateMatrix.

  Aktif bila  PIMPLE { speciesCoupling implicit; }  (createSpeciesCoupling.H).

  Batas implementasi (diperiksa saat start): satu region fluida, satu region
  padatan, dan kedua sisi tiap muka antarmuka di prosesor yang sama.
\*---------------------------------------------------------------------------*/

#ifndef speciesImplicitCouplerFluid_H
#define speciesImplicitCouplerFluid_H

namespace Foam
{

class speciesImplicitCouplerFluid
:
    public speciesImplicitCoupler
{
    // Referensi ke daftar region fluida solver
    PtrList<CombustionModel<rhoReactionThermo>>& reaction_;
    PtrList<volScalarField>& rho_;
    PtrList<surfaceScalarField>& phi_;
    PtrList<compressible::turbulenceModel>& turbulence_;
    PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>& fields_;
    PtrList<fv::options>& fvOptions_;
    const List<scalarList>& Dmol_;
    const scalar Sct_;

    const label fluidi_;           // region fluida (tunggal)
    labelList yIndex_;             // species padatan i -> indeks Y fluida (-1 = tak dikopel)
    PtrList<scalarField> Ykeep_;   // Y sebelum solve (untuk redaman)

    bool finalIter_;

    //- Sumber adveksi gas pori per species (tahap V2; fluidPoreSource.H). nullptr = tak dipakai.
    const List<PtrList<volScalarField>>* Sadv_ = nullptr;

public:

    void setAdvSource(const List<PtrList<volScalarField>>* s) noexcept { Sadv_ = s; }

    speciesImplicitCouplerFluid
    (
        PtrList<CombustionModel<rhoReactionThermo>>& reaction,
        PtrList<volScalarField>& rho,
        PtrList<surfaceScalarField>& phi,
        PtrList<compressible::turbulenceModel>& turbulence,
        PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>& fields,
        PtrList<fv::options>& fvOptions,
        const List<scalarList>& Dmol,
        const scalar Sct,
        const label fluidi,
        const wordList& solidSpecies
    )
    :
        reaction_(reaction),
        rho_(rho),
        phi_(phi),
        turbulence_(turbulence),
        fields_(fields),
        fvOptions_(fvOptions),
        Dmol_(Dmol),
        Sct_(Sct),
        fluidi_(fluidi),
        yIndex_(solidSpecies.size(), -1),
        Ykeep_(solidSpecies.size()),
        finalIter_(false)
    {
        const auto& comp = reaction_[fluidi_].thermo().composition();
        forAll(solidSpecies, i)
        {
            yIndex_[i] = comp.species().find(solidSpecies[i]);
        }
    }

    void setFinalIter(const bool f) { finalIter_ = f; }

    virtual bool coupled(const label i) const
    {
        return yIndex_[i] >= 0;
    }

    virtual void solve(const label i, fvScalarMatrix& CiEqn)
    {
        const label j = yIndex_[i];
        CombustionModel<rhoReactionThermo>& reaction = reaction_[fluidi_];
        const volScalarField& rho = rho_[fluidi_];
        const surfaceScalarField& phi = phi_[fluidi_];
        const fvMesh& mesh = rho.mesh();
        volScalarField& Yi = reaction.thermo().composition().Y()[j];

        Ykeep_.set(i, new scalarField(Yi.primitiveField()));

        // --- Persamaan Y_i: SAMA persis dengan YEqn.H ---
        tmp<fv::convectionScheme<scalar>> mvConvection
        (
            fv::convectionScheme<scalar>::New
            (
                mesh,
                fields_[fluidi_],
                phi,
                mesh.divScheme("div(phi,Yi_h)")
            )
        );

        const dimensionedScalar Dmol("Dmol", dimArea/dimTime, Dmol_[fluidi_][j]);
        volScalarField DEff
        (
            "DEff." + Yi.name(),
            rho*Dmol + turbulence_[fluidi_].mut()/Sct_
        );

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          - fvm::laplacian(DEff, Yi)
         ==
            reaction.R(Yi)
          + fvOptions_[fluidi_](rho, Yi)
        );

        // Tahap V2 + kopling implisit: sumber adveksi gas pori (sama dengan YEqn.H jalur partisi)
        if (Sadv_ && (*Sadv_)[fluidi_].size() > j && (*Sadv_)[fluidi_].set(j))
        {
            YiEqn -= (*Sadv_)[fluidi_][j];
        }

        YiEqn.relax();
        fvOptions_[fluidi_].constrain(YiEqn);

        // --- Rakit dan selesaikan bersama C_i ---
        fvMatrix<scalar> asmb(Yi, dimMass/dimTime);
        asmb.addFvMatrix(YiEqn);
        asmb.addFvMatrix(CiEqn);
        asmb.solve(mesh.solverDict(finalIter_ ? word("YiFinal") : word("Yi")));

        fvOptions_[fluidi_].correct(Yi);
        // TANPA clamp_min(0): clamp membuang/menambah massa. Nilai negatif
        // sekecil derau solver tetap ada dan dilaporkan neraca global.

        updateInert();
    }

    virtual void damp(const label i, const scalar alpha)
    {
        const label j = yIndex_[i];
        volScalarField& Yi =
            reaction_[fluidi_].thermo().composition().Y()[j];
        scalarField& Yp = Yi.primitiveFieldRef();
        Yp = Ykeep_[i] + alpha*(Yp - Ykeep_[i]);
        Yi.correctBoundaryConditions();
        updateInert();
    }

    //- Y inert = 1 - sum Y aktif (sama dgn YEqn.H)
    void updateInert()
    {
        rhoReactionThermo& thermo = reaction_[fluidi_].thermo();
        basicMultiComponentMixture& composition = thermo.composition();
        PtrList<volScalarField>& Y = composition.Y();
        if (!Y.size()) return;

        const label inertIndex =
            composition.species().find(thermo.get<word>("inertSpecie"));

        volScalarField Yt(0.0*Y[0]);
        forAll(Y, k)
        {
            if (k != inertIndex && composition.active(k))
            {
                Yt += Y[k];
            }
        }
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].clamp_min(0);
    }
};

} // End namespace Foam

#endif
