/*---------------------------------------------------------------------------*\
  dualSite.C
\*---------------------------------------------------------------------------*/

#include "dualSite.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace kineticClosures
{
    defineTypeNameAndDebug(dualSite, 0);
    addToRunTimeSelectionTable(kineticClosure, dualSite, dictionary);
}
}


Foam::kineticClosures::dualSite::dualSite
(
    const fvMesh& mesh,
    const word& specieName,
    const dictionary& dict
)
:
    kineticClosure(mesh, specieName, dict),
    eta_(dict.get<scalar>("eta")),
    k1_(dict.get<scalar>("k1")),
    k2_(dict.get<scalar>("k2")),
    alphaSL_(dict.getOrDefault<scalar>("alphaSL", 0)),
    qsSL_(dict.getOrDefault<scalar>("qsSL", 1))
{
    const dimensionSet dimLoad(dimMoles/dimMass);

    q1_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "q1." + specieName,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("q1", dimLoad, 0)
        )
    );
    q2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "q2." + specieName,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("q2", dimLoad, 0)
        )
    );

    q1_->storeOldTime();
    q2_->storeOldTime();

    Info<< "    dualSite[" << specieName << "]: eta=" << eta_
        << " k1=" << k1_ << " k2=" << k2_ << " 1/s" << endl;
    if (alphaSL_ > 0)
    {
        Info<< "    dualSite[" << specieName << "]: self-limiting k2 exp(-"
            << alphaSL_ << " q/" << qsSL_ << ")" << endl;
    }
}


void Foam::kineticClosures::dualSite::storeOldTime()
{
    q1_->storeOldTime();
    q2_->storeOldTime();
}


void Foam::kineticClosures::dualSite::update
(
    volScalarField& q,
    const volScalarField& qStar,
    const scalar dt
)
{
    const scalar k1dt = k1_*dt;
    const scalar k2dt = k2_*dt;
    const scalar d1 = 1.0 + k1dt;
    const scalar d2 = 1.0 + k2dt;

    volScalarField& q1 = q1_();
    volScalarField& q2 = q2_();

    // Internal
    {
        scalarField& q1c = q1.primitiveFieldRef();
        scalarField& q2c = q2.primitiveFieldRef();
        scalarField& qc  = q.primitiveFieldRef();
        const scalarField& q1o = q1.oldTime().primitiveField();
        const scalarField& q2o = q2.oldTime().primitiveField();
        const scalarField& qs  = qStar.primitiveField();
        forAll(qc, c)
        {
            q1c[c] = max((q1o[c] + k1dt*eta_*qs[c])/d1, scalar(0));
            if (alphaSL_ > 0)
            {
                const scalar k2l = k2dt*Foam::exp(-alphaSL_*(q1o[c] + q2o[c])/qsSL_);
                q2c[c] = max((q2o[c] + k2l*(1 - eta_)*qs[c])/(1 + k2l), scalar(0));
            }
            else
            {
                q2c[c] = max((q2o[c] + k2dt*(1 - eta_)*qs[c])/d2, scalar(0));
            }
            qc[c]  = q1c[c] + q2c[c];
        }
    }
    // Boundaries
    forAll(q.boundaryField(), patchi)
    {
        scalarField& q1p = q1.boundaryFieldRef()[patchi];
        scalarField& q2p = q2.boundaryFieldRef()[patchi];
        scalarField& qp  = q.boundaryFieldRef()[patchi];
        const scalarField& q1o = q1.oldTime().boundaryField()[patchi];
        const scalarField& q2o = q2.oldTime().boundaryField()[patchi];
        const scalarField& qs  = qStar.boundaryField()[patchi];
        forAll(qp, f)
        {
            q1p[f] = max((q1o[f] + k1dt*eta_*qs[f])/d1, scalar(0));
            if (alphaSL_ > 0)
            {
                const scalar k2l = k2dt*Foam::exp(-alphaSL_*(q1o[f] + q2o[f])/qsSL_);
                q2p[f] = max((q2o[f] + k2l*(1 - eta_)*qs[f])/(1 + k2l), scalar(0));
            }
            else
            {
                q2p[f] = max((q2o[f] + k2dt*(1 - eta_)*qs[f])/d2, scalar(0));
            }
            qp[f]  = q1p[f] + q2p[f];
        }
    }
}


// ************************************************************************* //
