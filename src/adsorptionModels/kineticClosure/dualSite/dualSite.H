/*---------------------------------------------------------------------------*\
  dualSite.H

  Closure dual-site dua-laju (Bab 7.3):
      q = q1 + q2
      dq1/dt = k1 ( eta q*     - q1 )
      dq2/dt = k2 ( (1-eta) q* - q2 )

  eta   = fraksi situs cepat (permukaan)
  Opsi difusi amina "self-limiting" (Grossmann & Mazzotti, IECR 64 (2025) 2339, pers. 40):
      k2 = k2 exp(-alphaSL q/qsSL),  q = q1 + q2 waktu lama (tertinggal satu langkah)
  alphaSL (bawaan 0 = laju tetap, jalur lama persis), qsSL (satuan sama dengan q).
  1-eta = fraksi situs lambat (bulk amine)  -> ekor breakthrough asimetris.

  q1, q2 adalah field state internal (di-write untuk restart). Integrasi Euler
  implisit per sel. Valid untuk: kemisorpsi amine (Bab 7.4).
\*---------------------------------------------------------------------------*/

#ifndef closure_dualSite_H
#define closure_dualSite_H

#include "kineticClosure.H"

namespace Foam
{
namespace kineticClosures
{

class dualSite
:
    public kineticClosure
{
    scalar eta_, k1_, k2_;
    scalar alphaSL_, qsSL_;

    autoPtr<volScalarField> q1_;
    autoPtr<volScalarField> q2_;

public:

    TypeName("dualSite");

    dualSite
    (
        const fvMesh& mesh,
        const word& specieName,
        const dictionary& dict
    );

    virtual bool isAlgebraic() const { return false; }

    virtual void update
    (
        volScalarField& q,
        const volScalarField& qStar,
        const scalar dt
    );

    // Faktor sink implisit (fix-b): dq/dq* dari Euler implisit per situs,
    // dijumlah berbobot eta. Tanpa ini closure finite-rate BOCOR massa di
    // regime kaku (saga P2: HSDM/LDF k-tinggi bocor 89-96%).
    virtual scalar dqImplicitFactor(const scalar dt) const
    {
        const scalar k1dt = k1_*dt;
        const scalar k2dt = k2_*dt;
        return eta_*k1dt/(1.0 + k1dt) + (1 - eta_)*k2dt/(1.0 + k2dt);
    }

    virtual void storeOldTime();
};

} // End namespace kineticClosures
} // End namespace Foam

#endif

// ************************************************************************* //
