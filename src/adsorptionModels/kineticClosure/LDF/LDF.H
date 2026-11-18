/*---------------------------------------------------------------------------*\
  LDF.H

  Closure Linear Driving Force tunggal (Bab 7.2):
      dq/dt = k (q* - q)

  Integrasi Euler implisit per sel (tak-bersyarat stabil):
      q^{n+1} = (q^n + k dt q*) / (1 + k dt)

  Tetapan laju k [1/s] diberikan langsung, atau diturunkan dari difusi
  intrakristalin: k = 15 Dc / rc^2 (entri kLDF, atau Dc & rc).

  Valid untuk: fisisorpsi mikropori / kinetika kimia lapisan aktif (Bab 7.4).
\*---------------------------------------------------------------------------*/

#ifndef closure_LDF_H
#define closure_LDF_H

#include "kineticClosure.H"

namespace Foam
{
namespace kineticClosures
{

class LDF
:
    public kineticClosure
{
    scalar k_;

public:

    TypeName("LDF");

    LDF
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

    //- g = dq/dq* = kdt/(1+kdt) (dari q=(q_old+kdt q*)/(1+kdt))
    virtual scalar dqImplicitFactor(const scalar dt) const
    {
        const scalar kdt = k_*dt;
        return kdt/(1.0 + kdt);
    }
};

} // End namespace kineticClosures
} // End namespace Foam

#endif

// ************************************************************************* //
