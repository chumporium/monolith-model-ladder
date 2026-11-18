/*---------------------------------------------------------------------------*\
  localEquilibrium.H

  Closure local-equilibrium / instan (Bab 7.1): q = q*(C,T). q bukan variabel
  PDE; storage transien ditangani regionAdsorption via Jacobian (effStorage).
  Valid untuk fisisorpsi makropori (laju dikendalikan transport teresolusi).
\*---------------------------------------------------------------------------*/

#ifndef closure_localEquilibrium_H
#define closure_localEquilibrium_H

#include "kineticClosure.H"

namespace Foam
{
namespace kineticClosures
{

class localEquilibrium
:
    public kineticClosure
{
public:

    TypeName("localEquilibrium");

    localEquilibrium
    (
        const fvMesh& mesh,
        const word& specieName,
        const dictionary& dict
    );

    virtual bool isAlgebraic() const { return true; }

    virtual void update
    (
        volScalarField& q,
        const volScalarField& qStar,
        const scalar dt
    );
};

} // End namespace kineticClosures
} // End namespace Foam

#endif

// ************************************************************************* //
