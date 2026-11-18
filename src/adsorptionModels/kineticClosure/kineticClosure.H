/*---------------------------------------------------------------------------*\
  kineticClosure.H

  Base class pluggable untuk CLOSURE KINETIK SOLID (Bab 7 panduan). Menentukan
  bagaimana loading q_i menuju q_i* (kesetimbangan).

  Tiga mode:
    - localEquilibrium : q = q*  (aljabar; storage via Jacobian di species eqn)
    - LDF              : dq/dt = k (q* - q)         (finite-rate)
    - dualSite         : dua laju, q = q1 + q2      (ekor asimetris)

  Antarmuka:
    isAlgebraic()   true untuk localEquilibrium (regionAdsorption memakai
                    bentuk effStorage); false untuk closure berbasis-laju
                    (species eqn memakai source eksplisit rho_s dq/dt).
    update(q,qStar,dt)  perbarui field loading q [mol/kg].
    storeOldTime()      simpan oldTime field state internal (mis. q1,q2).
\*---------------------------------------------------------------------------*/

#ifndef kineticClosure_H
#define kineticClosure_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

namespace Foam
{

// Deklarasi maju: hanya dipakai sbg referensi pada setIsothermContext, jadi
// tak perlu menyeret isothermModel.H ke setiap penyertaan closure.
class isothermModel;

class kineticClosure
{
protected:

    const fvMesh& mesh_;
    const word specieName_;
    const dictionary coeffs_;

public:

    TypeName("kineticClosure");

    declareRunTimeSelectionTable
    (
        autoPtr,
        kineticClosure,
        dictionary,
        (
            const fvMesh& mesh,
            const word& specieName,
            const dictionary& dict
        ),
        (mesh, specieName, dict)
    );

    kineticClosure
    (
        const fvMesh& mesh,
        const word& specieName,
        const dictionary& dict
    );

    kineticClosure(const kineticClosure&) = delete;

    static autoPtr<kineticClosure> New
    (
        const fvMesh& mesh,
        const word& specieName,
        const dictionary& dict
    );

    virtual ~kineticClosure() = default;


    // Member Functions

        //- true jika q = q* aljabar (local-equilibrium)
        virtual bool isAlgebraic() const = 0;

        //- Perbarui loading q [mol/kg] menuju qStar pada langkah dt
        virtual void update
        (
            volScalarField& q,
            const volScalarField& qStar,
            const scalar dt
        ) = 0;

        //- Simpan oldTime untuk field state internal (default: tidak ada)
        virtual void storeOldTime() {}

        //- Sediakan konteks isoterm+suhu+fActive kpd closure yang butuh.
        //
        //  Dipanggil SEKALI oleh regionAdsorption setelah closure dibuat.
        //  Hanya intraGrainDiffusion (koreksi Darken) yang memakainya; closure
        //  lain mengabaikannya.
        //
        //  Kenapa lewat setter, bukan lewat tabel RTS: menambah argumen ke
        //  tabel akan memaksa keempat closure + New() berubah walau tiga di
        //  antaranya tak butuh apa pun. Setter menjaga perubahan terbatas pada
        //  yang benar-benar memerlukan.
        //
        //  ⚠️ fActive WAJIB diteruskan: regionAdsorption menyimpan q dan qStar
        //  sbg fActive*q_isoterm (regionAdsorption.C:54), sedangkan isoterm
        //  merumuskan Gamma pada skalanya sendiri.
        virtual void setIsothermContext
        (
            const isothermModel& iso,
            const volScalarField& T,
            const scalar fActive
        )
        {}

        //- Faktor respons implisit g = d(q_)/d(qStar) pada langkah dt.
        //  Dipakai regionAdsorption (jalur finite-rate) utk melinearisasi sink
        //  secara implisit: dS/dC = g*cross_ii (cross_ii = rho_s M dq*/dC).
        //  Suku fvm::Sp(g*cross_ii/dt, C) + deferred-correction -> stabil (cegah
        //  C<0 pd closure KAKU) & KONSERVATIF (batal saat konvergen -> sink sekan
        //  eksak). Default 0 = sink eksplisit murni (perilaku lama, mis. dualSite).
        //  LDF & intraGrainDiffusion override. Lihat CATATAN BAGIAN 13/14.
        virtual scalar dqImplicitFactor(const scalar dt) const { return 0; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
