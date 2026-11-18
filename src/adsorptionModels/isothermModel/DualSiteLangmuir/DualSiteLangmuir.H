/*---------------------------------------------------------------------------*\
  DualSiteLangmuir.H

  Isoterm Dual-Site Langmuir (DSL) bergantung-suhu (komponen murni).

      q*(p,T) = qm1(T) b1 p/(1 + b1 p) + qm2(T) b2 p/(1 + b2 p)
      qmk(T)  = qm0k exp[ chi_k (1 - T/T0) ]
      bk(T)   = b0k  exp[ (dHk/(R T0)) (T0/T - 1) ]    [b0k dalam 1/Pa]

  Entri dictionary:
      T0 ; qm01 ; b01 ; deltaH1 ; chi1 (opsional) ;
           qm02 ; b02 ; deltaH2 ; chi2 (opsional) ;
\*---------------------------------------------------------------------------*/

#ifndef isotherm_DualSiteLangmuir_H
#define isotherm_DualSiteLangmuir_H

#include "isothermModel.H"

namespace Foam
{
namespace isotherms
{

class DualSiteLangmuir
:
    public isothermModel
{
    scalar T0_;
    scalar qm01_, b01_, dH1_, chi1_;
    scalar qm02_, b02_, dH2_, chi2_;

    inline scalar site
    (
        const scalar p, const scalar T,
        const scalar qm0, const scalar b0, const scalar dH, const scalar chi
    ) const
    {
        const scalar qm = qm0*Foam::exp(chi*(1 - T/T0_));
        const scalar b  = b0*Foam::exp(dH/(adsorption::Rgas*T0_)*(T0_/T - 1));
        const scalar bp = b*max(p, scalar(0));
        return qm*bp/(1 + bp);
    }
    inline scalar dSitedp
    (
        const scalar p, const scalar T,
        const scalar qm0, const scalar b0, const scalar dH, const scalar chi
    ) const
    {
        const scalar qm = qm0*Foam::exp(chi*(1 - T/T0_));
        const scalar b  = b0*Foam::exp(dH/(adsorption::Rgas*T0_)*(T0_/T - 1));
        const scalar bp = b*max(p, scalar(0));
        return qm*b/sqr(1 + bp);
    }

public:

    TypeName("DualSiteLangmuir");

    DualSiteLangmuir(const word& specieName, const dictionary& dict);

    virtual scalar qStar(const scalar p, const scalar T) const;
    virtual scalar dqStardp(const scalar p, const scalar T) const;
    virtual scalar dqStardT(const scalar p, const scalar T) const;

    //- SENGAJA mengembalikan 0 = "tidak menyatakan satu nilai".
    //
    //  DSL punya DUA energi situs (dH1_, dH2_) dengan afinitas dan kapasitas
    //  berbeda. Panas adsorpsi efektifnya bergantung pembagian muatan antar
    //  situs, jadi ia BERGANTUNG KEADAAN (p,T) -- bukan tetapan yang boleh
    //  dipulangkan lewat antarmuka ini.
    //
    //  Meringkasnya jadi satu angka (mis. rata-rata berbobot kapasitas) akan
    //  memberi nilai yang tampak masuk akal tetapi salah di sebagian besar
    //  rentang operasi, DAN menutup satu-satunya penanda bahwa penggunanya
    //  perlu berpikir. Karena itu 0 dikembalikan, sehingga `deltaH` tetap
    //  WAJIB ditulis di dict dan pilihannya terekam di sana.
    //
    //  Catatan: source energi solver memang memakai satu deltaH per species
    //  (bukan sum_k dH_k rho_s dq_k/dt). Itu aproksimasi yang sudah ada
    //  sebelumnya; keputusan ini tidak memperburuknya, hanya menolak
    //  menyembunyikannya.
    virtual scalar heatOfAdsorption() const { return scalar(0); }
};

} // End namespace isotherms
} // End namespace Foam

#endif

// ************************************************************************* //
