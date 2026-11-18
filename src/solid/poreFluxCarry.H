/*---------------------------------------------------------------------------*\
  poreFluxCarry.H  (Paper B, 2026-10-04, tahap V2)

  Fluks antarmuka pori -> kanal yang KONSERVATIF terhadap iterasi luar.

  Masalah (PB0-V2-1d/1e): fluida pada iterasi luar terakhir memakai fluks pori
  dari penyelesaian padatan SEBELUMNYA; padatan lalu diselesaikan sekali lagi.
  Yang dilepas pori (fluks akhir) != yang diterima kanal (fluks yang dipakai),
  dan selisihnya menumpuk (4 % pada evakuasi dengan relaksasi p_rgh 0,7).

  Obat: selisih langkah lalu dibawa ke langkah ini,
      terbit = mentah + bawaan,   bawaan = utang/dt,
      utang += (mentah_akhir - terbit_dipakai)_lalu x dt_lalu   (kumulatif)
  sehingga jumlah kumulatif yang diterima kanal = yang dilepas pori, tanpa
  bergantung pada relaksasi atau jumlah iterasi luar. Bawaan tidak ditulis ke
  disk: pada restart selisih satu langkah hilang.
\*---------------------------------------------------------------------------*/
#ifndef poreFluxCarry_H
#define poreFluxCarry_H

#include "surfaceFields.H"

namespace Foam
{

class poreFluxCarry
{
    autoPtr<surfaceScalarField> pub_;      // terdaftar; dibaca sisi fluida
    autoPtr<surfaceScalarField> carry_;
    autoPtr<surfaceScalarField> ledger_;   // utang kumulatif (fluks x waktu)
    autoPtr<surfaceScalarField> readLast_; // nilai terbit sebelum pembaruan terakhir
    autoPtr<surfaceScalarField> read2_;    // nilai terbit dua pembaruan sebelumnya
    autoPtr<surfaceScalarField> rawLast_;
    label tIndex_ = -1;

public:

    //- lag = 1: fluks dipakai fluida pada iterasi luar yang sama dengan pembacaannya (sumber massa, pEqn);
    //  lag = 2: dipakai satu iterasi luar kemudian (sumber species di YEqn).
    void update(const word& name, const surfaceScalarField& raw, const label lag = 1)
    {
        const Time& runTime = raw.mesh().time();

        if (!pub_)
        {
            pub_.reset
            (
                new surfaceScalarField
                (
                    IOobject(name, runTime.timeName(), raw.mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                    raw
                )
            );
            carry_.reset(new surfaceScalarField(name + "Carry", 0.0*raw));
            ledger_.reset(new surfaceScalarField(name + "Ledger", 0.0*raw));
            readLast_.reset(new surfaceScalarField(name + "Read", raw));
            read2_.reset(new surfaceScalarField(name + "Read2", raw));
            rawLast_.reset(new surfaceScalarField(name + "RawLast", raw));
            tIndex_ = runTime.timeIndex();
            return;
        }

        if (runTime.timeIndex() != tIndex_)
        {
            // buku utang [jumlah]: yang dilepas pori - yang diterima kanal, kumulatif
            ledger_() += (rawLast_() - (lag == 2 ? read2_() : readLast_()))*runTime.deltaT0Value();
            carry_() = ledger_()/runTime.deltaTValue();
            tIndex_ = runTime.timeIndex();
        }

        read2_() = readLast_();
        readLast_() = pub_();
        rawLast_() = raw;
        pub_() = raw + carry_();
    }
};

} // End namespace Foam

#endif
