// updatePoreInterfaceFlux.H (Paper B, tahap V2): terbitkan fluks antarmuka pori untuk sisi fluida,
// SESUDAH species padat diselesaikan (C muka = nilai yang dipakai div(phiDarcy,C)).
//   NtPore          [mol/s] fluks mol total;  FadvPore.<sp> [kg/s] adveksi species terlacak.
// Keduanya lewat poreFluxCarry (konservatif terhadap iterasi luar).
if (phiDarcyFields.set(i))
{
    regionAdsorption& adsI = adsorption[i];
    NtCarry[i].update("NtPore", NtPoreFields[i]);
    if (FadvCarry[i].size() != adsI.species().size())
    {
        FadvCarry[i].setSize(adsI.species().size());
    }
    forAll(adsI.species(), si)
    {
        if (!FadvCarry[i].set(si)) { FadvCarry[i].set(si, new poreFluxCarry()); }
        // keluar pori: C sel x phiCOut;  masuk pori: -|phiDarcy| (rho Y) sel kanal (poreInFlux)
        surfaceScalarField Fraw
        (
            "FadvRaw", fvc::interpolate(adsI.C(si))*phiCFields[i]
        );
        forAll(Fraw.boundaryField(), patchi)
        {
            Fraw.boundaryFieldRef()[patchi] +=
                adsI.C(si).boundaryField()[patchi].patchInternalField()
               *phiCOutFields[i].boundaryField()[patchi];
            if (poreInFluxFields[i].size() > si && poreInFluxFields[i].set(si))
            {
                Fraw.boundaryFieldRef()[patchi] += poreInFluxFields[i][si].boundaryField()[patchi];
            }
        }
        FadvCarry[i][si].update("FadvPore." + adsI.species()[si], Fraw, fadvLag);
    }
}

// porePressureFromSpecies (REGISTER 2026-10-18, PB2-N2; solid PIMPLE, bawaan MATI): bila SEMUA komponen gas pori diangkut sebagai species
// (N2 eksplisit), inventori pori ditentukan species, bukan persamaan tekanan. Tanpa penyelarasan, gas kanal yang masuk pori (phiDarcy < 0)
// menumpuk sebagai C.N2 tanpa menaikkan p pori, sehingga tak ada yang mendorongnya keluar (N3/N2_L6: C.N2 0,011 -> 0,59 kg/m3 dalam 1 s).
// Di sini, sesudah species diselesaikan, p pori sel := R T SUM_j C_j/M_j. Persamaan tekanan tetap memprediksi fluks Darcy/Knudsen tiap
// iterasi; simpanannya (oldTime) kini selalu = inventori species, jadi kelebihan species menaikkan p dan terdorong keluar.
if
(
    pPoreFields.set(i)
 && mesh.solutionDict().subDict("PIMPLE").getOrDefault<bool>("porePressureFromSpecies", false)
)
{
    regionAdsorption& adsP = adsorption[i];
    volScalarField& pPs = pPoreFields[i];
    const scalarField& Tc = thermo.T().primitiveField();
    scalarField cSum(mesh.nCells(), Zero);
    forAll(adsP.species(), sj)
    {
        cSum += adsP.C(sj).primitiveField()/adsP.M(sj);
    }
    pPs.primitiveFieldRef() = 8.314462618*Tc*cSum;
    volScalarField& pTh = thermo.p();
    pTh.primitiveFieldRef() = pPs.primitiveField();
    pTh.correctBoundaryConditions();
}
