// createAxialBins.H -- sekali saat start (Tahap 2.1, lihat axialBins.H).
// system/axialBins { nBins <label>; interval <s>; }  -- tanpa berkas: nonaktif, bit-identik.
label axialBinsN = 0;
scalar axialBinsInterval = GREAT;
scalar axialBinsZmin = 0, axialBinsZmax = 1;
{
    IOdictionary abDict
    (
        IOobject
        (
            "axialBins", runTime.system(), runTime,
            IOobject::READ_IF_PRESENT, IOobject::NO_WRITE, false
        )
    );
    if (abDict.headerOk())
    {
        if (fluidRegions.size() != 1 || solidRegions.size() != 1)
        {
            FatalIOErrorInFunction(abDict)
                << "axialBins saat ini untuk satu region fluida + satu padatan."
                << exit(FatalIOError);
        }
        axialBinsN = abDict.get<label>("nBins");
        axialBinsInterval = abDict.get<scalar>("interval");
        if (axialBinsN < 1 || axialBinsInterval <= 0)
        {
            FatalIOErrorInFunction(abDict)
                << "nBins >= 1 dan interval > 0" << exit(FatalIOError);
        }
        // Batas global (polyMesh::bounds tereduksi antar prosesor)
        axialBinsZmin = fluidRegions[0].bounds().min().z();
        axialBinsZmax = fluidRegions[0].bounds().max().z();
        Info<< "axialBins: " << axialBinsN << " potongan, z " << axialBinsZmin
            << " .. " << axialBinsZmax << " m, interval " << axialBinsInterval
            << " s" << nl << endl;
    }
}
