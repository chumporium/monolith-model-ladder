/*---------------------------------------------------------------------------*\
  testImplicitSpecies.C

  Uji kopling species fluida-padatan IMPLISIT vs PARTISI, tanpa aliran dan
  tanpa adsorpsi (difusi murni, dua region, BC speciesCoupledMixed):

      fluida  : ddt(rho Y)   = laplacian(rho D, Y)       rho seragam != 1
      padatan : eps ddt(C)   = laplacian(eps Ds, C)

  Jawaban yang harus dipenuhi:
    1. massa total  INT rho Y dV + INT eps C dV  kekal (mesin)
    2. keadaan akhir C = rho Y seragam = massa / (rho V_f + eps V_s)
    3. implisit dan partisi-konvergen memberi medan sama pada dt kecil
    4. implisit tetap stabil & kekal pada dt besar (Lambda >> 1)

  Kasus membaca constant/testImplicitProperties:
      rho D eps Ds Y0 dt nSteps mode(implicit|partitioned) nOuter tolOuter
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    fvMesh fluid
    (
        IOobject("fluid", runTime.timeName(), runTime, IOobject::MUST_READ)
    );
    fvMesh solid
    (
        IOobject("solid", runTime.timeName(), runTime, IOobject::MUST_READ)
    );

    IOdictionary props
    (
        IOobject
        (
            "testImplicitProperties", runTime.constant(), runTime,
            IOobject::MUST_READ
        )
    );
    const scalar rhoV = props.get<scalar>("rho");
    const scalar DV   = props.get<scalar>("D");
    const scalar epsV = props.get<scalar>("eps");
    const scalar DsV  = props.get<scalar>("Ds");
    const scalar Y0   = props.get<scalar>("Y0");
    const scalar dtV  = props.get<scalar>("dt");
    const label nSteps = props.get<label>("nSteps");
    const word mode(props.get<word>("mode"));
    const label nOuter = props.getOrDefault<label>("nOuter", 1);
    const scalar tolOuter = props.getOrDefault<scalar>("tolOuter", 1e-10);

    runTime.setDeltaT(dtV);

    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), fluid,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        fluid, dimensionedScalar(dimDensity, rhoV)
    );
    volScalarField Y
    (
        IOobject("C3H6", runTime.timeName(), fluid,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        fluid
    );

    volScalarField C
    (
        IOobject("C.C3H6", runTime.timeName(), solid,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        solid
    );

    // Ds.C3H6 terdaftar: speciesCoupledMixed sisi padatan memakai eps_p*Ds
    volScalarField DsF
    (
        IOobject("Ds.C3H6", runTime.timeName(), solid,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        solid, dimensionedScalar(dimArea/dimTime, DsV)
    );

    // useImplicit hanya pada mode implisit (solve tunggal dgn patch implisit
    // tanpa perakitan tidak sah)
    for (auto* fldPtr : {&Y, &C})
    {
        forAll(fldPtr->boundaryField(), patchi)
        {
            auto& pf = fldPtr->boundaryFieldRef()[patchi];
            if (pf.type() == "speciesCoupledMixed")
            {
                pf.useImplicit(mode == "implicit");
            }
        }
    }

    Y.primitiveFieldRef() = Y0;
    Y.correctBoundaryConditions();
    C.correctBoundaryConditions();

    const dimensionedScalar rhoD("rhoD", dimDensity*dimArea/dimTime, rhoV*DV);
    const dimensionedScalar epsD("epsD", dimArea/dimTime, epsV*DsV);
    const dimensionedScalar epsDim("eps", dimless, epsV);

    auto massa = [&]()
    {
        return gSum(rho.primitiveField()*Y.primitiveField()*fluid.V())
             + epsV*gSum(C.primitiveField()*solid.V());
    };

    const scalar m0 = massa();
    const scalar Vf = gSum(fluid.V());
    const scalar Vs = gSum(solid.V());
    const scalar Ceq = m0/(Vf + epsV*Vs);   // setimbang: C = rho Y
    Info<< "mode " << mode << " dt " << dtV << " Lambda_s ~ "
        << DsV*dtV*sqr(gMax(solid.surfaceInterpolation::deltaCoeffs().primitiveField()))
        << " m0 " << m0 << " C_eq " << Ceq << nl << endl;

    for (label step = 0; step < nSteps; ++step)
    {
        ++runTime;
        Y.storePrevIter();

        label nIt = 0;
        scalar sisa = GREAT;

        if (mode == "implicit")
        {
            fvScalarMatrix YEqn(fvm::ddt(rho, Y) - fvm::laplacian(rhoD, Y));
            fvScalarMatrix CEqn(epsDim*fvm::ddt(C) - fvm::laplacian(epsD, C));

            fvMatrix<scalar> asmb(Y, dimMass/dimTime);
            asmb.addFvMatrix(YEqn);
            asmb.addFvMatrix(CEqn);
            asmb.solve();
            nIt = 1;
        }
        else
        {
            for (nIt = 1; nIt <= nOuter; ++nIt)
            {
                const scalarField Cold(C.primitiveField());
                fvScalarMatrix YEqn(fvm::ddt(rho, Y) - fvm::laplacian(rhoD, Y));
                YEqn.solve();
                fvScalarMatrix CEqn(epsDim*fvm::ddt(C) - fvm::laplacian(epsD, C));
                CEqn.solve();
                sisa = gMax(mag(C.primitiveField() - Cold)())/max(Ceq, VSMALL);
                if (sisa < tolOuter) break;
            }
        }

        const scalar m = massa();
        Info<< "t " << runTime.value() << " iter " << nIt
            << " galatMassa " << (m - m0)/m0
            << " Ymin/max " << gMin(Y.primitiveField()) << ' ' << gMax(Y.primitiveField())
            << " rhoY min/max " << rhoV*gMin(Y.primitiveField()) << " " << rhoV*gMax(Y.primitiveField())
            << " C min/max " << gMin(C.primitiveField()) << ' ' << gMax(C.primitiveField())
            << " C_eq " << Ceq << endl;
    }

    runTime.writeNow();
    Info<< "End" << endl;
    return 0;
}
