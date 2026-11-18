# monolith-model-ladder

`adsorptionMultiRegionFoam` is a transient OpenFOAM solver for adsorption in structured contactors such as
monoliths. It solves laminar gas flow in the channels and heat and mass transfer in the porous sorbent wall as two
coupled regions. It extends `chtMultiRegionFoam` (OpenFOAM v2412) with pore diffusion in the wall, temperature-dependent
isotherms, uptake rate models, the heat of adsorption, and species coupling at the fluid-wall interface.

The solver was used for the three-dimensional channel and module simulations in

> E. Munandar, Nasruddin, A.P. Sasmito, *A model-ladder framework for monolith adsorption cycle models: derived
> parameters, reduction error and computing cost*, submitted to Chemical Engineering Science.

Release `v1.0-CES` is the version used in that article.

## Contents

| Path | Content |
|---|---|
| `src/` | solver, adsorption model library and coupled boundary conditions |
| `src/utilities/` | auxiliary solvers: `effCondFoam` (effective conductivity of a cell), `periodicDecayFoam`, `wallShapeFoam` |
| `tutorials/singleChannelBreakthrough/` | CO2 breakthrough in one square channel of an amine monolith |

## Requirements and build

OpenFOAM v2412 (openfoam.com) with a C++ compiler and MPI. `Allwmake` checks the OpenFOAM version and by default
installs the binaries in a separate user directory, `$HOME/OpenFOAM/bundle-v2412`, so that they do not overwrite other
user applications. Add that directory to the environment before running a case:

```
source /usr/lib/openfoam/openfoam2412/etc/bashrc    # or the path of your installation
cd src
./Allwmake -j
export PATH=$HOME/OpenFOAM/bundle-v2412/platforms/$WM_OPTIONS/bin:$PATH
export LD_LIBRARY_PATH=$HOME/OpenFOAM/bundle-v2412/platforms/$WM_OPTIONS/lib:$LD_LIBRARY_PATH
```

`BUNDLE_ISOLATED=0 ./Allwmake -j` installs into the standard `$FOAM_USER_APPBIN` and `$FOAM_USER_LIBBIN` instead,
which are already on the search path. The build produces the application `adsorptionMultiRegionFoam` and the libraries
`libadsorptionModels.so` and `libadsorptionCoupledBCs.so`. Cases load the boundary-condition library through `libs` in
`system/controlDict`. The utilities in `src/utilities` are built separately with `wmake` in their directories.

## Tutorial

`tutorials/singleChannelBreakthrough` is a unit cell of a square-channel monolith: a 2.0 mm channel with a 0.2 mm half-wall
(2.4 mm pitch) and a length of 80 mm. A dry feed of 400 ppm CO2 in N2 enters at 0.72 m/s and 298 K. The wall is a
porous mullite support with an amine-functionalized alumina wash-coat (Toth isotherm, heat of adsorption 70 kJ/mol).
The outer faces of the wall are symmetry planes, so the case represents a channel inside a monolith. The mesh has
about 1600 cells and is kept coarse for a short run: 300 s of simulated time take about 11 minutes on one core. Finer
meshes (more cells across the channel and the wall and along the axis) are needed for converged results; the
parameters `nz`, `nf` and `nw` in `system/blockMeshDict` set the resolution.

```
cd tutorials/singleChannelBreakthrough
./Allrun                     # blockMesh, splitMeshRegions, adsorptionMultiRegionFoam
python3 plotBreakthrough.py  # outlet c/c0, mean loading and mean wall temperature against time
./Allclean
```

The outlet CO2 mass fraction is written to `postProcessing/fluid/outletCO2` and the mean loading and temperature of the
wall to `postProcessing/solid/solidAverage`. The fields `q.CO2` (loading) and `qStar.CO2` (equilibrium loading), both
in mol per kg of wall, are written in the solid region and can be viewed with `paraFoam -region solid`.

With the mesh as supplied, about 2 % of the feed CO2 leaves the channel from the start, the outlet reaches half the
feed concentration after about 165 s, and the wall is saturated at 0.036 mol/kg after 300 s. The mean wall
temperature rises by about 0.7 K. These values serve to check an installation.

## Case set-up

A case has two regions, `fluid` and `solid`, listed in `constant/regionProperties`. The fluid region uses the standard
`heRhoThermo` with `multiComponentMixture`; the gas species and their molecular diffusivities are given in
`constant/fluid/thermophysicalProperties` (`species`, `inertSpecie`, `speciesDiffusivity`). The solid region uses
`heSolidThermo` and in addition reads `constant/solid/adsorptionProperties`.

### adsorptionProperties

| Keyword | Unit | Meaning |
|---|---|---|
| `species` | | list of adsorbing gas species; each needs a sub-dictionary |
| `mixingRule` | | `ideal`, `extendedLangmuir`, `IAST` or `coAdsorption` (multicomponent equilibrium) |
| `epsilonP` | - | transport (macro) porosity of the wall |
| `epsilonPtotal` | - | total porosity used for gas storage (default `epsilonP`) |
| `rhoS` | kg/m3 | density of the wall |
| `fActive` | - | mass fraction of active sorbent in the wall (default 1) |
| `isothermalAdsorption` | | `yes` sets the heat of adsorption to zero (default `no`) |
| `nNewton`, `newtonTol` | | iterations and tolerance of the local equilibrium solve (default 50, 1e-8) |

Per species:

| Keyword | Unit | Meaning |
|---|---|---|
| `molWeight` | g/mol | molar mass |
| `Ds` | m2/s | effective diffusivity of the species in the pores of the wall |
| `deltaH` | J/mol | heat of adsorption (default: `deltaH0` of the isotherm) |
| `isotherm { type ...; }` | | equilibrium loading q*, in mol per kg of active sorbent |
| `kinetics { type ...; }` | | uptake rate model |
| `bosanquet { active; Dbulk; Dknudsen or rPore; pRef; TRef; }` | | optional Bosanquet combination of bulk and Knudsen diffusion that replaces `Ds` |

Isotherms (p in Pa, T in K; `b0` and `Kref` in 1/Pa, heats in J/mol):

| `type` | Form | Keywords |
|---|---|---|
| `Langmuir` | q* = qm b p / (1 + b p) | `T0 qm0 b0 deltaH0 chi` |
| `DualSiteLangmuir` | sum of two Langmuir sites | `T0 qm01 b01 deltaH1 chi1 qm02 b02 deltaH2 chi2` |
| `Sips` | q* = qm (b p)^n / (1 + (b p)^n), n = n0 + alpha_n (1 - T0/T) | `T0 qm0 b0 n0 deltaH0 chi alpha_n` |
| `Toth` | q* = ns b p / (1 + (b p)^t)^(1/t), t = t0 + alpha_t (1 - T0/T) | `T0 ns0 b0 t0 deltaH0 chi alpha_t` |
| `NittaMultisite` | theta = K p (1 - theta)^a, q* = qmax theta | `qmax a Kref Tref deltaH0` |
| `asymptoticMatched` | q* = (q0^n + qp^n)^(1/n), q0 = K p/(R T), qp from sub-dictionary `base` | `base K n` |

In all forms b = b0 exp[deltaH0/(R T0) (T0/T - 1)] and the capacity is qm0 exp[chi (1 - T/T0)].

Uptake rate models:

| `type` | Model | Keywords |
|---|---|---|
| `localEquilibrium` | q = q*(p, T) in each cell; the rate is set by pore diffusion in the wall | none |
| `LDF` | dq/dt = k (q* - q) | `kLDF`, or `Dc` and `rc` with k = 15 Dc/rc^2 |
| `dualSite` | two sites with rates k1 and k2 and fraction `eta` of fast sites | `k1 k2 eta`, optional `alphaSL qsSL` |
| `intraGrainDiffusion` | diffusion in a sorbent film or grain, resolved on `nShells` shells | `Dc delta nShells geometry` and optional terms |

### Boundary conditions at the fluid-wall interface

| Field | Region | Type |
|---|---|---|
| species mass fraction (e.g. `CO2`) | fluid | `speciesCoupledMixed` with `phase fluid; nbrField C.<species>; D <molecular diffusivity>; rho rho;` |
| pore gas concentration `C.<species>` [kg/m3] | solid | `speciesCoupledMixed` with `phase solid; nbrField <species>; DField Ds.<species>; rho rho;` |
| `T` | both | `compressible::turbulentTemperatureRadCoupledMixed` |

`danckwertsInlet` and `danckwertsInletTemperature` are available for inlet conditions with axial dispersion.

### Solution controls

The coupling between the regions is controlled in the `PIMPLE` dictionary of `system/fvSolution`. `nOuterCorrectors`
sets the number of outer iterations per time step. Optional entries are `speciesCoupling implicit` and
`energyInterfaceCoupling implicit` (fluid-wall transfer solved in one system), `energySpeciesEnthalpy consistent`
(enthalpy carried by the transferred species), and `interfaceTolerance` and `interfaceRelaxation` for the outer loop.
The tutorial uses the explicit (default) coupling. The solid time step is limited by `maxDi` in
`system/controlDict` in addition to the Courant number of the fluid.

## Licence

GNU General Public License v3.0 (see `LICENSE`), as for OpenFOAM, from which the solver is derived.
