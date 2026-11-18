# monolith-model-ladder

Solver source code for the article *Effect of model reduction and derived parameters on the accuracy and computing
cost of monolith adsorption cycle models* (E. Munandar, Nasruddin, A.P. Sasmito; submitted to Chemical Engineering Science).

## Contents

`src/` contains `adsorptionMultiRegionFoam`, a transient solver for conjugate heat and mass transfer with adsorption
between fluid and sorbent regions, built on OpenFOAM v2412 (extends `chtMultiRegionFoam`). It is the solver used for
the three-dimensional channel and module simulations in the article.

## Building the solver

```
source <OpenFOAM v2412 installation>/etc/bashrc
cd src && ./Allwmake
```

## Licence

GNU General Public License v3.0 (see `LICENSE`), as required for code derived from OpenFOAM.
