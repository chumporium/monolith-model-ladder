#!/bin/bash
# bangun_wsl.sh -- bangun solver bundle di WSL (isolasi bawaan Allwmake -> ~/OpenFOAM/bundle-v2412); 4 core.
. /usr/lib/openfoam/openfoam2412/etc/bashrc 2>/dev/null
cd $(dirname "$0")/../../src || exit 1
WM_NCOMPPROCS=4 ./Allwmake > $HOME/OpenFOAM/bangun_terakhir.log 2>&1; rc=$?
grep -E "error|Error" $HOME/OpenFOAM/bangun_terakhir.log | grep -v "Werror\|error\.H\|errorManip" | head -20
tail -3 $HOME/OpenFOAM/bangun_terakhir.log; ls -la --time-style=+%F_%T $HOME/OpenFOAM/bundle-v2412/platforms/*/bin/adsorptionMultiRegionFoam; exit $rc
