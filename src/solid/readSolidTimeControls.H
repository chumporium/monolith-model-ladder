/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    readSolidTimeControls

Description
    Read the control parameters used in the solid

\*---------------------------------------------------------------------------*/

scalar maxDi = runTime.controlDict().getOrDefault<scalar>("maxDi", 10);

// Bilangan kopling adsorpsi -- lihat solid/adsorptionCouplingNo.H.
// BAWAAN GREAT = tidak mengikat, sehingga seluruh case lama bit-identik.
// Nilai yang disarankan 0,2 (terukur: lulus di 0,272, gagal di 0,327).
scalar maxLambda =
    runTime.controlDict().getOrDefault<scalar>("maxLambda", GREAT);

// ************************************************************************* //
