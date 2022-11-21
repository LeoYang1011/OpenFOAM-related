/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "spanwiseAverage.H"
#include "volFields.H"
#include "meshStructure.H"
#include "globalIndex.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::spanwiseAverage::meshAddressing()
{
    if (!meshStructurePtr_)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Count
        label sz = 0;
        sz += pbm[patchID_].size();

        // Fill
        labelList meshFaces(sz);
        sz = 0;
        label start = pbm[patchID_].start();
        label size = pbm[patchID_].size();
        for (label i = 0; i < size; ++i)
        {
            meshFaces[sz++] = start+i;
        }

        if (sz == 0)
        {
            // TODO: If source patch is a cyclic it may have have been
            // converted to a processorCyclic for parallel runs

            WarningInFunction
                << "Requested patches have zero faces"
                << endl;
        }

        uindirectPrimitivePatch uip
        (
            UIndirectList<face>(mesh_.faces(), meshFaces),
            mesh_.points()
        );

        globalFaces_.set(new globalIndex(uip.size()));
        globalEdges_.set(new globalIndex(uip.nEdges()));
        globalPoints_.set(new globalIndex(uip.nPoints()));
        meshStructurePtr_.reset
        (
            new meshStructure
            (
                mesh_,
                uip,
                globalFaces_(),
                globalEdges_(),
                globalPoints_()
            )
        );
    }
}


const Foam::word Foam::spanwiseAverage::averageName
(
    const word& fieldName
)
{
    return "spanwiseAverage(" + fieldName + ")";
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spanwiseAverage::spanwiseAverage
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    patchID_(mesh.boundaryMesh().findPatchID(dict.get<word>("patchName"))),
    fieldsName_(dict.get<wordList>("fields"))
{
    meshAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::spanwiseAverage::fields()
{
    return fieldsName_;
}

void Foam::spanwiseAverage::execute()
{

    for (const word& fieldName : fieldsName_)
    {
        spanwiseAverageField<scalar>(fieldName);
        spanwiseAverageField<vector>(fieldName);
        spanwiseAverageField<sphericalTensor>(fieldName);
        spanwiseAverageField<symmTensor>(fieldName);
        spanwiseAverageField<tensor>(fieldName);
    }

}


void Foam::spanwiseAverage::write()
{
    for (const word& fieldName : fieldsName_)
    {
        const regIOobject* obj =
            mesh_.cfindObject<regIOobject>(averageName(fieldName));

        if (obj)
        {
            obj->write();
        }
    }
}


// ************************************************************************* //
