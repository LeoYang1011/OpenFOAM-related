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
#include "meshStructure.H"
#include "globalIndex.H"
#include "polyPatch.H"
#include "pointField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshStructure& Foam::spanwiseAverage::meshAddressing(const polyMesh& mesh) const
{
    if (!meshStructurePtr_)
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Count
        label sz = 0;
        sz += pbm[sidePatchID_].size();

        // Fill
        labelList meshFaces(sz);
        sz = 0;
        label start = pbm[sidePatchID_].start();
        label size = pbm[sidePatchID_].size();
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
            UIndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        );

        globalFaces_.set(new globalIndex(uip.size()));
        globalEdges_.set(new globalIndex(uip.nEdges()));
        globalPoints_.set(new globalIndex(uip.nPoints()));
        meshStructurePtr_.reset
        (
            new meshStructure
            (
                mesh,
                uip,
                globalFaces_(),
                globalEdges_(),
                globalPoints_()
            )
        );
    }

    return *meshStructurePtr_;
}

const Foam::labelList& Foam::spanwiseAverage::patchAddressing(const primitivePatch& patch) const
{
    const edgeList& patchEdge = patch.edges();
    const pointField& patchPoint = patch.localPoints();
    const pointField& patchCentre = patch.faceCentres();
    
    scalar geoMin = min(patchPoint & averagePatchesDir_);
    patchSideEdgePtr_.reset(new edgeList());
    List<vector> patchSideEdgeCentre;
 
    forAll(patchEdge,edgei)
    {
        if (!patch.isInternalEdge(edgei))
	    {
            if (mag(patchEdge[edgei].unitVec(patchPoint) & averagePatchesDir_) < 1000*SMALL)
	        {
                vector edgeCentre = patchEdge[edgei].centre(patchPoint); 
        	    if (mag((edgeCentre & averagePatchesDir_) - geoMin) < 1000*SMALL)
	            {
    		        patchSideEdgePtr_().append(patchEdge[edgei]);
		            patchSideEdgeCentre.append(edgeCentre);
		        }
	        }	
	    }	
    }    

    patchFaceToEdgePtr_.reset(new labelList(patchCentre.size())); 
    forAll(patchCentre,facei)
    {
    	List<scalar> disToEdge = mag(patchCentre[facei] - patchSideEdgeCentre);
        label minIndex = findMin(disToEdge);
        patchFaceToEdgePtr_()[facei] = minIndex;
    }

    return *patchFaceToEdgePtr_;
}

const Foam::word Foam::spanwiseAverage::averageName 
(
    const word& fieldName
) const
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
    fieldsName_(dict.get<wordList>("fields")),
    includeInternal_(dict.subDict("Internal").getOrDefault<bool>("includeInternal",true)),
    includePatches_(dict.subDict("Patches").getOrDefault<bool>("includePatches",false)),
    sidePatchID_(),
    averagePatchesName_(), 
    averagePatchesDir_()
{
    if (includeInternal_)
    {
        sidePatchID_ = mesh.boundaryMesh().findPatchID(dict.subDict("Internal").get<word>("sidePatchName"));
    }
    
    if (includePatches_)
    {
        averagePatchesName_ = dict.subDict("Patches").get<wordList>("averagePatchesName");
        averagePatchesDir_ = dict.subDict("Patches").get<vector>("averagePatchesDir");
    }
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
