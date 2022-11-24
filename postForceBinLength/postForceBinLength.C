/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

Application
    postForceBinLength

Description
    Get the exact length of each bin, and the position of each face contained 
    in each bin. Valid for patches parallel to the x or y or z axes, skewed 
    patches have not been tested yet.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
#   include "createMesh.H"

    // Setup channel indexing for averaging over channel down to a line

    IOdictionary forceBinDict
    (
        IOobject
        (
            "postForceBinDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    label patchId = mesh.boundaryMesh().findPatchID(forceBinDict.get<word>("patchName"));

    label binNum = forceBinDict.get<label>("bins");
    scalar start = forceBinDict.get<scalar>("start");
    scalar delta = forceBinDict.get<scalar>("delta");
    vector binDir = forceBinDict.get<vector>("binDir");
    string writePath = forceBinDict.get<string>("writePath");

    const polyPatch& forcePatch = mesh.boundaryMesh()[patchId]; 
    const vectorField& patchFacesCentre = forcePatch.faceCentres();
    const vectorField& patchPoints = forcePatch.localPoints();
    const faceList& patchFaces = forcePatch.localFaces();

    scalarField facesCentreDis = patchFacesCentre & binDir;
    scalarField pointsDis = patchPoints & binDir;
    scalarField relDis = facesCentreDis - start;

    List<scalar> disMin(binNum,GREAT);
    List<scalar> disMax(binNum,-GREAT);
    List<List<scalar>> binFaceCentre(binNum);

    forAll(patchFaces,facei)
    {
        label bini = min(max(floor(relDis[facei]/delta), 0), binNum - 1); 

	// find face center for each bin
	scalar fCentreDis = facesCentreDis[facei];
        List<scalar>& binFC = binFaceCentre[bini];
	if (binFC.empty())
	{
	    binFC.append(fCentreDis);
	}
	else
        {
            List<scalar> relBinFC = binFC - fCentreDis;
	    forAll(relBinFC,reli)
            {
	        relBinFC[reli] = fabs(relBinFC[reli]); 
	    }
            if (min(relBinFC) > 0.001)
            {
	        binFC.append(fCentreDis);
	    }
	}
        
	// find the max and min point distance
	face pFace = patchFaces[facei];
        List<scalar> faceDis(pFace.size());
        forAll(pFace,pointi)
	{
            faceDis[pointi] = pointsDis[pFace[pointi]]; 
	}
        scalar faceMin = min(faceDis);	
        scalar faceMax = max(faceDis);
        disMin[bini] = min(faceMin,disMin[bini]);	
        disMax[bini] = max(faceMax,disMax[bini]);	
    }
   
    List<scalar> lengthBin = disMax - disMin;

    fileName writeFile = (writePath+"/forceBinLocate");
    OFstream OS(writeFile);

    Info<< "Writing the bin length and the bin faces center to file: " <<  writeFile  << endl;

    OS << "#" << " The first column is the length of each bin." << nl; 
    OS << "#" << " The data after the second column is the position of the faces centre in each bin." << nl; 

    forAll(lengthBin,bini)
    {
        OS << lengthBin[bini];
        forAll(binFaceCentre[bini],facei)
	{
            OS << tab << binFaceCentre[bini][facei];
	}	
        OS << nl;
    }
    OS << endl;

    // For each time step read all fields
    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
