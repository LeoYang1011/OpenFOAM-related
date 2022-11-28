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

#include "volFields.H"
#include "meshStructure.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::spanwiseAverage::spanwiseAverageField
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fieldType* fldPtr = mesh_.cfindObject<fieldType>(fieldName);

    if (fldPtr)
    {
        const fieldType& fld = *fldPtr;

        const word resultName(averageName(fieldName));

        fieldType* resPtr = mesh_.getObjectPtr<fieldType>(resultName);

        if (!resPtr)
        {
            resPtr = new fieldType
            (
                IOobject
                (
                    resultName,
                    fld.mesh().time().timeName(),
                    fld.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fld
            );
            mesh_.objectRegistry::store(resPtr);
        }
        fieldType& res = *resPtr;

	    if (includeInternal_)
        {

            Info << "Averaging of the field: " << fld.name() << " on internal mesh." << endl; 
            
            const meshStructure& ms = meshAddressing(fld.mesh());

     	    if (globalFaces_().empty())
            {
                return;
            }

            const labelList& cellToPatchFace = ms.cellToPatchFaceAddressing();

            Field<Type> regionField(globalFaces_().size(), Zero);
            labelList regionCount(globalFaces_().size(), 0);

            forAll(cellToPatchFace, celli)
            {
                const label regioni = cellToPatchFace[celli];
                regionField[regioni] += fld[celli];
                regionCount[regioni]++;
            }

            forAll(regionField, regioni)
            {
                regionField[regioni] /= regionCount[regioni];
            }

            // And send result back
            forAll(cellToPatchFace, celli)
            {
                const label regioni = cellToPatchFace[celli];
                res[celli] = regionField[regioni];
            }

	        res.correctBoundaryConditions();
        }

	    if (includePatches_)
	    {
            forAll(fld.boundaryField(),patchi)
	        {

		        if (averagePatchesName_.found(fld.mesh().boundaryMesh()[patchi].name()))
		        {

                    Info << "Averaging of the field: " << fld.name() << " on the patch: " 
                         << fld.mesh().boundaryMesh()[patchi].name() << endl;

                    const labelList& faceToEdge = patchAddressing(fld.mesh().boundaryMesh()[patchi]);
                    const Field<Type>& fldPatch = fld.boundaryField()[patchi];
                    Field<Type>& resPatch = res.boundaryFieldRef()[patchi];            		                  

		            Field<Type> edgeField(patchSideEdgePtr_().size(), Zero);
                    labelList edgeCount(patchSideEdgePtr_().size(), 0);

		            forAll(faceToEdge, facei)
                    {
                        const label edgei = faceToEdge[facei];
                        edgeField[edgei] += fldPatch[facei];
                        edgeCount[edgei]++;
                    }

                    forAll(edgeField, edgei)
                    {
                        edgeField[edgei] /= edgeCount[edgei];
                    }
    
                    // And send result back
                    forAll(faceToEdge, facei)
                    {
                        const label edgei = faceToEdge[facei];
                        resPatch[facei] = edgeField[edgei];
                    }
		        }	
	        }	    
	    }
    }
}


// ************************************************************************* //
