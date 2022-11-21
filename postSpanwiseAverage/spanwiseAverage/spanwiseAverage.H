/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | \\  /    A nd           | www.openfoam.com
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

Class
    Foam::functionObjects::spanwiseAverage

Group
    grpFieldFunctionObjects

Description
    Computes the arithmetic average of given quantities along columns of cells
    in a given direction for structured-like layered meshes. It is, for example,
    useful for channel-like cases where spanwise average of a field is desired.
    However, the \c spanwiseAverage function object does not operate on arbitrary
    unstructured meshes.

    For each patch face, calculates the average value of all cells attached in
    the patch face normal direction, and then pushes the average value back
    to all cells in the column.

    Operands:
    \table
      Operand           | Type               | Location
      input             | vol\<Type\>Field   | $FOAM_CASE/\<time\>/\<inpField\>
      output file       | -                  | -
      output field      | vol\<Type\>Field   | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

SourceFiles
    spanwiseAverage.C
    spanwiseAverageTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef spanwiseAverage_H
#define spanwiseAverage_H

#include "volFieldsFwd.H"
#include "polyMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;
class meshStructure;

/*---------------------------------------------------------------------------*\
                        Class spanwiseAverage Declaration
\*---------------------------------------------------------------------------*/

class spanwiseAverage
{
    // Private Data

        //- Patches on which to collapse the fields (in sorted order)
	const polyMesh& mesh_;

	label patchID_;

	wordList fieldsName_;

        autoPtr<globalIndex> globalFaces_;
        autoPtr<globalIndex> globalEdges_;
        autoPtr<globalIndex> globalPoints_;
        autoPtr<meshStructure> meshStructurePtr_;


    // Private Member Functions

        //- Create the column average field name
        const word averageName(const word& fieldName);

        //- Return the column-based addressing
        void meshAddressing();

        //- Calculate the averaged field and return true if successful
        template<class Type>
        void spanwiseAverageField(const word& fieldName);


public:

    // Constructors

        //- Construct from Time and dictionary
        spanwiseAverage
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

	//- No copy construct
        spanwiseAverage(const spanwiseAverage&) = delete;

        //- No copy assignment
        void operator=(const spanwiseAverage&) = delete;

    //- Destructor
    virtual ~spanwiseAverage() = default;

    // Member Functions

        const wordList& fields();

        //- Execute, currently does nothing
        void execute();

        //- Write the results
        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "spanwiseAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
