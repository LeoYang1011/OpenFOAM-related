/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    faceInterpolation

Description
    Face center to point interpolation. Included in faMesh.

SourceFiles
    faceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef faceInterpolation_H
#define faceInterpolation_H

#include "faMesh.H"
#include "scalarList.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceInterpolation Declaration
\*---------------------------------------------------------------------------*/
class faceInterpolation
{
    // Private data

        // Reference to faMesh
        const faMesh& faMesh_;

        // Demand-driven data

            //- Face-to-point weighting factors
            mutable scalarListList* faceToPointWeightsPtr_;

            //- Point-to-face weighting factors
            mutable scalarListList* pointToFaceWeightsPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        faceInterpolation(const faceInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const faceInterpolation&);

        //- Face-to-point weights
        const scalarListList& faceToPointWeights() const;

        //- Construct face-to-point weighting factors
        void makeFaceToPointWeights() const;

        //- Point-to-face weights
        const scalarListList& pointToFaceWeights() const;

        //- Construct Point-to-face weighting factors
        void makePointToFaceWeights() const;

        //- Clear all geometry and addressing
        void clearWeights();

public:

    // Constructors

        //- Construct given an faMesh
        faceInterpolation(const faMesh&);


    // Destructor

        ~faceInterpolation();


    // Member functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return faMesh_;
        }

        //- Interpolate from faces to points
        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Interpolate from faces to points
        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Do what is neccessary ifthe mesh has moved
        bool movePoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faceInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
