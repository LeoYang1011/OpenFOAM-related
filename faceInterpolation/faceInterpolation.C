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

Description
    Face to edge interpolation scheme. Included in faMesh.

\*---------------------------------------------------------------------------*/

#include "faceInterpolation.H"
#include "faceList.H"
#include "demandDrivenData.H"
#include "PstreamCombineReduceOps.H"
#include "processorFaPatch.H"
#include "processorFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
const scalarListList& faceInterpolation::faceToPointWeights() const
{
    if (!faceToPointWeightsPtr_)
    {
        makeFaceToPointWeights();
    }

    return *faceToPointWeightsPtr_;
}

void faceInterpolation::makeFaceToPointWeights() const
{
 
    const pointField& points = mesh().points();
    const faceList& faces = mesh().faces();
    const labelListList& pointFaces = mesh().patch().pointFaces();

    faceToPointWeightsPtr_ = new scalarListList(points.size());
    scalarListList& faceToPointWeights = *faceToPointWeightsPtr_;
    scalarField sumw(points.size(),0.0);

    forAll(pointFaces,pointI)
    {
        const labelList& curFaces = pointFaces[pointI];

        scalarList& pw = faceToPointWeights[pointI];
        pw.setSize(curFaces.size());

        forAll(curFaces,faceI)
        {
            pw[faceI] = 
                1.0/mag(faces[curFaces[faceI]].centre(points) - points[pointI]);
            sumw[pointI] += pw[faceI];
        }
    }

    for(const faPatch& fap : mesh().boundary())
    {
        if(Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            labelList patchPointLabels = procPatch.pointLabels(); 
   
            scalarField patchSumw
            (
                patchPointLabels.size(),
                0.0
            );

            forAll (patchSumw, pointI)
            {
                patchSumw[pointI] =
                    sumw[patchPointLabels[pointI]];
            }

            {
                OPstream::write
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(patchSumw.begin()),
                    patchSumw.byteSize()
                );
            }

            scalarField ngbPatchSumw
            (
                procPatch.neighbPoints().size(),
                0.0
            );

            {
                IPstream::read
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>(ngbPatchSumw.begin()),
                    ngbPatchSumw.byteSize()
                );
            }

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll (nonGlobalPatchPoints, pointI)
            {
                sumw[patchPointLabels[nonGlobalPatchPoints[pointI]]] +=
                    ngbPatchSumw
                    [
                        procPatch.neighbPoints()[nonGlobalPatchPoints[pointI]]
                    ];
            }
        }
    }

    if (mesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            mesh().globalData().sharedPointLabels();

        scalarField spSum(spLabels.size(), 0.0);
        forAll (spSum, pointI)
        {
            spSum[pointI] = sumw[spLabels[pointI]];
        }

        const labelList& addr = mesh().globalData().sharedPointAddr();

        scalarField gpSum
        (
            mesh().globalData().nGlobalPoints(),
            0.0
        );

        forAll (addr, i)
        {
            gpSum[addr[i]] += spSum[i];
        }

        combineReduce(gpSum, plusEqOp<scalarField>());

        // Extract local data
        forAll (addr, i)
        {
            spSum[i] = gpSum[addr[i]];
        }

        forAll (spSum, pointI)
        {
            sumw[spLabels[pointI]] = spSum[pointI];
        }
    }


    forAll(pointFaces,pointI)
    {
        const labelList& curFaces = pointFaces[pointI];

        scalarList& pw = faceToPointWeights[pointI];

        forAll(curFaces,faceI)
        {
            pw[faceI] /= sumw[pointI];
        }
    }

}

const scalarListList& faceInterpolation::pointToFaceWeights() const
{
    if (!pointToFaceWeightsPtr_)
    {
        makePointToFaceWeights();
    }

    return *pointToFaceWeightsPtr_;
}

void faceInterpolation::makePointToFaceWeights() const
{
 
    const pointField& points = mesh().points();
    const faceList& faces = mesh().faces();

    pointToFaceWeightsPtr_ = new scalarListList(faces.size());
    scalarListList& pointToFaceWeights = *pointToFaceWeightsPtr_;
    scalarField sumw(faces.size(),0.0);

    forAll(faces,faceI)
    {
        const labelList& curPoints = faces[faceI];

        scalarList& fw = pointToFaceWeights[faceI];
        fw.setSize(curPoints.size());

        forAll(curPoints,pointI)
        {
            fw[pointI] = 
                1.0/mag(faces[faceI].centre(points) - points[curPoints[pointI]]);
            sumw[faceI] += fw[pointI];
        }

        forAll(curPoints,pointI)
        {
            fw[pointI] /= sumw[faceI]; 
        }
    }
}

void faceInterpolation::clearWeights()
{
    deleteDemandDrivenData(faceToPointWeightsPtr_);
    deleteDemandDrivenData(pointToFaceWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

faceInterpolation::faceInterpolation(const faMesh& fam)
:
    faMesh_(fam),
    faceToPointWeightsPtr_(NULL),
    pointToFaceWeightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //
faceInterpolation::~faceInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
tmp<Field<Type> > faceInterpolation::faceToPointInterpolate
(
    const Field<Type>& ff
) const
{
    // Check size of the given field
    if (ff.size() != mesh().patch().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > faceInterpolation::"
            "faceToPointInterpolate(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << mesh().patch().size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            mesh().nPoints(), Zero
        )
    );

    Field<Type>& result = tresult.ref();

    const labelListList& pointFaces = mesh().patch().pointFaces();
    const scalarListList& weights = faceToPointWeights();

    forAll(pointFaces, pointI)
    {
        const labelList& curFaces = pointFaces[pointI];
        const scalarList& w = weights[pointI];

        forAll(curFaces, faceI)
        {
            result[pointI] += w[faceI]*ff[curFaces[faceI]];
        }
    }

    for(const faPatch& fap : mesh().boundary())
    {
        if(Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            labelList patchPointLabels = procPatch.pointLabels(); 
   
            Field<Type> patchValue
            (
                patchPointLabels.size(),
                Zero
            );

            forAll (patchValue, pointI)
            {
                patchValue[pointI] =
                    result[patchPointLabels[pointI]];
            }
 
            {
                OPstream::write
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(patchValue.begin()),
                    patchValue.byteSize()
                );
            }

            Field<Type> ngbPatchValue
            (
                procPatch.neighbPoints().size(),
                Zero
            );

            {
                IPstream::read
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>(ngbPatchValue.begin()),
                    ngbPatchValue.byteSize()
                );
            }
 
            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll (nonGlobalPatchPoints, pointI)
            {
                result[patchPointLabels[nonGlobalPatchPoints[pointI]]] +=
                    ngbPatchValue
                    [
                        procPatch.neighbPoints()[nonGlobalPatchPoints[pointI]]
                    ];
            }
        }
    }

    if (mesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            mesh().globalData().sharedPointLabels();

        Field<Type> spResult(spLabels.size(), Zero);
        forAll (spResult, pointI)
        {
            spResult[pointI] = result[spLabels[pointI]];
        }

        const labelList& addr = mesh().globalData().sharedPointAddr();

        Field<Type> gpResult
        (
            mesh().globalData().nGlobalPoints(),
            Zero
        );

        forAll (addr, i)
        {
            gpResult[addr[i]] += spResult[i];
        }

        combineReduce(gpResult, plusEqOp<Field<Type> >());

        // Extract local data
        forAll (addr, i)
        {
            spResult[i] = gpResult[addr[i]];
        }

        forAll (spResult, pointI)
        {
            result[spLabels[pointI]] = spResult[pointI];
        }
    }

    return tresult;
}

template<class Type>
tmp<Field<Type> > faceInterpolation::faceToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = faceToPointInterpolate(tff());
    tff.clear();
    return tint;
}

template<class Type>
tmp<Field<Type> > faceInterpolation::pointToFaceInterpolate
(
    const Field<Type>& pf
) const
{
    // Check size of the given field
    if (pf.size() != mesh().patch().nPoints())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > faceInterpolation::"
            "pointToFaceInterpolate(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << mesh().patch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            mesh().nFaces(), Zero
        )
    );

    Field<Type>& result = tresult.ref();

    const faceList& faces = mesh().faces();
    const scalarListList& weights = pointToFaceWeights();

    forAll(faces, faceI)
    {
        const labelList& curPoints = faces[faceI];
        const scalarList& w = weights[faceI];

        forAll(curPoints, pointI)
        {
            result[faceI] += w[pointI]*pf[curPoints[pointI]];
        }
    }

    return tresult;
}

template<class Type>
tmp<Field<Type> > faceInterpolation::pointToFaceInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = faceToPointInterpolate(tff());
    tff.clear();
    return tint;
}

// Do what is neccessary if the mesh has moved
bool faceInterpolation::movePoints() const
{
    deleteDemandDrivenData(faceToPointWeightsPtr_);
    deleteDemandDrivenData(pointToFaceWeightsPtr_);
 
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
