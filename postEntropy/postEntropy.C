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
    postSpanwiseAverage

Description
    Performs averaging of a field in one homogeneous direction (e.g. spanwise
    direction of an airfoil). Is only applicable to structered grid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    word outputDir = "./postProcessing/enstrophies/";

    if(!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    // Setup channel indexing for averaging over channel down to a line

    // For each time step read all fields
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

	Info<< "Collapsing fields for time " << runTime.timeName() << endl;

        volVectorField vorticity
        (
            IOobject
            (
                "vorticity",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
       );

       volScalarField enstrophiesx
       (
            IOobject
            (
                "enstrophiesx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("enstrophiesx",dimless,0.0)
       );

       volScalarField enstrophiesy
       (
            IOobject
            (
                "enstrophiesy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("enstrophiesy",dimless,0.0)
       );

       volScalarField enstrophiesz
       (
            IOobject
            (
                "enstrophiesz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("enstrophiesz",dimless,0.0)
       );

       volScalarField enstrophiesp
       (
            IOobject
            (
                "enstrophiesp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("enstrophiesp",dimless,0.0)
       );


       const scalarField& V = mesh.V();

       scalar enstrophiesxSum = 0.0; 
       scalar enstrophiesySum = 0.0;  
       scalar enstrophieszSum = 0.0;  
       scalar enstrophiespSum = 0.0;  

       forAll(V,celli)
       {
           enstrophiesx[celli] = 0.5*Foam::sqr(vorticity[celli].x())*V[celli];
           enstrophiesy[celli] = 0.5*Foam::sqr(vorticity[celli].y())*V[celli];
           enstrophiesz[celli] = 0.5*Foam::sqr(vorticity[celli].z())*V[celli];
           enstrophiesp[celli] = 0.5*(Foam::sqr(vorticity[celli].x())+Foam::sqr(vorticity[celli].y()))*V[celli];

           enstrophiesxSum += enstrophiesx[celli]; 
           enstrophiesySum += enstrophiesy[celli]; 
           enstrophieszSum += enstrophiesz[celli]; 
           enstrophiespSum += enstrophiesp[celli]; 
       }  

       fileName writeFile = ("postProcessing/enstrophies/"+std::to_string(timeDirs[timeI].value()));
       OFstream OS(writeFile);

       Info<< "Writing the intergation of enstrophies to file: " <<  writeFile  << endl;

       OS << "enstrophiesxSum: " << enstrophiesxSum << nl;
       OS << "enstrophiesySum: " << enstrophiesySum << nl;
       OS << "enstrophieszSum: " << enstrophieszSum << nl;
       OS << "enstrophiespSum: " << enstrophiespSum << nl;

       OS << endl;     
 
       enstrophiesx.write();
       enstrophiesy.write();
       enstrophiesz.write();
       enstrophiesp.write();

       runTime.write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
