/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvIOoptionList.H"
#include "simpleControl.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
	#include "CourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"  //correct timestep
	
	volScalarField DEff("DEff", turbulence->nu()/Sc + turbulence->nut()/Sct);

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(C)
              + fvm::div(phi, C)
              - fvm::laplacian(DEff, C)
             ==
                fvOptions(C)
            );
        }

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(Cp)
              + fvm::div(phi, Cp)
              - fvm::laplacian(DEff, Cp)
	     ==
                fvOptions(Cp)
            );
        }
//if (oCorr == nOuterCorr-1)
//{
	forAll(C,celli)
	{
		if(C[celli] > Cs[celli])	//is Cs>Ccurrent?
			{
			Cp[celli] = Cp[celli] + ( C[celli] - Cs[celli] ); //if yes, lets define precipitation from Ccurrent
			C[celli] = Cs[celli];				  //C goes without precipitation
			}	
	}
	
	forAll(Cp,celli)
	{
		if( ( Cp[celli] > 0) & (C[celli] < Cs[celli]) )
			{
			C[celli] = C[celli] + Cp[celli];
			Cp[celli] = 0;
				
				if (C[celli] > Cs[celli])
					{
						Cp[celli] = (C[celli] - Cs[celli]);	
						C[celli] = Cs[celli];	
					}						  
			}	
	}
//}
////////////////////////////////////////////////////////
    scalar DiNum = 0.0;
    scalar meanDiNum = 0.0;

    surfaceScalarField DeffInterpolated
    (
        fvc::interpolate(DEff)
      / mesh.surfaceInterpolation::deltaCoeffs()
    );

    DiNum = gMax(DeffInterpolated.internalField())*runTime.deltaT().value();

    meanDiNum = (average(DeffInterpolated)).value()*runTime.deltaT().value();

    Info<< "Diffusion Number mean: " << meanDiNum << " max: " << DiNum << endl; //    return DiNum;

////////////////////////////////////////////////////////

Info<< "Calculating kc and Sh" << endl;

label patchi = mesh.boundaryMesh().findPatchID("lowerWall");

    surfaceScalarField CpInterpolated
    (
        fvc::interpolate(Cp)
      / mesh.surfaceInterpolation::deltaCoeffs()
    );

// kc = 1 / Tb * D * (dT/dy)_wall
kc.boundaryField()[patchi] = 1/(average(CpInterpolated)).value()*(average(DeffInterpolated)).value()*-Cp.boundaryField()[patchi].snGrad();

// Calculates average kc

		scalar area = gSum(mesh.magSf().boundaryField()[patchi]);
                scalar sumField = 0;

                if (area > 0)
                {
                    sumField = gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                      * kc.boundaryField()[patchi]
                    )/ area;
                }


                Info<< "    Average kc over patch heater = "   << sumField << endl;

           volScalarField Sh_new
           (
                IOobject
                (
                    "Sh_new",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
// Sh = kc * d / DT
		kc*d/(average(DeffInterpolated)).value()
           );

/////////////////////////////////////////////////

	if(runTime.outputTime())
	{
		C.write();
		DEff.write();
    	}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
