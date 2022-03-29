/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.
    Consistent formulation without time-step and relaxation dependence by
    Jasak and Tukovic.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "CPLSocketFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Check if coupled based on cpl/COUPLER.in input file                                                                                                                                     
    bool coupled;                                                                                                                                                                             
    if (file_exists("./cpl/COUPLER.in")) {
    	Info<< "Assuming coupled run as cpl/COUPLER.in exists\n" << endl;
     	coupled=true;                                                                                                                                                                         
 	} else {                                                                                                                                                                                  
        Info<< "Assuming uncoupled run as cpl/COUPLER.in does not exist\n" << endl;
        coupled=false;                                                                                                                                                                        
    }                                                                                                                                                                                         
    // Create a CPL object (not used if uncoupled) and intialises MPI                                                                                                                         
    CPLSocketFOAM CPL;                                                                                                                                                                        
    if (coupled)                                                                                                                                                                              
        CPL.initComms(argc, argv);                                                                                                                                                            


#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

    // MPI_Init is called somewhere in the PStream library                                                                                                                                    
    if (coupled)                                                                                                                                                                              
        CPL.initCFD(runTime, mesh); 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	// Initial communication to initialize domains
    if (coupled){
        CPL.pack(U, p, nu, mesh, CPL.VEL);
        CPL.send();
        CPL.recvVelocity();
        CPL.unpackVelocity(U, mesh);
    //    CPL.recvVelocityPressure();
    //    CPL.unpackVelocityPressure(U, p, mesh);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		if (coupled){
            std::cout << "CPL.VEL is on " << std::endl;
            CPL.pack(U, p, nu, mesh, CPL.VEL);
            //CPL.pack(U, p, nu, mesh, CPL.STRESS);
            CPL.send();
            CPL.recvVelocity();
            CPL.unpackVelocity(U, mesh);
        //    CPL.recvVelocityPressure();
        //    CPL.unpackVelocityPressure(U, p, mesh);
        }

#       include "CourantNo.H"

        // Time-derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(U));

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(ddtUEqn + HUEqn == -fvc::grad(p));
        }

        // Prepare clean 1/a_p without time derivative contribution
        volScalarField rAU = 1.0/HUEqn.A();

        // --- PISO loop

        while (piso.correct())
        {
            // Calculate U from convection-diffusion matrix
            U = rAU*HUEqn.H();

            // Consistently calculate flux
            piso.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

            adjustPhi(phi, U, p);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        fvc::interpolate(rAU)/piso.aCoeff(U.name()),
                        p,
                        "laplacian(rAU," + p.name() + ')'
                    )
                 ==
                    fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Consistently reconstruct velocity after pressure equation
            piso.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

	CPL::finalize();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
