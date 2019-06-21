/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "CPLSocketFOAM.H"
#include "StressOutgoingField.H"
#include "VelocityIncomingField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Check args first

    CPLSocketFOAM socket;
    socket.initComms(argc, argv);
    socket.loadParamFile();
    CPL::get_file_param("constrain-enabled", "", socket.sendEnabled);
    CPL::get_file_param("bc-enabled", "", socket.recvEnabled);
    // Currently density is not used here
    double density;
    CPL::get_file_param("initial-conditions", "density", density);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"

    // Conversion factors (Units come in S.I)
    double den_Si2den_Md, eta_Si2eta_Md;
    CPL::get_file_param("conversion-factors", "density", den_Si2den_Md);
    CPL::get_file_param("conversion-factors", "dyn-viscosity", eta_Si2eta_Md);
    density *= den_Si2den_Md;
    dimensionedScalar eta(density*nu*eta_Si2eta_Md);

    //Create stress field
    volSymmTensorField sigma
            (
                IOobject
                (
                    "sigma",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                eta*2*dev(symm(fvc::grad(U)))
            );

    #include "initContinuityErrs.H"
    socket.setOpenFOAM(runTime, mesh);
    socket.init();
    CPL::OutgoingFieldPool cnstPool(socket.cnstPortionRegion, socket.cnstRegion);
    CPL::IncomingFieldPool bcPool(socket.bcPortionRegion, socket.bcRegion); 

    
    if (socket.sendEnabled) {
        (new StressOutgoingField("stresscnst", socket.cnstPortionRegion, 
                                 socket.cnstRegion, sigma, mesh))->addToPool(&cnstPool);
        cnstPool.setupAll();
        if (!socket.sendBuffAllocated)
            socket.allocateSendBuffer(cnstPool);
    }
    if (socket.recvEnabled) {
        (new VelIncomingField("velocitybc", socket.bcPortionRegion, socket.bcRegion, 
                              U, mesh, density))->addToPool(&bcPool);
        bcPool.setupAll();
        if (!socket.recvBuffAllocated)
            socket.allocateRecvBuffer(bcPool);
    }

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	// Initial communication to initialize domains
    socket.communicate(cnstPool, bcPool);

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }
        
        // Update sigma value after solving for U
        //NOTE: The first row of cells in the y direction will not contain
        //the correct stress but it does not matter as that is the BC region.
        // On the other hand if this is used for plotting stress instead of fields
        // generated from stressComponents utility, be careful!
        sigma = eta*2*dev(symm(fvc::grad(U)));
        // Pack/unpack and communicate after computing fields but before writing to file.
        socket.communicate(cnstPool, bcPool);
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;
	CPL::finalize();

    return 0;
}


// ************************************************************************* //
