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
#include "IcoFOAMOutgoingField.H"
#include "IcoFOAMIncomingField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    CPLSocketFOAM socket;
    socket.initComms(argc, argv);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"
    socket.setOpenFOAM(runTime, mesh);
    socket.init();
    CPL::OutgoingFieldPool cnstPool(socket.cnstPortionRegion, socket.cnstRegion);
    CPL::IncomingFieldPool bcPool(socket.bcPortionRegion, socket.bcRegion); 
    // Define constrains and BCs
    (new StressOutgoingField("stresscnst", socket.cnstPortionRegion, socket.cnstRegion,U, nu, mesh, 1.0))->addToPool(cnstPool);
    (new VelIncomingField("velocitybc", socket.bcPortionRegion, socket.bcRegion, U, nu, mesh, 1.0))->addToPool(bcPool);

    // Setup of the Fields
    cnstPool.setupAll();
    bcPool.setupAll();
    // Allocation of internal send/receive buffers
    socket.allocateBuffers(cnstPool, bcPool);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	// Initial communication to initialize domains
    // cnstPool.packAll();
    // socket.send();
    // socket.receive();
    // bcPool.unpackAll();
    socket.communicate(cnstPool, bcPool);

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        socket.communicate(cnstPool, bcPool);
        // cnstPool.updateAll();
        // cnstPool.packAll();
        // socket.send();
        // socket.receive();
        // bcPool.updateAll();
        // bcPool.unpackAll();
        //
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
