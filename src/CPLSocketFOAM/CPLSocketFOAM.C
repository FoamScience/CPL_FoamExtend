/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________        
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________       
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________      
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________     
        _\/\\\_____________\/\\\/////////____\/\\\_____________    
         _\//\\\____________\/\\\_____________\/\\\_____________   
          __\///\\\__________\/\\\_____________\/\\\_____________  
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_ 
            _______\/////////__\///______________\///////////////__


                         C P L  -  L I B R A R Y 

           Copyright (C) 2012-2017 Edward Smith, Eduardo Ramos-Fernandez  & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    See CPLSocketFOAM.H 

*/

#include "blockMesh.H"
#include <sstream>
#include <unistd.h>
#include <bitset>
#include "interpolation.H"
#include "PstreamGlobals.H" 
#include "CPLSocketFOAM.H"

#include "cpl/CPL_cartCreate.h"

// Initialise CFD realm communicator
void CPLSocketFOAM::initComms (int& argc, char**& argv) {

    Foam::Info << "CPLSocketFOAM::initComms "<< Foam::endl;

    int flag = 0;
    int ierr = MPI_Initialized(&flag);

    Foam::Info << "MPI_Initialized(&flag) " << flag << Foam::endl;
    if (flag == 0)
		MPI_Init(&argc, &argv);
    CPL::init (CPL::cfd_realm, realmComm);
    ::Foam::PstreamGlobals::MPI_COMM_FOAM = realmComm;
    MPI_Comm_rank (realmComm, &rankRealm);

}

void CPLSocketFOAM::finalize() {

   CPL::finalize();

}

// Analyse mesh topology and perform CFD-side CPL_init.
void CPLSocketFOAM::initCFD(const Foam::Time &runTime, const Foam::fvMesh &mesh) {

	meshSearcher = new Foam::meshSearch(mesh);

    Foam::Info << "CPLSocketFOAM: Analysing processor and mesh topology"
               << Foam::endl;

    // Read from decomposePar dictionary the number of processors in each
    // direction. Must be decomposed with the "simple" method.
	//TODO: THis should not be necessary if the cell grid is defined in the config
	//file for the overlap region
    Foam::IOdictionary decomposeDict(Foam::IOobject ("decomposeParDict", 
									 runTime.time().system(), runTime,
                        			 IOobject::MUST_READ, IOobject::NO_WRITE,
									 false));
    Foam::dictionary simpleCoeffs = decomposeDict.subDict("simpleCoeffs");
    Foam::Vector<int> np = simpleCoeffs.lookup("n");
    nprocs = np.x() * np.y() * np.z();

    // Define arrays needed by MPI & CPL cart create routines 
    npxyz[0] = np.x();
    npxyz[1] = np.y();
    npxyz[2] = np.z();
    periods[0] = 1;
    periods[1] = 0;
    periods[2] = 1;

    Foam::Info << "CPLSocketFOAM: Defining new MPI Cartesian communicator"
               << Foam::endl;

    // Create custom cartesian communicator (cartComm) based on myCoords	
	// Assume xyz ordering (Simple decomposition)
	myCoords[2] = Pstream::myProcNo() / (np.x()*np.y());
	int mod_aux = Pstream::myProcNo() % (np.x()*np.y());
	myCoords[1] = mod_aux / np.x();
	myCoords[0] = mod_aux % np.x();
	
    CPL::Cart_create (realmComm, 3, npxyz, periods, 
                      myCoords.data(), &cartComm);

    MPI_Comm_rank (cartComm, &rankCart);


    // Prepare inputs for CPL::cfd_init 
    double dt_cfd = runTime.deltaTValue();
    double st = runTime.startTime().value();
    double et = runTime.endTime().value();
    if (((et-st)/dt_cfd)>1e8){
        FatalErrorIn("CPLSocketFOAM::initCFD()")
            << " Number of steps " << (et-st)/dt_cfd
            << " too large for single integer, Aborting." << exit(FatalError);
    } 
    int nsteps = nint((et-st)/dt_cfd);

    Foam::IOdictionary blockMeshDict(Foam::IOobject ("polyMesh/blockMeshDict", 
									 runTime.time().constant(), runTime,
                        			 IOobject::MUST_READ, 
									 IOobject::NO_WRITE, false));

	Foam::List<Foam::Vector<double>> vertices(blockMeshDict.lookup("vertices"));
    Foam::scalar convertToMeters(readScalar(blockMeshDict.lookup("convertToMeters")));

    //Read optional interpolate boundary flag from blockmesh
	interp_BC.readIfPresent("interp_BC", blockMeshDict);
    Foam::Info << "interp_BC" << interp_BC << Foam::endl;
    //if (interp_BC == false)
    //    FatalErrorIn("CPLSocketFOAM::initCFD()") << exit(FatalError);
    
    //Foam::dictionary boundary = blockMeshDict.subDict("boundary");
    //Foam::dictionary CPLReceiveMD = boundary.subDict("CPLReceiveMD");
    //CPLReceiveMD.readIfPresent("interp_BC");


    // Domain dimensions
    xyzL[0] = (vertices[1][0] - vertices[0][0])*convertToMeters;
    xyzL[1] = (vertices[3][1] - vertices[0][1])*convertToMeters;
    xyzL[2] = (vertices[4][2] - vertices[0][2])*convertToMeters;  
  
    double dummyDensity = -666.0;

    Foam::word dummyRegionName("dummy");
    Foam::blockMesh blocks(blockMeshDict, dummyRegionName);
    Foam::Vector<int> meshDensity = blocks[0].meshDensity();
   
    // Global number of cells
    ncxyz[0] = meshDensity.x();
    ncxyz[1] = meshDensity.y();
    ncxyz[2] = meshDensity.z();

    // Origin of the domain
    //double xyz_orig[3] = {0.0, 0.0, 0.0};
    double xyz_orig[3] = {vertices[0][0],vertices[0][1], vertices[0][2]}; //Hanqiao-Jun. 1 2020

    // Initialise CPL library
	//CPL::set_timing(0, nsteps, dt_cfd);
    CPL::setup_cfd (cartComm, xyzL, xyz_orig, ncxyz);
    getCellTopology();

    // Store some values from CPL that are useful later
    CPLDensity = CPL::density_cfd();

    Foam::Info << "OpenFOAM CPL topology initialisation complete" << Foam::endl;
    return;

}

//void CPLSocketFOAM::createFields(const Foam::Time &runTime, 
//                                 const Foam::fvMesh &mesh)
//{

//    dPdr = std::make_shared<Foam::volVectorField>
//    (
////    std::make_shared Foam::volVectorField dPdr
////    (
//        Foam::IOobject
//        (
//            "dPdr",
//            runTime.timeName(),
//            mesh,
//            Foam::IOobject::NO_READ,
//            Foam::IOobject::AUTO_WRITE
//        ),
//        mesh
//    );

//}

void CPLSocketFOAM::getCellTopology() {                                                                                                                                                                                                                       
    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
   
    // Get overlap extents
    CPL::get_olap_limits(olapRegion.data());
    
    // Processor cell bounds for the overlap region
    CPL::my_proc_portion(olapRegion.data(), olapPortion.data());
    CPL::get_no_cells(olapPortion.data(), olapCells);

    // Processor cell bounds for velocity BCs region
    CPL::get_bnry_limits(velBCRegion.data());
    CPL::my_proc_portion(velBCRegion.data(), velBCPortion.data());
    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Processor cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion(cnstFRegion.data(), cnstFPortion.data());
    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);
}

void CPLSocketFOAM::allocateBuffers(int sendtype) {

    //Check what is to be packed and sent
    int packsize=0;
    if ((sendtype & VEL) == VEL){
        packsize += VELSIZE;
    }
    if ((sendtype & PRESSURE) == PRESSURE){
        packsize += PRESSURESIZE;
    }
    if ((sendtype & GRADPRESSURE) == GRADPRESSURE){
        packsize += GRADPRESSURESIZE;
    }
    if ((sendtype & STRESS) == STRESS){
        packsize += STRESSSIZE;
    }
    if ((sendtype & DIVSTRESS) == DIVSTRESS){
        packsize += DIVSTRESSSIZE;
    }

    // Components for every local cell
    int sendShape[4] = {packsize, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
    sendBuf.resize(4, sendShape);

    if (sendtype > 31){
        FatalErrorIn
        (
            "CPLSocketFOAM::allocateBuffers()"
        )
            << " sendtype bit flag unknown type "
            << sendtype << " Aborting." << exit(FatalError);
    }
}  
    
// Packs the 9 components of the stress-tensor to 
// the socket's CPL::ndArray storage.
void CPLSocketFOAM::pack(volVectorField &U, 
                         volScalarField &p, 
                         dimensionedScalar &nu, 
                         fvMesh &mesh, 
                         int sendtype)
{

    // It would be nice to write these CPL fields directly here
    // but write interval is not accessible from Foam::Time object
    // So instead I've added all of them to fields

    // Evaluate gradient of pressure
    //if ((sendtype & GRADPRESSURE) == GRADPRESSURE)
        Foam::volVectorField gradP("gradP", fvc::grad(p));

    // Evaluate the stress tensor sigma at all local cells
    //if (((sendtype & STRESS) == STRESS) | 
    //    ((sendtype & DIVSTRESS) == DIVSTRESS)) {
        Foam::dimensionedScalar mu(CPLDensity*nu);
	    Foam::volSymmTensorField sigma(nu*2*dev(symm(fvc::grad(U))));
    //}

    // Evaluate divergence of sigma
    //if ((sendtype & DIVSTRESS) == DIVSTRESS)
        Foam::volVectorField divsigma(fvc::div(sigma));

    //Foam::volTensorField Ei = tensor::one;

    //Reallocate buffer depending on send type
    allocateBuffers(sendtype);
	if (CPL::is_proc_inside(cnstFPortion.data())) {

		// Loop over socket cells
		Foam::label cell;
		Foam::point globalPos;
		int glob_cell[3], loc_cell[3];
		for (int ix = cnstFPortion[0]; ix <= cnstFPortion[1]; ix++) {
			for (int iy = cnstFPortion[2]; iy <= cnstFPortion[3]; iy++) {
				for (int iz = cnstFPortion[4]; iz <= cnstFPortion[5]; iz++) {
					// Global position at cell center
					glob_cell[0] = ix;
					glob_cell[1] = iy;
					glob_cell[2] = iz;
					CPL::map_glob2loc_cell(cnstFPortion.data(), glob_cell, loc_cell);
                    // Minus one here as constraint is convention of starting at 1
					globalPos = Foam::point((glob_cell[0] + 0.5) * dx,
											(glob_cell[1] + 0.5) * dy, 
											(glob_cell[2] + 0.5) * dz);

                    //Foam::Info << "CPLSocketFOAM globalPos " << ix << " " << iy 
                    //           << " " << iz << " " << globalPos << Foam::endl;
					cell = meshSearcher->findNearestCell(globalPos);

					if (cell != -1) {
                        int npack = 0;
                        if ((sendtype & VEL) == VEL)
				        {
				            // Get value of velocity 3D
				            sendBuf(npack+0,loc_cell[0],loc_cell[1],loc_cell[2]) = U[cell].x();
				            sendBuf(npack+1,loc_cell[0],loc_cell[1],loc_cell[2]) = U[cell].y();
				            sendBuf(npack+2,loc_cell[0],loc_cell[1],loc_cell[2]) = U[cell].z();
                            npack += VELSIZE;

#if DEBUG
			                Foam::Info << "CPLSocketFOAM::pack vel " << ix << " " << iy << " " << iz
                                       << " " << cell << " " << " " << globalPos << " " <<
                                       sendBuf(npack-3,loc_cell[0],loc_cell[1],loc_cell[2]) << " " <<
                                       sendBuf(npack-2,loc_cell[0],loc_cell[1],loc_cell[2]) << " " <<
                                       sendBuf(npack-1,loc_cell[0],loc_cell[1],loc_cell[2]) << " " << Foam::endl;
#endif
				        }

                        if ((sendtype & PRESSURE) == PRESSURE)
				        {
				            // Store value of gradient pressure
				            sendBuf(npack+0,loc_cell[0],loc_cell[1],loc_cell[2]) = p[cell];
                            npack += PRESSURESIZE;
				        }

                        if ((sendtype & GRADPRESSURE) == GRADPRESSURE)
				        {
				            // Store value of gradient pressure
				            sendBuf(npack+0,loc_cell[0],loc_cell[1],loc_cell[2]) = gradP[cell].x();
				            sendBuf(npack+1,loc_cell[0],loc_cell[1],loc_cell[2]) = gradP[cell].y();
				            sendBuf(npack+2,loc_cell[0],loc_cell[1],loc_cell[2]) = gradP[cell].z();
                            npack += GRADPRESSURESIZE;
				        }

                        if ((sendtype & STRESS) == STRESS)
				        {
				            // Get value of stress 9D by interpolating
				            sendBuf(npack+0,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xx();
				            sendBuf(npack+1,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xy();
				            sendBuf(npack+2,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xz();
				            sendBuf(npack+3,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xy();
				            sendBuf(npack+4,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].yy();
				            sendBuf(npack+5,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].yz();
				            sendBuf(npack+6,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xz();
				            sendBuf(npack+7,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].yz();
				            sendBuf(npack+8,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].zz();
                            npack += STRESSSIZE;
				        }

                        if ((sendtype & DIVSTRESS) == DIVSTRESS)
				        {
				            // Get value of velocity 3D
				            sendBuf(npack+0,loc_cell[0],loc_cell[1],loc_cell[2]) = divsigma[cell].x();
				            sendBuf(npack+1,loc_cell[0],loc_cell[1],loc_cell[2]) = divsigma[cell].y();
				            sendBuf(npack+2,loc_cell[0],loc_cell[1],loc_cell[2]) = divsigma[cell].z();
                            npack += DIVSTRESSSIZE;

				        }

				        if (sendtype == DEBUG)
				        {
				            sendBuf(0,loc_cell[0],loc_cell[1],loc_cell[2]) = globalPos[0];
				            sendBuf(1,loc_cell[0],loc_cell[1],loc_cell[2]) = globalPos[1];
				            sendBuf(2,loc_cell[0],loc_cell[1],loc_cell[2]) = globalPos[2];
                            npack = 3;
				        }

                        //Check we have packed array correctly
                        if (npack != sendBuf.shapeData()[0]){ 
		                    FatalErrorIn ( "CPLSocketFOAM::pack()")
		                        << " sendtype mask error " << sendtype << ". "
		                           " Aborting."
		                        << exit(FatalError);
                        }
					}
					else
						std::cerr << "Warning: The point (" << (glob_cell[0]+0.5)*dx << "," << 
									(glob_cell[1]+0.5)*dy << "," << (glob_cell[2]+0.5)*dz << 
									") is outside the mesh for whatever reason! - Cell: " << 
									cell << std::endl;

				}
			}
		}
	}
}

// Unpacks the 3 components of the velocity-tensor from the socket's
// recvVelocity (cpl::ndArray) storage into a boundary condition.
double CPLSocketFOAM::unpackVelocity(volVectorField &U, fvMesh &mesh) 
{

//    Foam::Info << " unpackVelocity " << rankRealm << " " << 
//                CPL::is_proc_inside(velBCPortion.data()) << " " << 
//                velBCPortion.data() << Foam::endl;

    //bool interp_BC = false;

	if (CPL::is_proc_inside(velBCPortion.data())) {

		// TODO: Make this a utility general function that can be used on buffers
		if (CPL::get<int>("cpl_cfd_bc_slice")) {

			Foam::Info << "CPL_CFD_BC_SLICE is on: averaging CFD recvVelocity "
						  "in the x-z plane" << Foam::endl;

			// Number of cells in the local processor in x-z plane
			int N = recvVelocityBuff.shape(1) * recvVelocityBuff.shape(3);

			// For every component and y-value 
			for (int j = 0; j < recvVelocityBuff.shape(2); ++j) {
				for (int c = 0; c < recvVelocityBuff.shape(0); ++c) {
				// Sum across the x-z plane 
				double total = 0.0;
				for (int k = 0; k < recvVelocityBuff.shape(3); ++k)
					for (int i = 0; i < recvVelocityBuff.shape(1); ++i)
						total += recvVelocityBuff(c, i, j, k);

			// Find mean by dividing sum by number of cells 
				for (int k = 0; k < recvVelocityBuff.shape(3); ++k)
					for (int i = 0; i < recvVelocityBuff.shape(1); ++i)
						recvVelocityBuff(c, i, j, k) = total / static_cast<double> (N);
				}
			}
		}

		// Apply BCs only in certain directions.
		int applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
		int applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
		int applyBCz = CPL::get<int> ("cpl_cfd_bc_z");

		// Patch receiving B.Cs
		Foam::word receivePatchName = "CPLReceiveMD";
		Foam::label rvPatchID = mesh.boundaryMesh().findPatchID(receivePatchName);

		if (rvPatchID == -1) {
			FatalErrorIn ( "CPLSocketFOAM::unpack()")
				<< " Could not find patch ID " << receivePatchName << ". "
				   " Aborting."
				<< exit(FatalError);
		}

		Foam::fvPatchVectorField& rvPatch = U.boundaryField()[rvPatchID];
		const Foam::vectorField BoundaryfaceCntr = mesh.boundary()[rvPatchID].Cf();

		Foam::label cell;
		Foam::point closestCellCentre;
		for (int faceI = 0; faceI != BoundaryfaceCntr.size(); ++faceI) {
			double facex = BoundaryfaceCntr[faceI].x();
			double facey = BoundaryfaceCntr[faceI].y();
			double facez = BoundaryfaceCntr[faceI].z();
			// Find the cell indices for this position recvVelocity(:, ix, iy, iz)
        	int glob_cell[3]; int loc_cell[3];
			CPL::map_coord2cell(facex, facey, facez, glob_cell);
            //glob_cell[1] += 1; // Add one as boundary outside overlap by construction
	        bool valid_cell = CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

//            Foam::Info << "recvBuf " << facex << " " << facey << " " << facez << " "
//                        << glob_cell[0] << " " << glob_cell[1] << " " << glob_cell[2] << " "
//                        << loc_cell[0] << " " << loc_cell[1] << " " << loc_cell[2] << " "
//                        << valid_cell << Foam::endl;

            if (valid_cell) {

                double m = recvVelocityBuff(3, loc_cell[0], loc_cell[1], loc_cell[2]);
                double recvvx, recvvy, recvvz;
                if (m < 1e-5) {
//                    FatalErrorIn ( "CPLSocketFOAM::unpackVelocity()")
//                        << "Zero molecules in boundary region or averaging error"
//                        << exit(FatalError);
				    recvvx = 0.0;
        			recvvy = 0.0;
        			recvvz = 0.0;
                } else {
    				recvvx = recvVelocityBuff(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
        			recvvy = recvVelocityBuff(1, loc_cell[0], loc_cell[1], loc_cell[2])/m;
        			recvvz = recvVelocityBuff(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;
                }

                //Note here velocity is set straight to MD average value
                if (interp_BC == false) {
				    if (applyBCx) rvPatch[faceI].x() = recvvx;
				    if (applyBCy) rvPatch[faceI].y() = recvvy;
				    if (applyBCz) rvPatch[faceI].z() = recvvz;

//                    Foam::Info << "recvBuf with no interp " << interp_BC << " " 
//                                << facex << " "  << facey << " "  << facez << " " 
//                                << recvvx << " " << recvvy << " " << recvvz << " "
//                                << rvPatch[faceI].x() << " " << rvPatch[faceI].y() 
//                                << " " << rvPatch[faceI].z() << " " << Foam::endl;


		        //or use interpolation assuming specified cell is the one outside the domain
                } else {

                    Foam::point closestCellCentre((glob_cell[0]+0.5)*dx, 
                                                  (glob_cell[1]+0.5)*dy, 
                                                  (glob_cell[2]+0.5)*dz);
                    Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);
                    //Foam::label cell = mesh.findCell(closestCellCentre);
                    if (applyBCx) rvPatch[faceI].x() = (recvvx + U[cell].x()) / 2.0;
                    if (applyBCy) rvPatch[faceI].y() = (recvvy + U[cell].y()) / 2.0;
                    if (applyBCz) rvPatch[faceI].z() = (recvvz + U[cell].z()) / 2.0;

//                    Foam::Info << "recvBuf " << facex << " " << facey << " " << facez << " " << cell << " "
//                                << recvvx << " " << recvvy << " " << recvvz << " "
//                                << rvPatch[faceI].x() << " " << rvPatch[faceI].y() << " " << rvPatch[faceI].z() << " "
//                                << rvPatchP[faceI] << " " << Foam::endl;

                }

                //Note here velocity is set straight to MD average value
		        //which may be correct or you may need to use interpolation
                // with something like (recvvx + U[cell].x()) / 2.0; where
		        //cell = mesh.findCell (closestCellCentre);
//				if (applyBCx) rvPatch[faceI].x() = recvvx;
//				if (applyBCy) rvPatch[faceI].y() = recvvy;
//				if (applyBCz) rvPatch[faceI].z() = recvvz;
//                if (interp_BC) {
//                    Foam::point closestCellCentre((glob_cell[0]+0.5)*dx, 
//                                                  (glob_cell[1]+0.5)*dy, 
//                                                  (glob_cell[2]+0.5)*dz);
//                    Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);
//                    //Foam::label cell = mesh.findCell(closestCellCentre);
//                    if (applyBCx) rvPatch[faceI].x() = (recvvx + U[cell].x()) / 2.0;
//                    if (applyBCy) rvPatch[faceI].y() = (recvvy + U[cell].y()) / 2.0;
//                    if (applyBCz) rvPatch[faceI].z() = (recvvz + U[cell].z()) / 2.0;


//                    Foam::Info << "recvBuf " << facex << " " << facey << " " << facez << " " << cell << " "
//                                << recvvx << " " << recvvy << " " << recvvz << " "
//                                << rvPatch[faceI].x() << " " << rvPatch[faceI].y() << " " << rvPatch[faceI].z() << " "
//                                << Foam::endl;

//                }
            }

        }
    }
    return 0;
}


// Unpacks the 3 components of the velocity-tensor from the socket's
// recvVelocity (cpl::ndArray) storage into a boundary condition.
double CPLSocketFOAM::unpackVelocityPressure(volVectorField &U, volScalarField &p, fvMesh &mesh) 
{

    //bool interp_BC = false;

	if (CPL::is_proc_inside(velBCPortion.data())) {

		// TODO: Make this a utility general function that can be used on buffers
		if (CPL::get<int>("cpl_cfd_bc_slice")) {
            FatalErrorIn ( "CPLSocketFOAM::unpackVelocityPressure()")
                << " CPL_CFD_BC_SLICE not supported for velocity and pressure. Aborting."
                << exit(FatalError);
		} 

		// Apply BCs only in certain directions.
		int applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
		int applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
		int applyBCz = CPL::get<int> ("cpl_cfd_bc_z");

		// Patch receiving B.Cs
		Foam::string receivePatchName ("CPLReceiveMD");
		Foam::label rvPatchID = mesh.boundaryMesh().findPatchID(receivePatchName);

		if (rvPatchID == -1) {
			FatalErrorIn ( "CPLSocketFOAM::unpack()")
				<< " Could not find patch ID " << receivePatchName << ". "
				   " Aborting."
				<< exit(FatalError);
		}

		Foam::fvPatchVectorField& rvPatch = U.boundaryField()[rvPatchID];
		Foam::fvPatchScalarField& rvPatchP = p.boundaryField()[rvPatchID];
		const Foam::vectorField BoundaryfaceCntr = mesh.boundary()[rvPatchID].Cf();

		Foam::label cell;
		Foam::point closestCellCentre;
		for (int faceI = 0; faceI != BoundaryfaceCntr.size(); ++faceI) {
			double facex = BoundaryfaceCntr[faceI].x();
			double facey = BoundaryfaceCntr[faceI].y();
			double facez = BoundaryfaceCntr[faceI].z();
			// Find the cell indices for this position recvVelocity(:, ix, iy, iz)
        	int glob_cell[3]; int loc_cell[3];
			CPL::map_coord2cell(facex, facey, facez, glob_cell);
            //glob_cell[1] += 1; // Add one as boundary outside overlap by construction
	        bool valid_cell = CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

//            Foam::Info << "recvBuf " << facex << " " << facey << " " << facez << " " << cell << " "
//                        << glob_cell[0] << " " << glob_cell[1] << " " << glob_cell[2] << " "
//                        << loc_cell[0] << " " << loc_cell[1] << " " << loc_cell[2] << " " << valid_cell 
//                        << Foam::endl;

            if (valid_cell) {

                double m = recvVelocityBuff(3, loc_cell[0], loc_cell[1], loc_cell[2]); 
                double recvvx, recvvy, recvvz, recvP;
                if (m < 1e-5) {
//                    FatalErrorIn ( "CPLSocketFOAM::unpackVelocity()")
//                        << "Zero molecules in boundary region or averaging error"
//                        << exit(FatalError);
				    recvvx = 0.0;
        			recvvy = 0.0;
        			recvvz = 0.0;
                    recvP = 0.0;
                } else {
				    recvvx = recvVelocityBuff(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
				    recvvy = recvVelocityBuff(1, loc_cell[0], loc_cell[1], loc_cell[2])/m;
				    recvvz = recvVelocityBuff(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;
				    recvP  = recvVelocityBuff(4, loc_cell[0], loc_cell[1], loc_cell[2]);
                }

                //Note here velocity is set straight to MD average value
                if (interp_BC == false) {
				    if (applyBCx) rvPatch[faceI].x() = recvvx;
				    if (applyBCy) rvPatch[faceI].y() = recvvy;
				    if (applyBCz) rvPatch[faceI].z() = recvvz;
				    rvPatchP[faceI] = recvP;

//                    Foam::Info << "recvBuf with no interp " << interp_BC << " " << facex << " " << facey << " " << facez << " " << cell << " "
//                                << recvvx << " " << recvvy << " " << recvvz << " "
//                                << rvPatch[faceI].x() << " " << rvPatch[faceI].y() << " " << rvPatch[faceI].z() << " "
//                                << rvPatchP[faceI] << " " << Foam::endl;


		        //or use interpolation assuming specified cell is the one outside the domain
                } else {

                    Foam::point closestCellCentre((glob_cell[0]+0.5)*dx, 
                                                  (glob_cell[1]+0.5)*dy, 
                                                  (glob_cell[2]+0.5)*dz);
                    Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);
                    //Foam::label cell = mesh.findCell(closestCellCentre);
                    if (applyBCx) rvPatch[faceI].x() = (recvvx + U[cell].x()) / 2.0;
                    if (applyBCy) rvPatch[faceI].y() = (recvvy + U[cell].y()) / 2.0;
                    if (applyBCz) rvPatch[faceI].z() = (recvvz + U[cell].z()) / 2.0;
                    rvPatchP[faceI] = (recvP + p[cell]) / 2.0;

//                    Foam::Info << "recvBuf " << facex << " " << facey << " " << facez << " " << cell << " "
//                                << recvvx << " " << recvvy << " " << recvvz << " "
//                                << rvPatch[faceI].x() << " " << rvPatch[faceI].y() << " " << rvPatch[faceI].z() << " "
//                                << rvPatchP[faceI] << " " << Foam::endl;

                }
            }

        }
    }
}


// Unpacks the components from the socket's
double CPLSocketFOAM::unpackPorousForce(volVectorField &F, volScalarField &eps, fvMesh &mesh) 
{

    for (int ix=0; ix<recvBuf.shape(1); ix++) {
        for (int iy=0; iy<recvBuf.shape(2); iy++) {
            for (int iz=0; iz<recvBuf.shape(3); iz++) {

                double Fx = recvBuf(0, ix, iy, iz);
                double Fy = recvBuf(1, ix, iy, iz);
                double Fz = recvBuf(2, ix, iy, iz);
                double e  = recvBuf(3, ix, iy, iz);

                double glob_pos[3];
                CPL::map_cell2coord(ix, iy, iz, glob_pos);
                Foam::point closestCellCentre(glob_pos[0]+0.5*dx, glob_pos[1]+0.5*dy, glob_pos[2]+0.5*dz);
                Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);
                eps[cell] = e;
                F[cell].x() = Fx;
                F[cell].y() = Fy;
                F[cell].z() = Fz;

            }
        }
    }
}


// Unpacks the components from the socket's
double CPLSocketFOAM::unpackPorousVelForceCoeff(volVectorField &U, 
                                                volVectorField &F, 
                                                volScalarField &Fcoeff, 
                                                volScalarField &eps,
                                                scalar maxPossibleAlpha,
                                                fvMesh &mesh) 
{

    for (int ix=0; ix<recvBuf.shape(1); ix++) {
        for (int iy=0; iy<recvBuf.shape(2); iy++) {
            for (int iz=0; iz<recvBuf.shape(3); iz++) {

                //The recieved values are the raw sums collected by particle solver
                double Uxsum = recvBuf(0, ix, iy, iz);
                double Uysum = recvBuf(1, ix, iy, iz);
                double Uzsum = recvBuf(2, ix, iy, iz);
                double Fxsum = recvBuf(3, ix, iy, iz);
                double Fysum = recvBuf(4, ix, iy, iz);
                double Fzsum = recvBuf(5, ix, iy, iz);
                double Cdsum = recvBuf(6, ix, iy, iz);
                double volSum  = recvBuf(7, ix, iy, iz);

                //Get OpenFOAM cell corresponding to CPL grid
                double glob_pos[3];
                CPL::map_cell2coord(ix, iy, iz, glob_pos);
                Foam::point closestCellCentre(glob_pos[0]+0.5*dx, glob_pos[1]+0.5*dy, glob_pos[2]+0.5*dz);
                Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);
                double Vcell = mesh.V()[cell];

//                Foam::Info << "recvBuf " << ix << " " << iy << " " << iz << " " << cell << " "
//                            << Fxsum << " " << Fysum << " " << Fzsum << " "
//                            << Uxsum << " " << Uysum << " " << Uzsum << " "
//                            << Cdsum << " " << volSum << " " << maxPossibleAlpha << Foam::endl;

                //N.B we use phi/eps instead of alpha/beta instead here (in keeping with granular literature)
                double phi = volSum/Vcell;
                //Foam::Info << "CPLSocketFOAM::unpackPorousVelForceCoeff recvBuf " << ix << " " << iy << " " << iz << " " << cell << " "
                //            << phi << " " << volSum << " " << Vcell << " " << maxPossibleAlpha << Foam::endl;
                if (phi > maxPossibleAlpha) {
                    //Default value set in createFields or read from transportProperties
                    eps[cell] = 1.0 - maxPossibleAlpha;
                    phi = 1 - eps[cell];
                } else {
                    eps[cell] = 1.0 - phi;
                }

                // F and Fcoeff are multiplied by instantaneous eps when sent
                // from MD and are divided here by the time-average eps to
                // ensure momentum conservation.
                Fcoeff[cell] = Cdsum/(eps[cell]*Vcell);
                F[cell].x() = Fxsum/(eps[cell]*Vcell);
                F[cell].y() = Fysum/(eps[cell]*Vcell);
                F[cell].z() = Fzsum/(eps[cell]*Vcell);
                
                // Previous incorrect version of Ua
                // U[cell].x() = Uxsum/(eps[cell]*Vcell);
                // U[cell].y() = Uysum/(eps[cell]*Vcell);
                // U[cell].z() = Uzsum/(eps[cell]*Vcell);

                // Corrected Ua expression. For fluid only cells, set Ua = 0.
                // To avoid hard-wiring a limit, use the maxPossibleAlpha
                // value.
                if (phi < (1-maxPossibleAlpha))
                {
                    U[cell].x() = 0.0;
                    U[cell].y() = 0.0;
                    U[cell].z() = 0.0;
                }
                else
                {
                    U[cell].x() = Uxsum/(phi*Vcell);
                    U[cell].y() = Uysum/(phi*Vcell);
                    U[cell].z() = Uzsum/(phi*Vcell);
                }


#if DEBUG
                if (eps[cell] != 1){
                    Foam::Info << "CPLSocketFOAM::unpackPorousVelForceCoeff recv " << ix << " " << iy << " " << iz << " " << cell << " "
                            << F[cell].x() << " " << F[cell].y() << " " << F[cell].z() << " "
                            << U[cell].x() << " " << U[cell].y() << " " << U[cell].z() << " "
                            << Fcoeff[cell] << " " << eps[cell] << " " << maxPossibleAlpha << Foam::endl;
                }
#endif


//                double Ux = recvBuf(0, ix, iy, iz);
//                double Uy = recvBuf(1, ix, iy, iz);
//                double Uz = recvBuf(2, ix, iy, iz);
//                double Fx = recvBuf(3, ix, iy, iz);
//                double Fy = recvBuf(4, ix, iy, iz);
//                double Fz = recvBuf(5, ix, iy, iz);
//                double Cd = recvBuf(6, ix, iy, iz);
//                double e  = recvBuf(7, ix, iy, iz);

//                double glob_pos[3];
//                CPL::map_cell2coord(ix, iy, iz, glob_pos);
//                Foam::point closestCellCentre(glob_pos[0]+0.5*dx, glob_pos[1]+0.5*dy, glob_pos[2]+0.5*dz);
//                Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);
//                eps[cell] = e;
//                Fcoeff[cell] = Cd;
//                F[cell].x() = Fx;
//                F[cell].y() = Fy;
//                F[cell].z() = Fz;
//                U[cell].x() = Ux;
//                U[cell].y() = Uy;
//                U[cell].z() = Uz;

//                if (eps[cell] == 0.26){
//                    Foam::Info << "recvBuf " << ix << " " << iy << " " << iz << " " 
//                                << F[cell].x() << " " << F[cell].y() << " " << F[cell].z() << " "
//                                << U[cell].x() << " " << U[cell].y() << " " << U[cell].z() << " "
//                                << eps[cell] << " " << Fcoeff[cell] << Foam::endl;
 //               }


            }
        }
    }
}


volVectorField CPLSocketFOAM::divideFieldsVectorbyScalar(volVectorField &V, volScalarField &S, fvMesh &mesh) 
{

    for (int ix=0; ix<recvBuf.shape(1); ix++) {
        for (int iy=0; iy<recvBuf.shape(2); iy++) {
            for (int iz=0; iz<recvBuf.shape(3); iz++) {

                double glob_pos[3];
                CPL::map_cell2coord(ix, iy, iz, glob_pos);
                Foam::point closestCellCentre(glob_pos[0]+0.5*dx, glob_pos[1]+0.5*dy, glob_pos[2]+0.5*dz);
                Foam::label cell = meshSearcher->findNearestCell(closestCellCentre);

//                Foam::Info << "Divide F " << cell << " " << F[cell].x()
//                           << " " << F[cell].y() << " " << F[cell].z() 
//                           << " " << S[cell] << Foam::endl;

                if (S[cell] > 1e-6) {
                    V[cell].x() = V[cell].x()/S[cell];
                    V[cell].y() = V[cell].y()/S[cell];
                    V[cell].z() = V[cell].z()/S[cell];
                } else {
		            FatalErrorIn ( "CPLSocketFOAM::divideFieldsVectorbyScalar()")
		                     << "Divide by approx zero "<< cell << " " << S[cell] 
                            << " " << V[cell].x() << " " << V[cell].y() 
                            << " " << V[cell].z() << exit(FatalError);
                }

            }
        }
    }

    return V;
}


// Sends buffer to overlapping MD processes.
void CPLSocketFOAM::send() 
{
    CPL::send(sendBuf.data(), sendBuf.shapeData(), cnstFRegion.data());
}

// Receives buffer from overlapping MD processes.
void CPLSocketFOAM::recv()
{
    // LAMMPS olap size field
    int recvShape[4] = {8, olapCells[0], olapCells[1], olapCells[2]};
    recvBuf.resize(4, recvShape);

    CPL::recv(recvBuf.data(), recvBuf.shapeData(), olapPortion.data());

//    for (int ix=0; ix<recvBuf.shape(1); ix++) {
//    for (int iy=0; iy<recvBuf.shape(2); iy++) {
//    for (int iz=0; iz<recvBuf.shape(3); iz++) {
//    for (int n = 0; n < 8; n++) {
//        std::cout << "RECV " << ix << " " << iy << " " << iz << " " << n << " " 
//                  << recvBuf(n, ix, iy, iz) << std::endl;
//    }}}}

}

// Sends 9 components of the stress-tensor to overlapping MD processes.
void CPLSocketFOAM::sendStress() 
{
    CPL::send (sendStressBuff.data(), sendStressBuff.shapeData(), cnstFRegion.data());
}

// Receives 3 components of the velocity vector & 1 no. of particles from overlapping MD processes.
void CPLSocketFOAM::recvVelocity()
{
    // LAMMPS computed fields
    int recvVelocityShape[4] = {4, velBCCells[0], velBCCells[1], velBCCells[2]};
    recvVelocityBuff.resize(4, recvVelocityShape);
    CPL::recv(recvVelocityBuff.data(), recvVelocityBuff.shapeData(), 
              velBCPortion.data());

//    for (int ix=0; ix<recvVelocityBuff.shape(1); ix++) {
//    for (int iy=0; iy<recvVelocityBuff.shape(2); iy++) {
//    for (int iz=0; iz<recvVelocityBuff.shape(3); iz++) {
//    for (int n = 0; n < 4; n++) {
//        std::cout << "RECV i j k n buf " << ix << " " << iy << " " << iz << " " << n << " " 
//                  << recvVelocityBuff(n, ix, iy, iz) << std::endl;
//    }}}}
}
  
// Receives 3 components of the velocity vector, 1 no. partciles and 1 pressure from overlapping MD processes.
void CPLSocketFOAM::recvVelocityPressure()
{
    // LAMMPS computed fields
    int recvVelocityShape[4] = {5, velBCCells[0], velBCCells[1], velBCCells[2]};
    recvVelocityBuff.resize(4, recvVelocityShape);
    CPL::recv(recvVelocityBuff.data(), recvVelocityBuff.shapeData(), 
              velBCPortion.data());
//    for (int ix=0; ix<recvVelocityBuff.shape(1); ix++) {
//    for (int iy=0; iy<recvVelocityBuff.shape(2); iy++) {
//    for (int iz=0; iz<recvVelocityBuff.shape(3); iz++) {
//    for (int n = 0; n < 5; n++) {
//        std::cout << "RECV i j k n buf " << ix << " " << iy << " " << iz << " " << n << " " 
//                  << recvVelocityBuff(n, ix, iy, iz) << std::endl;
//    }}}}
}

