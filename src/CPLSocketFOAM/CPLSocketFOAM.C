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

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

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
#include "CPLSocketFOAM.H"
#include "blockMesh.H"
#include <sstream>
#include <unistd.h>
#include "interpolation.H"
#include "PstreamGlobals.H" 


// Initialise CFD realm communicator
void CPLSocketFOAM::initComms (int& argc, char**& argv) {
	MPI_Init(&argc, &argv);
    CPLSocket::initComms();
	Foam::PstreamGlobals::CPLRealmComm = realmComm;

}


void CPLSocketFOAM::setTimingInfo() {
    dt = foamRunTime->deltaTValue();
    nSteps = nint((foamRunTime->endTime().value() - \
                   foamRunTime->startTime().value()) / dt);
}


void CPLSocketFOAM::setCartCommInfo() {
    Foam::IOdictionary decomposeDict(Foam::IOobject ("decomposeParDict", 
									 foamRunTime->time().system(), *foamRunTime,
                        			 IOobject::MUST_READ, IOobject::NO_WRITE,
									 false));
    Foam::dictionary simpleCoeffs = decomposeDict.subDict ("simpleCoeffs");
    Foam::Vector<int> np = simpleCoeffs.lookup ("n");
    procGrid = std::vector<int>({np.x(), np.y(), np.z()});

    // Create custom cartesian communicator (cartComm) based on myCoords	
	// Assume xyz ordering (Simple decomposition)
	myProcCoords[2] = Pstream::myProcNo() / (np.x()*np.y());
	int mod_aux = Pstream::myProcNo() % (np.x()*np.y());
	myProcCoords[1] = mod_aux / np.x();
	myProcCoords[0] = mod_aux % np.x();
}

void CPLSocketFOAM::setRealmDomainInfo() {
    Foam::IOdictionary blockMeshDict(Foam::IOobject ("polyMesh/blockMeshDict", 
                                     foamRunTime->time().constant(), *foamRunTime,
                        			 IOobject::MUST_READ, 
									 IOobject::NO_WRITE, false));
	Foam::List<Foam::Vector<double>> vertices(blockMeshDict.lookup("vertices"));

    // Domain dimensions
    std::valarray<double> domain_length({vertices[1][0] - vertices[0][0],
                                         vertices[3][1] - vertices[0][1],
                                         vertices[4][2] - vertices[0][2]});
    Foam::word dummyRegionName("dummy");
    Foam::blockMesh blocks(blockMeshDict, dummyRegionName);
    Foam::Vector<int> meshDensity = blocks[0].meshDensity();
    // Global number of cells
    int ncx, ncy, ncz;
    CPL::get_file_param("coupling-mesh", "ncx", ncx);
    CPL::get_file_param("coupling-mesh", "ncy", ncy);
    CPL::get_file_param("coupling-mesh", "ncz", ncz);
    cfdCells = std::vector<int>({ncx, ncy, ncz});
    // cfdCells = std::vector<int>({meshDensity.x(), meshDensity.y(),
    //                              meshDensity.z()});
    std::valarray<double> domain_orig({0.0, 0.0, 0.0});
    realmDomain = CPL::Domain(domain_orig, domain_length);
}

void CPLSocketFOAM::setOpenFOAM(const Foam::Time &run_time, const Foam::fvMesh &mesh) {
    foamMesh = &mesh;
    foamRunTime = &run_time;
}


void CPLSocketFOAM::init() {
    CPLSocket::init();
}

// double CPLSocketFOAM::\
// unpackVelocity(volVectorField &U, fvMesh &mesh) {
//
//
// 	if (CPL::is_proc_inside(bcPortion.data())) { 
//
// 		// TODO: Make this a utility general function that can be used on buffers
// 		if (CPL::get<int>("cpl_cfd_bc_slice")) {
//
// 			Foam::Info << "CPL_CFD_BC_SLICE is on: averaging CFD recvVelocity "
// 						  "in the x-z plane" << Foam::endl;
//
// 			// Number of cells in the local processor in x-z plane
// 			int N = recvVelocityBuff.shape(1) * recvVelocityBuff.shape(3);
//
// 			// For every component and y-value 
// 			for (int j = 0; j < recvVelocityBuff.shape(2); ++j) {
// 				for (int c = 0; c < recvVelocityBuff.shape(0); ++c) {
// 				// Sum across the x-z plane 
// 				double total = 0.0;
// 				for (int k = 0; k < recvVelocityBuff.shape(3); ++k)
// 					for (int i = 0; i < recvVelocityBuff.shape(1); ++i)
// 						total += recvVelocityBuff(c, i, j, k);
//
// 			// Find mean by dividing sum by number of cells 
// 				for (int k = 0; k < recvVelocityBuff.shape(3); ++k)
// 					for (int i = 0; i < recvVelocityBuff.shape(1); ++i)
// 						recvVelocityBuff(c, i, j, k) = total / static_cast<double> (N);
// 				}
// 			}
// 		} 
//         
//         //TODO:SETUP
// 		// Apply BCs only in certain directions.
// 		int applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
// 		int applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
// 		int applyBCz = CPL::get<int> ("cpl_cfd_bc_z");
//
// 		// Patch receiving B.Cs
// 		Foam::string receivePatchName ("CPLReceiveMD");
// 		Foam::label rvPatchID = mesh.boundary().findPatchID(receivePatchName);
//
// 		if (rvPatchID == -1) {
// 			FatalErrorIn ( "CPLSocketFOAM::unpack()")
// 				<< " Could not find patch ID " << receivePatchName << ". "
// 				   " Aborting."
// 				<< exit(FatalError);
// 		}
//
// 		Foam::fvPatchVectorField& rvPatch = U.boundaryField()[rvPatchID];
// 		const Foam::vectorField faceCenters = mesh.boundary()[rvPatchID].Cf();
//         //TODO:SETUP
//
//
//         std::valarray<double> tol({1e-6, 1e-6, 1e-6});
//         std::valarray<double> face_center(3);
//         std::vector<int> glob_cell(3), loc_cell(3);
//         std::valarray<double> coord(3);
// 		for (int i = bcPortion[0]; i <= bcPortion[1]; i++)
//         for (int j = bcPortion[2]; j <= bcPortion[3]; j++)
//         for (int k = bcPortion[4]; k <= bcPortion[5]; k++) {
//             glob_cell = std::vector<int>({i, j, k});
// 			CPL::map_glob2loc_cell(bcPortion.data(), glob_cell.data(), 
//                                    loc_cell.data());
//             CPL::map_cell2coord(glob_cell[0], glob_cell[1], glob_cell[2], &coord[0]); 
//             //Center the face
//             coord[0] += dx/2.0;
//             coord[2] += dz/2.0;
//             forAll(faceCenters, faceI) {
//                 face_center[0] = faceCenters[faceI].x(); 
//                 face_center[1] = faceCenters[faceI].y(); 
//                 face_center[2] = faceCenters[faceI].z(); 
//                 if ((std::abs(face_center - coord) < tol).min()) {
//                     double m = 1.0;
//                     double recvvx = recvVelocityBuff(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
//                     double recvvy = recvVelocityBuff(1, loc_cell[0], loc_cell[1], loc_cell[2]) /m;
//                     double recvvz = recvVelocityBuff(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;
//         			if (applyBCx) rvPatch[faceI].x() = recvvx;
// 	        		if (applyBCy) rvPatch[faceI].y() = recvvy;
// 		        	if (applyBCz) rvPatch[faceI].z() = recvvz;
//                 }
//                     
//             }
//         }
//
// 		if (applyBCx) Foam::Info << "MD->CFD BC x-velocity applied." << Foam::endl;
// 		if (applyBCy) Foam::Info << "MD->CFD BC y-velocity applied." << Foam::endl;
// 		if (applyBCz) Foam::Info << "MD->CFD BC z-velocity applied." << Foam::endl;
// 	}
//
// }
