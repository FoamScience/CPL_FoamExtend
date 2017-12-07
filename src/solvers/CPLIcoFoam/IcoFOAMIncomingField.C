#include "IcoFOAMIncomingField.H"
#include "cpl/cpl.h"

void VelIncomingField::unpack_(const std::vector<int>& glob_cell,
                               const std::vector<int>& loc_cell,
                               const std::valarray<double>& coord) {
    std::valarray<double> tol({1e-6, 1e-6, 1e-6});
    std::valarray<double> face_center(3);
    std::valarray<double> coord_center = coord;
    coord_center[0] += dx/2.0;
    coord_center[2] += dz/2.0;
    forAll(bcFaceCenters, faceI) {
        face_center[0] = bcFaceCenters[faceI].x(); 
        face_center[1] = bcFaceCenters[faceI].y(); 
        face_center[2] = bcFaceCenters[faceI].z(); 
        if ((std::abs(face_center - coord_center) < tol).min()) {
            double m = 1.0;
            double bcvx = buffer(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
            double bcvy = buffer(1, loc_cell[0], loc_cell[1], loc_cell[2]) /m;
            double bcvz = buffer(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;
            if (applyBCx) (*bcPatch)[faceI].x() = bcvx;
            if (applyBCy) (*bcPatch)[faceI].y() = bcvy;
            if (applyBCz) (*bcPatch)[faceI].z() = bcvz;
        }

    }
}
void VelIncomingField::setup() {
        data_size = 4;
        // Apply BCs only in certain directions.
        applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
        applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
        applyBCz = CPL::get<int> ("cpl_cfd_bc_z");

        // Patch receiving B.Cs
        Foam::string bcPatchName ("CPLReceiveMD");
        Foam::label bcPatchID = foamMesh->boundary().findPatchID(bcPatchName);

        if (bcPatchID == -1) {
            FatalErrorIn ( "CPLSocketFOAM::unpack()")
                << " Could not find patch ID " << bcPatchName << ". "
                   " Aborting."
                << exit(FatalError);
        }
        bcPatch = &velocityField->boundaryField()[bcPatchID];
        bcFaceCenters = foamMesh->boundary()[bcPatchID].Cf();

}
void VelIncomingField::update() {
 
}
