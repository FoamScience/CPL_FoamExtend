#include "VelocityIncomingField.H"

//TODO: Create a function to average DataFields
/**
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
        }
    }
}
**/

bool VelIncomingField::unpackCustom_() {
    std::valarray<double> face_center(3);
    std::vector<int> glob_cell(3), loc_cell(3);
    bool valid_cell;
    double bcvx, bcvy, bcvz;
    forAll(bcFaceCenters, faceI) {
        face_center[0] = bcFaceCenters[faceI].x(); 
        face_center[1] = bcFaceCenters[faceI].y(); 
        face_center[2] = bcFaceCenters[faceI].z(); 
        CPL::map_coord2cell(face_center[0], face_center[1], face_center[2], glob_cell.data());
        loc_cell = getLocalCell(glob_cell, valid_cell);
        if (valid_cell) {
            bcvx = buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]);
            bcvy = buffer(1, loc_cell[0], loc_cell[1], loc_cell[2]);
            bcvz = buffer(2, loc_cell[0], loc_cell[1], loc_cell[2]);
            if (applyBc[0]) (*bcPatch)[faceI].x() = bcvx;
            if (applyBc[1]) (*bcPatch)[faceI].y() = bcvy;
            if (applyBc[2]) (*bcPatch)[faceI].z() = bcvz;
        }
    }
    return true;
}

// void VelIncomingField::unpack_(const std::vector<int>& glob_cell,
//                                const std::vector<int>& loc_cell,
//                                const std::valarray<double>& coord) {
//     std::valarray<double> tol({1e-6, 1e-6, 1e-6});
//     std::valarray<double> face_center(3);
//     std::valarray<double> coord_center = coord;
//     coord_center[0] += dx/2.0;
//     coord_center[2] += dz/2.0;
//     double bcvx, bcvy, bcvz;
//     forAll(bcFaceCenters, faceI) {
//         face_center[0] = bcFaceCenters[faceI].x(); 
//         face_center[1] = bcFaceCenters[faceI].y(); 
//         face_center[2] = bcFaceCenters[faceI].z(); 
//         if ((std::abs(face_center - coord_center) < tol).min()) {
//             bcvx = buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]);
//             bcvy = buffer(1, loc_cell[0], loc_cell[1], loc_cell[2]);
//             bcvz = buffer(2, loc_cell[0], loc_cell[1], loc_cell[2]);
//             if (applyBc[0]) (*bcPatch)[faceI].x() = bcvx;
//             if (applyBc[1]) (*bcPatch)[faceI].y() = bcvy;
//             if (applyBc[2]) (*bcPatch)[faceI].z() = bcvz;
//         }
//
//     }
// }
void VelIncomingField::setup() {
        data_size = 4;
        // Apply BCs only in certain directions.
        applyBc = std::vector<int>(3);
        CPL::get_file_param("bc.velocity", "components", applyBc);

        // Patch receiving B.Cs
        Foam::string bcPatchName ("CPLReceiveMD");
        Foam::label bcPatchID = foamMesh->boundary().findPatchID(bcPatchName);

        if (bcPatchID == -1) {
            FatalErrorIn ("CPLSocketFOAM::unpack()")
                << " Could not find patch ID " << bcPatchName << ". "
                   " Aborting."
                << exit(FatalError);
        }
        bcPatch = &velocityField->boundaryField()[bcPatchID];
        bcFaceCenters = foamMesh->boundary()[bcPatchID].Cf();

}
void VelIncomingField::update() {
 
}
