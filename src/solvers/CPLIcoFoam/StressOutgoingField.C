#include "StressOutgoingField.H"
#include <sstream> 


void StressOutgoingField::pack_(const std::vector<int>& glob_cell,
                                const std::vector<int>& loc_cell,
                                const std::valarray<double>& coord) {

    Foam::label cell;
    Foam::point globalPos = Foam::point((glob_cell[0] + 0.5) * dx,
									    (glob_cell[1] + 0.5) * dy, 
									    (glob_cell[2] + 0.5) * dz);

    if (cell_array_init < nof_cells) {
        cell_array_init += 1;
        cell = meshSearcher->findNearestCell(globalPos);
        // std::stringstream ss;
        if (cell != -1) {
            // Do interpolation 
            if (compute_mode == "plane") {
                const Foam::cell& faces = foamMesh->cells()[cell];
                forAll( faces, i )  {
                    if (foamMesh->isInternalFace(faces[i])) {
                        Foam::vector faceINormal = foamMesh->Sf()[faces[i]] / foamMesh->magSf()[faces[i]] ; 
                        Foam::label cell_neigh = foamMesh->neighbour()[faces[i]];
                        // OpenFOAM set faces'  normal sign using cell numbering. From lower to higher is positive.
                        if (faceINormal == Foam::vector(0, 1, 0) && cell_neigh > cell) {
                            cell = faces[i];
                            // ss << " i = " << i << "- stress ("<< cell << "," << cell_neigh <<") :" << stressField[cell].xy() <<" - "<< ((*sigma)[cell2].xy() + (*sigma)[cell1].xy())/2.0;
                            // std::cout << ss.str() << std::endl;
                            break;
                        }
                    }
                }
            }
            cell_array(loc_cell[0], loc_cell[1], loc_cell[2]) = cell;
        } 
        else
            std::cerr << "Warning: The point (" << (glob_cell[0]+0.5)*dx << "," << 
                        (glob_cell[1]+0.5)*dy << "," << (glob_cell[2]+0.5)*dz << 
                        ") is outside the mesh for whatever reason! - Cell: " << 
                        cell << std::endl;
    }
    else {
        cell = cell_array(loc_cell[0], loc_cell[1], loc_cell[2]);
    }
    if (cell != -1) {
        // Update buffer
        buffer(0,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].xx();
        buffer(1,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].xy();
        buffer(2,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].xz();
        buffer(3,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].xy();
        buffer(4,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].yy();
        buffer(5,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].yz();
        buffer(6,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].xz();
        buffer(7,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].yz();
        buffer(8,loc_cell[0],loc_cell[1],loc_cell[2]) = stressField[cell].zz();
    }
}

void StressOutgoingField::setup() {
    CPL::get_file_param("constrain.momentum", "compute-mode", compute_mode);
    if (compute_mode != "plane" && compute_mode != "cell")
        std::cout<< "Error mode compute stress." << "mode " << compute_mode <<std::endl;
    cell_array =  CPL::ndArray<Foam::label>(3, nCells.data());
    cell_array_init = 0;
    nof_cells = nCells[0]*nCells[1]*nCells[2];
    meshSearcher = new Foam::meshSearch(*foamMesh);
    data_size = 9;
    update();
}

void StressOutgoingField::update() {
	if (compute_mode == "plane")
		stressField = (Foam::symmTensorField) fvc::interpolate(*sigma);
    else if (compute_mode == "cell")
		stressField = (Foam::symmTensorField) *sigma;
}
 
