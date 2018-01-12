#include "IcoFOAMOutgoingField.H"


void StressOutgoingField::pack_(const std::vector<int>& glob_cell,
                                const std::vector<int>& loc_cell,
                                const std::valarray<double>& coord) {

    //TODO: Check for mode=0 to work letting one more cell in overalp region
    Foam::label cell;
    // Default compute_mode=="cell"
    double ycenter_dy = 0.5;
	if (compute_mode == "plane")
		ycenter_dy = 1.0;
    Foam::point globalPos = Foam::point((glob_cell[0] + 0.5) * dx,
									    (glob_cell[1] + ycenter_dy) * dy, 
									    (glob_cell[2] + 0.5) * dz);

	cell = meshSearcher->findNearestCell(globalPos);
    if (cell != -1) {
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
    else
        std::cerr << "Warning: The point (" << (glob_cell[0]+0.5)*dx << "," << 
                    (glob_cell[1]+0.5)*dy << "," << (glob_cell[2]+0.5)*dz << 
                    ") is outside the mesh for whatever reason! - Cell: " << 
                    cell << std::endl;

}

void StressOutgoingField::setup() {
    CPL::get_file_param("constrain.momentum", "compute-mode", compute_mode);
    if (compute_mode != "plane" && compute_mode != "cell")
        std::cout<< "Error mode compute stress." << "mode " << compute_mode <<std::endl;
    meshSearcher = new Foam::meshSearch(*foamMesh);
    data_size = 9;
    update();
}

void StressOutgoingField::update() {
    Foam::dimensionedScalar mu(foamDensity*kViscosity);
    //TODO: Should be mu not nu : Tests with density equals 1.0
	Foam::volSymmTensorField sigma_vol(kViscosity*2*dev(symm(fvc::grad(*velocityField))));
	if (compute_mode == "plane")
		stressField = (Foam::symmTensorField) fvc::interpolate(sigma_vol);
    else if (compute_mode == "cell")
		stressField = (Foam::symmTensorField) sigma_vol;
}
 
