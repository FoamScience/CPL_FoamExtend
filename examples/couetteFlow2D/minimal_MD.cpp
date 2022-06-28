// This is a minimal MD "script" to exchange U with CPLTestFoam

#include "mpi.h"
#include <iostream>

#include "cpl.h"
#include "CPL_ndArray.h"

using namespace std;

int main() {

   // MPI initialization
   MPI_Init(NULL, NULL); 
   int MD_realm = 2;
   MPI_Comm MD_COMM;
   CPL::init(MD_realm, MD_COMM);
   int nprocs_realm;
   MPI_Comm_size(MD_COMM, &nprocs_realm);


   // Parameters of the CPU topology (cartesian grid)
   // To work with the sample OpenFOAM case; in a 1-to-1 CPU setting
   int npxyz[3] = {1, 1, 1};
   int periods[3] = {1, 1, 1};
   int ncxyz[3] = {8, 8, 8};
   double xyzL[3] = {ncxyz[0]*2.099495, ncxyz[0]*5.668637, ncxyz[0]*2.099495};
   double xyz_orig[3] = {0.0, 0.0, 0.0};

   // Check that number of processors is consistent
   if (nprocs_realm != (npxyz[0] * npxyz[1] * npxyz[2])) {
      cout << "Non-coherent number of processes." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // Setup cartesian topology
   int rank;
   MPI_Comm_rank(MD_COMM, &rank);
   MPI_Comm CART_COMM;
   MPI_Cart_create(MD_COMM, 3, npxyz, periods, true, &CART_COMM);

   // Coupler setup
   CPL::setup_md(CART_COMM, xyzL, xyz_orig);

   // Get grid info
   int Ncells[3];
   int olap_limits[6], portion[6];
   CPL::get_olap_limits(olap_limits);
   CPL::my_proc_portion(olap_limits, portion);
   CPL::get_no_cells(portion, Ncells);

   // Receive the packed velocity from OpenFOAM
   int recv_shape[4] = {3, Ncells[0], Ncells[1], Ncells[2]};
   CPL::ndArray<double> recv_array(4, recv_shape);
   CPL::recv(recv_array.data(), recv_array.shapeData(), olap_limits);

   // TODO: check velocity values

   // Send artifical velocity to OpenFOAM
   // This tests  with 'valid cell indices'
   int send_shape[4] = {4, Ncells[0], Ncells[1], Ncells[2]};
   CPL::ndArray<double> send_array(4, send_shape);
   short ncells=0;
   for(int i=0; i<Ncells[0]; i++)
   for(int j=0; j<Ncells[1]; j++)
   for(int k=0; k<Ncells[2]; k++)
   {
   	    send_array(0,i,j,k) = 1;
   	    send_array(1,i,j,k) = 2;
   	    send_array(2,i,j,k) = 3;
   	    send_array(3,i,j,k) += ncells;
	    ncells++;
   }
   CPL::send(send_array.data(), send_array.shapeData(), olap_limits);

   // Release all MD comms 
   CPL::finalize();
   MPI_Comm_free(&MD_COMM);
   MPI_Comm_free(&CART_COMM);
   MPI_Finalize();
}
