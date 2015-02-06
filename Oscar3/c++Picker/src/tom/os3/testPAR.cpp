#include <iostream>
#include <fstream>
#include <assert.h>
#include <complex>
#include <sstream>

//-------------------------
#include "tom/os3/os3_volume.hpp"
#include "tom/os3/os3_functions.hpp"
#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_results.hpp"
#include "tom/os3/os3_options.hpp"
#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_picker.hpp"
//-------------------------
#include "tom/core/volume.hpp"
#include "tom/core/fftw_plan.hpp"
#include "tom/core/volume_fcn.hpp"	
#include "tom/mpi/mpi_fcn.hpp"




int main(int argc,char** argv){

	if(*argc < 2)
		std::cout << "spcify options file!" << std::endl;

	MPI_Init(&argc,&argv);

	const int root = 0;

	tom::fftw::setFcnMemFFTW<TFLOAT>();

	int numprocs_world, my_rank_world;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs_world);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_world);


	//tom::mpi::startup(vargv, argc, argv);

	MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

	
	
	if (my_rank_world == root){
		os3_server(comm,std::string(argv[1]));
	}else{
		os3_client(comm,my_rank_world);
	}
 
    MPI_Comm_free(&comm);

    //tom::fftw::save_wisdom();

    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

	return 0;
}





