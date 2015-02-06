/*
 *  os3_picker_starter.h
 *  
 *
 *  Created by Thomas Hrabe on 09.03.08.
 *
 */


#include "tom/mpi/mpi_fcn.hpp"
#include <string>

namespace tom{
	void os3_printHelp();
	void os3_server(MPI_Comm comm,std::string optionPath);
	void os3_client(MPI_Comm comm,int rank);
}



