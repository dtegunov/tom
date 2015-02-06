/*
 *  os3_picker_starter.cpp
 *  
 *
 *  Created by Thomas Hrabe on 09.03.08.
 */

#include "tom/os3/os3_picker_starter.hpp"
#include "tom/os3/os3_options.hpp"
#include "tom/os3/os3_picker.hpp"
#include "tom/os3/os3_functions.hpp"
#include "tom/mpi/mpi_fcn.hpp"
#include "tom/core/fftw_plan.hpp"

#include "helper/filesystem.hpp"

#include <iostream>
#include <vector>



namespace tom{
/****************************************************************************//**
 * \brief Prints help.
 *
 * Prints help.
 *******************************************************************************/
	void os3_printHelp(){
		std::cout << "OSCAR 3 start options: "<< std::endl;
		
		std::cout << " -h Print this help" << std::endl;
		std::cout << " -cf --create_folders <path> Creates folders structure: " << std::endl;
		std::cout << " \t volumes - folder with images / volumes" << std::endl;
		std::cout << " \t patterns - folder with the patterns searched" << std::endl;
		std::cout << " \t masks - masks for the patterns " << std::endl;
		std::cout << " \t psf - point spread function for each pattern " << std::endl;
		std::cout << " \t results - folder storing each result " << std::endl;  
		std::cout << " \t options - folder with an options file templage " << std::endl;
		std::cout <<  std::endl;
		std::cout <<  std::endl;
		std::cout << " Note: if you want to start the picker with lam use \n\t mpirun -np <number of nodes> os3par <options file>" << std::endl;
	}









/****************************************************************************//**
 * \brief Job server for parallel processing.
 *
 * Job server for parallel processing.
 *******************************************************************************/
	void os3_server(MPI_Comm comm,std::string optionsPath){
		tom::os3_options opt;

		opt.readFromFile(optionsPath);
		std::vector<tom::os3_job> jobList = opt.createJobList();
		std::cout << "Number of jobs =  " << jobList.size() << std::endl;
		std::vector<tom::os3_job>::iterator iterator;
		
		std::vector<tom::os3_job> finishedJobList;
			
		int clientID = MPI_ANY_SOURCE;
		iterator = jobList.begin();

		while(iterator < jobList.end()){
	
			tom::mpi::os3_recieve_job_request(comm,clientID);
			std::cout << "* Client : " << clientID << " starts job " << (*iterator).id << std::endl; 
			tom::mpi::send_os3_job(clientID, 1, comm, *(iterator));
	
			tom::os3_job* newJob = tom::checkForFinishedVolumes(finishedJobList);
	
			if(newJob->jobType == __TOM_OS3_PICKLIST_JOB__){
			//insert the stack job into the job queue
			//set the correlation filter settings of the stack job
			//setJobFilterWeight(newJob,opt.getFilterWeight());
			
				jobList.insert(iterator+1,*newJob);
			}

			if((*iterator++).jobType == __TOM_OS3_CONTAINER_JOB__){
				std::cout << "WAIT\n";
				tom::mpi::os3_recieve_job_request(comm,clientID);
				std::cout << "* Client : " << clientID << " starts job " << (*iterator).id << std::endl; 
				tom::mpi::send_os3_job(clientID, 1, comm, *(iterator++));
			}
			clientID= MPI_ANY_SOURCE;

			


		}
	
		int shutOutClients = 0;
		int numberClients;
		MPI_Comm_size(MPI_COMM_WORLD, &numberClients);
		tom::os3_job j;
		j.jobType = __TOM_OS3_NO_JOB__;
		while(shutOutClients< numberClients-1){
			clientID = MPI_ANY_SOURCE;
			tom::mpi::os3_recieve_job_request(comm,clientID);
			tom::mpi::send_os3_job(clientID, 1, comm, j);
			std::cout << "* Client : " << clientID << " is shutting down. " << std::endl; 
			shutOutClients++;
		}
		std::cout << "* " << "Server is shutting down. " << std::endl;


	}


/****************************************************************************//**
 * \brief Client for parallel processing.
 *
 * Client for parallel processing.
 *******************************************************************************/

	void os3_client(MPI_Comm comm,int rank){
	
		bool finish = false;
		bool display = false;
		tom::os3_job job;
	
		while(!finish){
			tom::mpi::os3_send_job_request(0,0,comm,rank);
			tom::mpi::recv_os3_job(0, 1,comm,job);
			finish = job.jobType == __TOM_OS3_NO_JOB__;
			tom::os3_picker<TFLOAT>(job,display);
		}
	
	}
}


int main(int argc,char** argv){
// start parallel picker

std::cout << "OSCAR 3 v.03 - 08.03.2008" << std::endl;

if(argc == 2 && std::string(argv[1]).find(std::string("-h")) != std::string::npos){
	tom::os3_printHelp();
	return 0;
}

if(argc == 3 && (std::string(argv[1]).find(std::string("-cf")) != std::string::npos || std::string(argv[1]).find(std::string("--create_folders")) != std::string::npos)){
	
	tom::os3_options::createDirectories(std::string(argv[2]));
	return 0;
}

if(argc > 3){
	std::cerr << "Specify a valid options file when starting OSCAR 3 or start with the following options" << std::endl;
	tom::os3_printHelp();
	return 0;
}

	
std::string argument(argv[1]);
	
if(argument.find(std::string("-")) != std::string::npos){
	std::cerr << "Please specify an option. Choose from the following options"<< std::endl;
	tom::os3_printHelp();
	return 0;
}

if(! helper::fs::is_regular(argument)){
	std::cerr << "The specified file does not exist!" << std::endl;
	std::cerr << "Either change the path to the options file or choose from the following options: " << std::endl;
	tom::os3_printHelp();
	return 0;
}

MPI_Init(&argc,&argv);
const int root = 0;

tom::fftw::setFcnMemFFTW<TFLOAT>();

int numprocs_world, my_rank_world;

MPI_Comm_size(MPI_COMM_WORLD, &numprocs_world);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_world);

MPI_Comm comm;
MPI_Comm_dup(MPI_COMM_WORLD, &comm);


if (my_rank_world == root){
	tom::os3_server(comm,argument);
}else{
	tom::os3_client(comm,my_rank_world);
}

MPI_Comm_free(&comm);

//tom::fftw::save_wisdom();

MPI_Finalize();

return 0;	
}
