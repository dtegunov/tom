
#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_volume.hpp"
#include "tom/os3/os3_functions.hpp"
#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_results.hpp"
#include "tom/os3/os3_picklist.hpp"
/*----------------------------------*/
#include "tom/core/volume.hpp"
#include "tom/core/volume_fcn.hpp"
/*----------------------------------*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


namespace tom{


/****************************************************************************//**
 * \brief The correlation based picker.
 *
 * \param[in] job The job structure. Specifies volume, pattern, result directory, ...
 * \param[in] printStatus If true, debug output will appear on the screen.
 * 
 * The correlation based picker. Applies the correlation filter to the volume specified in job.  
 *******************************************************************************/
template<typename T> void os3_correlationPicker(tom::os3_job &job,bool printStatus){
	

		tom::os3_volume<T> image(job.volumePath,job.subregion,NULL,job.binning);
		tom::os3_volume<T> patternOriginal(job.patternPath,NULL,NULL,job.binning);
		tom::os3_volume<T> mask(job.maskPath,NULL,NULL,job.binning);
		
		tom::os3_results<T> res;
	
		patternOriginal.calculateStatisticVolumes(&mask);
		//pattern will store the result of the flcf!!!
		tom::os3_volume<T> pattern(&patternOriginal);
		tom::os3_volume<T> soc(&patternOriginal);
		tom::os3_volume<T> psrV(&patternOriginal);
		double angles[3] ={0,0,0};
		
		//std::cout << "(" << patternOriginal.getSizeX() << "," << patternOriginal.getSizeY() << ","<< patternOriginal.getSizeZ() << ")\n";
		//std::cout << "(" << pattern.getSizeX() << "," << pattern.getSizeY() << ","<< pattern.getSizeZ() << ")\n";
	
		std::vector<tom::angleTupel> angleList = job.angleList;
		std::vector<tom::angleTupel>::iterator iterator;
		int iterationCounter = 0;
		
		if(printStatus){
			std::cout << "START" << std::endl;
			std::cout << "Number rotations: " << angleList.size() << std::endl;

		}
		for( iterator = angleList.begin(); iterator != angleList.end(); iterator++ ) {
		
			if(printStatus)
				std::cout << "Rotate" << std::endl;
			tom::angleTupelTOarray(*iterator,&angles[0]);
			//rotate the copied 
			pattern.rotate(angles);
			
			//std::cout << "(" << angles[0] << "," << angles[1] << "," << angles[2] << ")" << std::endl;
			//std::cout << std::endl << __FILE__ << " " << __LINE__ << std::endl;
			//copy the rotated pattern
			soc.setVolume(&pattern);
			//calculate the autocorrelation of the pattern - store it in soc
			if(printStatus)
				std::cout << "SOC1" << std::endl;
			tom::os3_flcf(&pattern,&soc,&mask);
			if(printStatus)
				std::cout << "FLCF" << std::endl;
			//classical flcf - store it in pattern
			tom::os3_flcf(&image,&pattern,&mask);
			if(printStatus)
				std::cout << "SOC2" << std::endl;
			//second order correlation is stored in soc after the flcf
			tom::os3_soc(&pattern,&soc,&mask);
			if(printStatus)
				std::cout << "PSR" << std::endl;
			//calculate the psr
			tom::os3_psr(&pattern,&psrV);
			if(printStatus)
				std::cout << "UPDATE" << std::endl;
			//update the result volume
			try{
				res.update(pattern.getVolume(),psrV.getVolume(),soc.getVolume(),iterationCounter);
			}catch(...){
				std::cerr << " Error during update of result object.  " << std::endl;
			}
			pattern.setVolume(&patternOriginal);
			if(printStatus)
				std::cerr << ".";
			iterationCounter++;
		}
		if(printStatus)
			std::cout << "END" << std::endl;

		try{
			if(job.filterSave){
				res.writeToEm(job);
				if(printStatus)
					std::cout << "Filter save" << std::endl;

			}
			if(job.listSave){
				if(printStatus)
					std::cout << "list save" << std::endl;
				res.setPicklist(job,printStatus);
				tom::os3_picklist<T> picklist = res.getPicklist();
				picklist.savePicklist(job.resultPath +std::string("pl.txt"));
			}
			
		}catch(...){
			std::cerr << "An error occured while saving results of job " << job.id<< ". Check your options!" << std::endl;
		}
		if(printStatus)
			std::cerr << "|" << job.id <<"|";
	}

/****************************************************************************//**
 * \brief Creates a particle stack based on correlation results.
 *
 * \param[in] job The job structure. Specifies correlation result path, filter weighting ...
 * \param[in] printStatus If true, debug output will appear on the screen.
 * \param[in] picklist
 * 
 * 
 * The first parameter (picklist) will be filled with particles extracted out of the volume specified by job.
 *******************************************************************************/
	template<typename T> 	
	void os3_createPicklist(tom::os3_job &job,bool printStatus){
		if(printStatus)
			jobInfo(job);

		tom::os3_results<T> results(job);
		results.setPicklist(job,printStatus);
		tom::os3_picklist<T> picklist;
		picklist = results.getPicklist();
		picklist.savePicklist(job.resultPath +std::string("pl.txt"));
		
	}

	template<typename T> void os3_picker(tom::os3_job &job,bool printStatus){
	
		if(printStatus)
			tom::jobInfo(job);

		try{
			switch(job.jobType){
			
				case(__TOM_OS3_CONTAINER_JOB__):{
					createResultContainers<T>(job);
					break;
				}
				case(__TOM_OS3_CORRELATION_JOB__):{
					os3_correlationPicker<T>(job,printStatus);
					break;
				}
				case(__TOM_OS3_PICKLIST_JOB__):{
					//tom::os3_picklist<T> picklist;
					os3_createPicklist<T>(job,printStatus);
					break;
				}
				case(__TOM_OS3_CLASSIFICATION_JOB__):{
					break;
				}
				case(__TOM_OS3_NO_JOB__):{
					break;
				}
				default:
					break;
			}
		}catch(...){
			throw;
		}
	}

}

template void tom::os3_correlationPicker<TFLOAT>(tom::os3_job &job,bool printStatus);
template void tom::os3_correlationPicker<TDOUBLE>(tom::os3_job &job,bool printStatus);

template void tom::os3_createPicklist<TFLOAT>(tom::os3_job &job,bool printStatus);
template void tom::os3_createPicklist<TDOUBLE>(tom::os3_job &job,bool printStatus);

template void tom::os3_picker<TFLOAT>(tom::os3_job& job,bool printStatus);
template void tom::os3_picker<TDOUBLE>(tom::os3_job& job,bool printStatus);




