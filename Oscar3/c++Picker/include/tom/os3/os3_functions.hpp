#ifndef ___INLCUDE_TOM_OS3_FUNCTIONS_H__
#define ___INLCUDE_TOM_OS3_FUNCTIONS_H__

#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_volume.hpp"

#include "tom/core/volume.hpp"

namespace tom{
	/*from os3_correlation.cpp*/
	template<typename T> void os3_correlate(tom::os3_volume<T>* searchV, tom::os3_volume<T>* patternV,const bool normalise);
	template<typename T> void os3_flcf(tom::os3_volume<T>* searchV,tom::os3_volume<T>* patternV, tom::os3_volume<T>* maskV);
	template<typename T> void os3_psr(tom::os3_volume<T>* xcfV,tom::os3_volume<T>* patternV);
	template<typename T> void os3_soc(tom::os3_volume<T>* corrV,tom::os3_volume<T>* autoV, tom::os3_volume<T>* maskV);
	template<typename T> void os3_apply_weight_function(tom::os3_volume<T>& volume,const tom::os3_volume<T>& psf);	

	/*from os3_organisation.cpp*/
	tom::os3_job* checkForFinishedVolumes(std::vector<tom::os3_job>& finishedJobList);
	std::vector<tom::angleTupel> generateAngleList(tom::angleTupel angleStart,tom::angleTupel angleIncrement,tom::angleTupel angleEnd);
	
	std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > split2d(tom_io_em_header* volumeHeader,tom_io_em_header* patternHeader,tom::st_idx splitFactor);
	std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > split3d(tom_io_em_header* volumeHeader,tom_io_em_header* patternHeader,tom::st_idx splitFactor);


	template<typename T> void createResultContainers(tom::os3_job & job);
	std::string appendStringToFilename(std::string filename,std::string str);

	void setJobSubregion(tom::os3_job& job,tom::subVolumeCoordinates& coordinates,tom::subVolumeCoordinates& writeBackCoordinates);
	void setJobSampling(tom::os3_job& job,uint32_t* sampling);
	void setJobBinning(tom::os3_job& job,uint32_t* binning);
	void setJobFilterWeight(tom::os3_job& job,double* weight);
	void jobInfo(tom::os3_job& job);
}

#endif




