#include "tom/os3/os3_structures.hpp"

#include <iostream>

namespace tom{
/****************************************************************************//**
 * \brief Sets the subregion coordinates of a particular job. 
 *
 * \param[in] job job reference.
 * \param[in] coordinates Coordinates of the subregion in the file, containing an extension of the pattern size.
 * \param[in] writeBackCoordinates Coordinates of the subregion without the extension. Point to the coordinates in the file.
 *
 * Sets the subregion coordinates of a particular job. 
 *******************************************************************************/
void setJobSubregion(tom::os3_job& job,tom::subVolumeCoordinates& coordinates,tom::subVolumeCoordinates& writeBackCoordinates){
	   	job.subregion[0] = (uint32_t) coordinates.startX;
		job.subregion[1] = (uint32_t) coordinates.startY;
		job.subregion[2] = (uint32_t) coordinates.startZ;
		job.subregion[3] = (uint32_t) coordinates.noVoxelsX;
		job.subregion[4] = (uint32_t) coordinates.noVoxelsY;
		job.subregion[5] = (uint32_t) coordinates.noVoxelsZ;

		job.writeBackSubregion[0] = (uint32_t) writeBackCoordinates.startX;
		job.writeBackSubregion[1] = (uint32_t) writeBackCoordinates.startY;
		job.writeBackSubregion[2] = (uint32_t) writeBackCoordinates.startZ;
		job.writeBackSubregion[3] = (uint32_t) writeBackCoordinates.noVoxelsX;
		job.writeBackSubregion[4] = (uint32_t) writeBackCoordinates.noVoxelsY;
		job.writeBackSubregion[5] = (uint32_t) writeBackCoordinates.noVoxelsZ;
}

void setJobSampling(tom::os3_job& job,uint32_t* extSampling){
	job.sampling[0] = extSampling[0];
	job.sampling[1] = extSampling[1];
	job.sampling[2] = extSampling[2];
}

void setJobBinning(tom::os3_job& job,uint32_t* extBinning){
    job.binning[0] = extBinning[0];
	job.binning[1] = extBinning[1];
	job.binning[2] = extBinning[2];
}

void setJobFilterWeight(tom::os3_job& job,double* weight){
	job.weight[0] = weight[0];
	job.weight[1] = weight[1];
	job.weight[2] = weight[2];
}

void jobInfo(tom::os3_job& job){
	std::cout << std::endl<< "*-----------------------------------------------------------*" << std::endl;
	std::cout << "Job ID : " << job.id << std::endl;
	std::cout << "Number of total jobs : " << job.numberJobs << std::endl;
	std::cout << "Volume: " << job.volumePath << std::endl;
	std::cout << "Pattern: " << job.patternPath << std::endl;
	std::cout << "maskPath: " << job.maskPath<< std::endl;
	std::cout << "psfPath: " << job.psfPath<< std::endl;
	std::cout << "resultPath: " << job.resultPath<< std::endl;
	std::cout << "number of rotations: " << job.angleList.size() << std::endl;
    std::cout << "Subregion:" << std::endl;
	std::cout << job.subregion[0] << "," << job.subregion[1] << ","<< job.subregion[2] << ","<< 
				 job.subregion[0]+ job.subregion[3] << ","<< job.subregion[1] +job.subregion[4] << ","<< job.subregion[2] + job.subregion[5] << std::endl;
	std::cout << "Write Back Subregion:" << std::endl;
	std::cout << job.writeBackSubregion[0] << "," << job.writeBackSubregion[1] << ","<< job.writeBackSubregion[2] << ","<< 
				 job.writeBackSubregion[3] << "," << job.writeBackSubregion[4] << ","<< job.writeBackSubregion[5] << std::endl;
	std::cout << "*-----------------------------------------------------------*" << std::endl;
}
}

void tom::setJobSubregion(tom::os3_job& job,tom::subVolumeCoordinates& coordinates,tom::subVolumeCoordinates& writeBackCoordinates);
void tom::setJobSampling(tom::os3_job& job,uint32_t* sampling);
void tom::setJobBinning(tom::os3_job& job,uint32_t* binning);
void setJobFilterWeight(tom::os3_job& job,double* weight);
void jobInfo(tom::os3_job& job);


