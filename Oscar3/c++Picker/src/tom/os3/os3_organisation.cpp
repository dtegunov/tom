/***********************************************************************//**
 * \file tom_os3_organisation.cpp
 * \brief Some helper functions used to generate jobs.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    05.02.2008
 **************************************************************************/

#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_functions.hpp"


#include "tom/core/io.h"




namespace tom{

/****************************************************************************//**
 * \brief Generates a list of tom::angleTupel.
 *
 * \param[in] angleStart The start destination. 
 * \param[in] angleIncrement The increment along each angle.
 * \param[in] angleEnd The end value until each angle is incremented.
 * 
 * Generates a list of tom::angleTupel used for rotating the pattern searched.
 *******************************************************************************/
std::vector<tom::angleTupel> generateAngleList(tom::angleTupel angleStart,tom::angleTupel angleIncrement,tom::angleTupel angleEnd){

	std::vector<tom::angleTupel> angleList;

	tom::angleTupel angleTupel;

	for(double phi=angleStart.phi;phi<angleEnd.phi;phi += angleIncrement.phi)
		for(double psi=angleStart.psi;psi<angleEnd.psi;psi += angleIncrement.psi)
			for(double  theta=angleStart.theta;theta<angleEnd.theta;theta += angleIncrement.theta){

				angleTupel.phi = PI * phi/180;
				angleTupel.psi = PI * psi/180;
				angleTupel.theta = PI * theta/180;

				angleList.push_back(angleTupel);

			}
	return angleList;
}


/****************************************************************************//**
 * \brief Returns a list of coordinates used for splitting the image.
 *
 * \param[in] volumeHeader The header of the volume.
 * \param[in] patternHeader The header of the pattern searched.
 * \param[in] splitFactor The number of subvolumes along the x and y axis.
 * 
 * Returns a list of coordinates used for splitting the image (refered to as volume)
 *******************************************************************************/
std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > split2d(tom_io_em_header* volumeHeader,tom_io_em_header* patternHeader,tom::st_idx split){

	uint32_t volumeX = volumeHeader->dims[0];
	uint32_t volumeY = volumeHeader->dims[1];

	uint32_t patternX = patternHeader->dims[0]-1;
	uint32_t patternY = patternHeader->dims[1]-1;	

	std::size_t splitFactor = split.x;

	uint32_t xIncrement = volumeX/splitFactor;
	uint32_t yIncrement = volumeY/splitFactor;
	uint32_t lastX = 0;
	uint32_t lastY = 0;
	uint32_t offsetX = 0;
	uint32_t offsetY = 0;
	//list of subregion coordinates with extension (size of pattern volume)
	std::vector<tom::subVolumeCoordinates> coordinateList;
	//list of subregion coordinates without extension
	std::vector<tom::subVolumeCoordinates> writeBackList;	

	for(uint32_t x =(xIncrement-1);x<volumeX;x+=xIncrement){
		for(uint32_t  y=(yIncrement-1);y<volumeY;y+=yIncrement){
			tom::subVolumeCoordinates coordinates;
			tom::subVolumeCoordinates writeBackCoordinates;		

			offsetX = 0;
			offsetY = 0;
			
			if(lastX > 0){
				lastX = lastX - patternX;
				offsetX = patternX;
			}

			if(lastY > 0){
				lastY = lastY - patternY;
				offsetY = patternY;
			}

			writeBackCoordinates.startX = offsetX;
			writeBackCoordinates.startY = offsetY;
			writeBackCoordinates.startZ = 0;

			coordinates.startX = lastX;
			coordinates.startY = lastY;
			coordinates.startZ = 0;
			
			lastX += offsetX;
			lastY += offsetY;

			if(x + patternX > volumeX-1)
				coordinates.noVoxelsX = offsetX + xIncrement; 
			else
				coordinates.noVoxelsX = offsetX + xIncrement + patternX;
				
			if(y + patternY> volumeY-1)
				coordinates.noVoxelsY = offsetY + yIncrement; 
			else
				coordinates.noVoxelsY = offsetY + yIncrement + patternY; 

			coordinates.noVoxelsZ = 1;


			writeBackCoordinates.noVoxelsX = xIncrement;
			writeBackCoordinates.noVoxelsY = yIncrement;
			writeBackCoordinates.noVoxelsZ = 1;
	
			writeBackList.push_back(writeBackCoordinates);
			coordinateList.push_back(coordinates);
	
			lastY = y+1;

		}
		lastY = 0;
		lastX = x+1;
	}

	std::pair<std::vector<tom::subVolumeCoordinates> ,std::vector<tom::subVolumeCoordinates> > p(coordinateList,writeBackList);
	
	return p;
}


/****************************************************************************//**
 * \brief Returns a list of coordinates used for splitting the volume.
 *
 * \param[in] volumeHeader The header of the volume.
 * \param[in] patternHeader The header of the pattern searched.
 * \param[in] splitFactor The number of subvolumes along the x and y axis.
 * 
 * Returns a list of coordinates used for splitting the image (refered to as volume)
 *******************************************************************************/
std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > split3d(tom_io_em_header* volumeHeader,tom_io_em_header* patternHeader,tom::st_idx subVolumeSize){
	
	uint32_t volumeX = volumeHeader->dims[0];
	uint32_t volumeY = volumeHeader->dims[1];
	uint32_t volumeZ = volumeHeader->dims[2];
//  	std::cout << "vol(" << volumeX << "," << volumeY<< "," << volumeZ << ")"<< std::endl;

	uint32_t patternX = patternHeader->dims[0]-1;
	uint32_t patternY = patternHeader->dims[1]-1;	
	uint32_t patternZ = patternHeader->dims[2]-1;	
//  	std::cout << "pat(" << patternX << "," << patternY<< "," << patternZ<< ")"<< std::endl;

	uint32_t xIncrement = subVolumeSize.x;
	uint32_t yIncrement = subVolumeSize.y;
	uint32_t zIncrement = subVolumeSize.z;
//  	std::cout << "inc(" << xIncrement << "," << yIncrement<< "," << zIncrement << ")"<< std::endl;	
	uint32_t lastX = 0;
	uint32_t lastY = 0;
	uint32_t lastZ = 0;
/*	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;	*/
	uint32_t offsetX = 0;
	uint32_t offsetY = 0;
	uint32_t offsetZ = 0;
	
	//list of subregion coordinates with extension (size of pattern volume)
	std::vector<tom::subVolumeCoordinates> coordinateList;
	//list of subregion coordinates without extension
	std::vector<tom::subVolumeCoordinates> writeBackList;	
	tom::subVolumeCoordinates coordinates;
	tom::subVolumeCoordinates writeBackCoordinates;		
	//
	for(uint32_t x=0;x<volumeX;x+=xIncrement){
// 		std::cout << "(X)" << x << std::endl;
		for(uint32_t y=0;y<volumeY;y+=yIncrement){
// 			std::cout << "(Y)" << y << std::endl;
			for(uint32_t z=0;z<volumeZ;z+=zIncrement){
// 				std::cout << "(Z)" << z << std::endl;
	


				offsetX = 0;
				offsetY = 0;
				offsetZ = 0;
				
				if(lastX > 0){
					lastX = lastX - patternX;
					offsetX = patternX;
				}
				if(lastY > 0){
					lastY = lastY - patternY;
					offsetY = patternY;
				}
				if(lastZ > 0){
					lastZ = lastZ - patternZ;
					offsetZ= patternZ;
				}
				
				//if an overlap to the left is possible, then save the offset of the "real" volume area in writeBackCoordinates
				//offset is the number of voxels of the "real" volume from (0,0,0)
				writeBackCoordinates.startX = offsetX;
				writeBackCoordinates.startY = offsetY;
				writeBackCoordinates.startZ = offsetZ;
				
				coordinates.startX = lastX;
				coordinates.startY = lastY;
				coordinates.startZ = lastZ;
				
				lastX += offsetX;
				lastY += offsetY;
				lastZ += offsetZ;
				
				if(lastX + xIncrement > volumeX-1){
					uint32_t rest = volumeX - lastX;
					coordinates.noVoxelsX = offsetX + rest; 
					writeBackCoordinates.noVoxelsX = rest;
				}
				else if(lastX + xIncrement + patternX > volumeX-1){
					uint32_t rest = volumeX - (lastX + xIncrement);
					coordinates.noVoxelsX = offsetX + xIncrement + rest; 
					writeBackCoordinates.noVoxelsX = xIncrement;
				}
				else{
					coordinates.noVoxelsX = offsetX + xIncrement + patternX;
					writeBackCoordinates.noVoxelsX = xIncrement;
				}


				if(lastY + patternY> volumeY-1){
					uint32_t rest = volumeY - lastY;
					coordinates.noVoxelsY = offsetY + rest; 
					writeBackCoordinates.noVoxelsY = rest;
				}
				else if(lastY + yIncrement + patternY > volumeY-1){
					uint32_t rest = volumeY - (lastY + yIncrement);
					coordinates.noVoxelsY = offsetY + yIncrement + rest; 
					writeBackCoordinates.noVoxelsY = yIncrement;
				}
				else{
					coordinates.noVoxelsY = offsetY + yIncrement + patternY; 
					writeBackCoordinates.noVoxelsY = yIncrement;
				}


				if(lastZ + patternZ> volumeZ-1){
					uint32_t rest = volumeZ - lastZ;
					coordinates.noVoxelsZ = offsetZ + rest; 
					writeBackCoordinates.noVoxelsZ = rest;
				}
				else if(lastZ + zIncrement + patternZ > volumeZ){
					uint32_t rest = volumeZ - (lastZ + zIncrement);
					coordinates.noVoxelsZ = offsetZ + zIncrement+rest; 
					writeBackCoordinates.noVoxelsZ = zIncrement;
				}
				else{
					coordinates.noVoxelsZ = offsetZ + zIncrement + patternZ; 
					writeBackCoordinates.noVoxelsZ = zIncrement;
				}
				
				writeBackList.push_back(writeBackCoordinates);
				coordinateList.push_back(coordinates);

//   				std::cout << "last(" << lastX << "," << lastY<< "," << lastZ << ")"<< std::endl;
//   				std::cout << "c(" << coordinates.startX << "," << coordinates.startY<< "," << coordinates.startZ << ")";
//   				std::cout << "(" <<coordinates.noVoxelsX << "," << coordinates.noVoxelsY<< "," << coordinates.noVoxelsZ << ")" << std::endl;
//   				std::cout << "wb(" << writeBackCoordinates.startX << "," << writeBackCoordinates.startY<< "," << writeBackCoordinates.startZ << ")";
//   				std::cout <<  "(" <<writeBackCoordinates.noVoxelsX << "," << writeBackCoordinates.noVoxelsY<< "," << writeBackCoordinates.noVoxelsZ << ")" << std::endl<< std::endl;

				lastZ = coordinates.startZ + writeBackCoordinates.startZ + writeBackCoordinates.noVoxelsZ;
			}
			lastZ = 0;
			lastY = coordinates.startY + writeBackCoordinates.startY + writeBackCoordinates.noVoxelsY;
		}
		lastY = 0;
		lastX = coordinates.startX + writeBackCoordinates.startX + writeBackCoordinates.noVoxelsX;
	}

	std::pair<std::vector<tom::subVolumeCoordinates> ,std::vector<tom::subVolumeCoordinates> > p(coordinateList,writeBackList);
	
	return p;

}

std::string appendStringToFilename(std::string filename,std::string str){

	std::size_t index = filename.find_last_of(".em");

	filename = filename.substr(0,index-2);
	
	return std::string(filename + str);

}

template<typename T> void createResultContainers(tom::os3_job& job){
	
	std::string resultName = job.resultPath;
	
	tom_io_em_header volumeHeader;
	
	if(tom_io_em_read_header(job.volumePath.c_str(), NULL, &volumeHeader,NULL) != TOM_ERR_OK)
			throw std::runtime_error("\n tom::os3_options::applyReadoutOptions() - Error reading volume file.");

	tom::Volume<T> volume(volumeHeader.dims[0],volumeHeader.dims[1],volumeHeader.dims[2],NULL,NULL);
	volume.setValues(0);
	
	volume.write_to_em(appendStringToFilename(resultName,"XCF.em"),NULL);
	volume.write_to_em(appendStringToFilename(resultName,"PSR.em"),NULL);
	volume.write_to_em(appendStringToFilename(resultName,"SOC.em"),NULL);
	tom::Volume<int> volume2(volumeHeader.dims[0],volumeHeader.dims[1],volumeHeader.dims[2],NULL,NULL);
	volume2.write_to_em(appendStringToFilename(resultName,"ANG.em"),NULL);

	
}


/****************************************************************************//**
 * \brief Determines whether a volume has been entirely correlated. 
 *
 * \param[in] finishedJobList List of finished jobs.
 *
 * Determines whether a volume has been entirely correlated. If so, the sum of all job ids must be equal to the 
 * gaussian sum of the particular job numbers.
 *******************************************************************************/
tom::os3_job* checkForFinishedVolumes(std::vector<tom::os3_job>& finishedJobList){
	std::vector<tom::os3_job>::iterator jobIterator = finishedJobList.begin();
	std::vector<tom::os3_job> volumeList;
	std::vector<std::string> visitedList;
	
	tom::os3_job* returnJob = NULL;
	std::auto_ptr<tom::os3_job> autJob(returnJob);
	returnJob = new tom::os3_job();
	returnJob->jobType =__TOM_OS3_NO_JOB__;
	while((*jobIterator).jobType != __TOM_OS3_CORRELATION_JOB__ && jobIterator < finishedJobList.end())
		jobIterator++;

	for(;jobIterator<finishedJobList.end();jobIterator++){

		std::string currentVolume = (*jobIterator).volumePath;
		std::vector<std::string>::iterator visitedIterator;
		
		/*check whether the job has been visited already*/
		visitedIterator=visitedList.begin();
		bool found = false;
		while(visitedIterator < visitedList.end() && !found){
			found = currentVolume.compare(*visitedIterator++) == 0 || (*jobIterator).jobType == __TOM_OS3_CORRELATION_JOB__;
		}

		if(!found){
			visitedList.push_back(currentVolume);
			
			std::size_t currentSum = 0;
			
			std::vector<tom::os3_job>::iterator sumIterator;
// 			std::cout << currentVolume << std::endl;	
	
			for(sumIterator = finishedJobList.begin();sumIterator<finishedJobList.end();sumIterator++){
				if((currentVolume.compare((*sumIterator).volumePath) ==0))
					currentSum += (*sumIterator).jobType;
			}
			
			if((*jobIterator).numberJobs == currentSum){	
				returnJob->volumePath=((*jobIterator).volumePath);
				returnJob->patternPath=((*jobIterator).patternPath);
				returnJob->resultPath = ((*jobIterator).resultPath);
				returnJob->jobType = __TOM_OS3_PICKLIST_JOB__;
			}
		}
	}

	//remove from finished list
	autJob.release();
	return returnJob;
}

}

std::vector<tom::angleTupel> tom::generateAngleList(tom::angleTupel angleStart,tom::angleTupel angleIncrement,tom::angleTupel angleEnd);
std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > tom::split2d(tom_io_em_header* volumeHeader,tom_io_em_header* patternHeader,tom::st_idx split);
std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > tom::split3d(tom_io_em_header* volumeHeader,tom_io_em_header* patternHeader,tom::st_idx split);
tom::os3_job* tom::checkForFinishedVolumes(std::vector<tom::os3_job>& finishedJobList);
template void tom::createResultContainers<TFLOAT>(tom::os3_job& job);
template void tom::createResultContainers<TDOUBLE>(tom::os3_job& job);
tom::os3_job* checkForFinishedVolumes(std::vector<tom::os3_job>& finishedJobList);
std::string tom::appendStringToFilename(std::string filename,std::string str);





