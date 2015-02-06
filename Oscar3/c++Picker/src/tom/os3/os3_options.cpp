/****************************************************************************//**
 * \file tom_os3_options.cpp
 * \brief The source file for the class os3_Options.
 * \author  Thomas Hrabe 
 * \version 0.1
 * \date    08.01.2008
 *******************************************************************************/ 
#include "tom/os3/os3_functions.hpp"
#include "tom/os3/os3_options.hpp"
#include "tom/os3/os3_types.hpp"
/*----------------------------*/
#include "tom/core/io.h"
/*----------------------------*/
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>
#include <helper/filesystem.hpp>
#include <boost/lexical_cast.hpp>

namespace tom {

/****************************************************************************//**
 * \brief The standard constructor
 *
 * The standard constructor 
 *******************************************************************************/
	tom::os3_options::os3_options(){
	
	}


/****************************************************************************//**
 * \brief The standard destructor
 *
 * The standard destructor 
 *******************************************************************************/

	tom::os3_options::~os3_options(){
	
	}

/****************************************************************************//**
 * \brief Read an options file and return an options object.
 *
 * \param[in] *file absolute path to file
 * \param[in] *pointer to os3_options object
 *
 *******************************************************************************/

void tom::os3_options::readDirectories(){
	this->volumeList   = this->readFileDirectory(this->volumeDir);
	this->patternList  = this->readFileDirectory(this->patternDir);
	this->maskList     = this->readFileDirectory(this->maskDir);
	this->psfList 	   = this->readFileDirectory(this->psfDir);
}


/****************************************************************************//**
 * \brief Read a direcory containing a some volumes / images.
 *
 * \param[in] dir The absolute path to the directory.
 *
 * Returns a list of all em files in the directory.
 *******************************************************************************/
inline std::vector<std::string> tom::os3_options::readFileDirectory(std::string dir){

		//if(dir.length() == 0){
		//	throw std::runtime_error("os3_options::readFileDirectory - string to directory is empty.\n");
		//}
	
		std::cout << dir << std::endl;

		struct dirent **dirEntries;
		int count;
	
		count = scandir(dir.c_str(), &dirEntries, 0, 0);
	
		std::vector<std::string> files;
		files.clear();
		std::string file;
	
		while(count--) {
			file = std::string((*dirEntries)->d_name);
			
			if(! ( file.length() == 1 || file.length() == 2 ) && file.find(std::string(".em"))){
				files.push_back(dir + "/" + file);
			}
	
			free(*dirEntries++);
		}

    return files;

}


/****************************************************************************//**
 * \brief Sets the volume read out values.
 *
 * Sets subregion read and binning to reasonable values.
 *******************************************************************************/
void tom::os3_options::setVolumeProperties(short split_x,short split_y,short split_z,uint32_t bin_x,uint32_t bin_y,uint32_t bin_z){

	this->split.x = split_x;
	this->split.y = split_y;
	this->split.z = split_z;

	this->binning[0] = bin_x;
	this->binning[1] = bin_y;
	this->binning[2] = bin_z;

}

/****************************************************************************//**
 * \brief Creates a joblist.
 *
 * Creates a list of jobs for the picker.
 *******************************************************************************/
std::vector<tom::os3_job> tom::os3_options::createJobList(){

	//for each volume in the volumeList 

	std::vector<std::string>::iterator volumeIterator;
	std::vector<tom::os3_job> jobs; 
	if(volumeList.size() == 0)
		this->volumeList   = this->readFileDirectory(this->volumeDir);

	if(patternList.size() == 0)
		this->patternList   = this->readFileDirectory(this->patternDir);

	if(maskList.size() == 0)
		this->maskList   = this->readFileDirectory(this->maskDir);

	if(psfList.size() == 0)
		this->psfList   = this->readFileDirectory(this->psfDir);

	for(volumeIterator = volumeList.begin();volumeIterator<volumeList.end();volumeIterator++){
		tom::os3_job job;
		std::cout << *volumeIterator << std::endl;
		job.volumePath = *volumeIterator;
	
		std::size_t pos = (*volumeIterator).find_last_of("/",(*volumeIterator).length());
		std::string  volumeName = (*volumeIterator).substr(pos+1,(*volumeIterator).length());
		job.resultPath = this->resultDir+volumeName;
		
		std::vector<std::string>::iterator patternIterator;	

		for(patternIterator= patternList.begin();patternIterator<patternList.end();patternIterator++){

			job.patternPath = *patternIterator; 

			pos = (*patternIterator).find_last_of("/",(*patternIterator).length());

			std::string patternName = (*patternIterator).substr(pos+1,(*patternIterator).length());
			
			bool found = false;
			
			std::vector<std::string>::iterator maskIterator = maskList.begin();	
			
			while(!found && maskIterator < maskList.end()){
				found = (*maskIterator).find(patternName);
			}
		
			if(!found)
				throw std::runtime_error("\ntom::os3_options::createJobList() - Could not find mask for pattern.");
			
			job.maskPath = *maskIterator;
			
			std::vector<std::string>::iterator psfIterator = psfList.begin();	
			found = false;
			while(!found && psfIterator < psfList.end()){
				found = (*psfIterator).find(patternName);
			}
			
			if(!found){
				std::cout << std::endl<< "\ntom::os3_options::createJobList() - Warning: could not find a point spread function for the pattern " << patternName << "!\nNo modification will be applied to the pattern during picking." << std::endl<< std::endl;
				job.psfPath = std::string("NONE");
			}
			else
				job.psfPath = *psfIterator;
		}

		job.filterSave = this->filterSave;
		job.listSave = this->listSave;
		job.stackSave = this->stackSave;
		
		//is the prototype for a real job. a container job calls function xyz and allocates result containers on disk for the correlation
		std::vector<tom::os3_job> subJobList = applyReadoutOptions(job);
		std::vector<tom::os3_job>::iterator subJobIterator;

		for(subJobIterator = subJobList.begin();subJobIterator < subJobList.end();subJobIterator++)
			jobs.push_back(*subJobIterator);
	}

	return jobs;

}

/****************************************************************************//**
 * \brief Creates an oscar3 directory structure at given path.
 *
 * Creates an oscar3 directory structure at given path.
 *******************************************************************************/
void tom::os3_options::createDirectories(std::string path){

	boost::filesystem::path p(path);

	if(!exists(p))
		throw std::runtime_error("\n tom::os3_picker::createDirectories() - Wrong path?");
	
	boost::filesystem::create_directory( p / "volumes" );
	boost::filesystem::create_directory( p / "pattern" );
	boost::filesystem::create_directory( p / "mask" );
	boost::filesystem::create_directory( p / "results" );
	boost::filesystem::create_directory( p / "psf" );
	boost::filesystem::create_directory( p / "options" );

	createOptionsFile(p / "options" /"options.txt");
}

/****************************************************************************//**
 * \brief Creates a template for a options file.
 *
 * Creates a template for a options file.
 *******************************************************************************/
void tom::os3_options::createOptionsFile(boost::filesystem::path path){

	std::ofstream fileOut;

	fileOut.open(path.file_string().c_str());

	if(fileOut.bad())
		throw std::runtime_error("\n tom::os3_picker::os3_savePicklist() - Error saving picklist! Wrong path?");
	
	fileOut << "#OSCAR 3 options file" << std::endl;
	fileOut << "#set every value. Should you not need it, set it to zero." << std::endl;
	fileOut << "#volumes directory: (string)" << std::endl;
	fileOut << "???" <<std::endl;
	fileOut << "#pattern directory: (string)" << std::endl;
	fileOut << "???" <<std::endl;
	fileOut << "#mask directory: (string)" << std::endl;
	fileOut << "???" <<std::endl;
	fileOut << "#psf directory: (string)" << std::endl;
	fileOut << "???" <<std::endl;
	fileOut << "#result directory: (string)" << std::endl;
	fileOut << "???" <<std::endl;

	fileOut << "#jobType: (string) - 2D / 3D" << std::endl;
	fileOut << "???" << std::endl;

	fileOut << "#filter weighting: (scalar) - XCF, PSR, SOC" << std::endl;
	fileOut << "0,0,0" << std::endl;

	fileOut << "#binning: (scalar) - x , y , z" << std::endl;
	fileOut << "0,0,0"<< std::endl;
	
	fileOut << "#subregion: (scalar) - x , y, z - [if 2D, set x value only and the image will be split into x^2 subimages. Set y and z to 0!]" << std::endl;
	fileOut << "0,0,0" << std::endl;

	fileOut << "#angles start: (scalar) - phi , psi , theta" << std::endl;
	fileOut << "0,0,0" << std::endl;
	
	fileOut << "#angles increment: (scalar) - phi , psi , theta" << std::endl;
	fileOut << "0,0,0" << std::endl;

	fileOut << "#angles end: (scalar) - phi , psi , theta" << std::endl;
	fileOut << "0,0,0" << std::endl;

	fileOut << "#number of particles per image: (scalar) " << std::endl;
	fileOut << "0" << std::endl;

	fileOut << "#save picklist to disk?: (boolean) - 1 for yes, 0 for no" << std::endl;
	fileOut << "0" << std::endl;

	fileOut << "#save filter results to disk?: (boolean) - 1 for yes, 0 for no" << std::endl;
	fileOut << "0" << std::endl;

	fileOut << "#save particle stack to disk?: (boolean) - 1 for yes, 0 for no" << std::endl;
	fileOut << "0" << std::endl;

	fileOut.close();

}


namespace{
	
	std::vector<std::string> splitString(std::string str, std::string splitChar){

		std::vector<std::string> stringV;
		std::string::size_type loc;
		while( (loc = str.find_first_of(splitChar)) != std::string::npos){
			
			stringV.push_back(str.substr(0,loc));
			str = str.substr(loc+1);
		}
		stringV.push_back(str);
		return stringV;

	}
}
/****************************************************************************//**
 * \brief Reads options from file.
 *
 * Reads options from file.
 *******************************************************************************/
void tom::os3_options::readFromFile(std::string path){

	std::ifstream fileIn;
			

	fileIn.open(path.c_str());

	if(fileIn.bad())
		throw std::runtime_error("\n tom::os3_options::readFromFile() - Error reading options! Wrong path?");

	std::string line;

	while(! fileIn.eof()){
		
		std::getline(fileIn,line);
		
		if( line.find("#volumes directory") != std::string::npos ){
			std::getline(fileIn,line);
			this->volumeDir = line;
		}
		
		if( line.find("#pattern directory") != std::string::npos ){
			std::getline(fileIn,line);
			this->patternDir = line;
		}
		
		if( line.find("#mask directory") != std::string::npos ){
			std::getline(fileIn,line);
			this->maskDir = line;
		}
		
		if( line.find("#psf directory") != std::string::npos ){
			std::getline(fileIn,line);
			this->psfDir = line;
		}
		
		if( line.find("#result directory") != std::string::npos ){
			std::getline(fileIn,line);
			this->resultDir= line;
		}
		
		if( line.find("#filter weighting") != std::string::npos ){
			std::getline(fileIn,line);
			std::vector<std::string> numbers = splitString(line, ",");
			this->weight[0] = boost::lexical_cast<double>(numbers[0]);
			this->weight[1] = boost::lexical_cast<double>(numbers[1]);
			this->weight[2] = boost::lexical_cast<double>(numbers[2]);
		}

		if( line.find("#binning") != std::string::npos ){
			std::getline(fileIn,line);
			std::vector<std::string> numbers = splitString(line, ",");
			this->binning[0] = boost::lexical_cast<uint32_t>(numbers[0]);
			this->binning[1] = boost::lexical_cast<uint32_t>(numbers[1]);
			this->binning[2] = boost::lexical_cast<uint32_t>(numbers[2]);
		}
		

		if( line.find("#subregion") != std::string::npos ){
			std::getline(fileIn,line);
			std::vector<std::string> numbers = splitString(line, ",");
			
			this->split.x = boost::lexical_cast<std::size_t>(numbers[0]);
			this->split.y = boost::lexical_cast<std::size_t>(numbers[1]);
			this->split.z = boost::lexical_cast<std::size_t>(numbers[2]);
		}

		if( line.find("#angles start") != std::string::npos ){
			std::getline(fileIn,line);
			std::vector<std::string> numbers = splitString(line, ",");
			
			this->anglesStart.phi = boost::lexical_cast<double>(numbers[0]);
			this->anglesStart.psi = boost::lexical_cast<double>(numbers[1]);
			this->anglesStart.theta = boost::lexical_cast<double>(numbers[2]);
		}

		if( line.find("#angles increment") != std::string::npos ){
			std::getline(fileIn,line);
			std::vector<std::string> numbers = splitString(line, ",");
			
			this->anglesIncrement.phi = boost::lexical_cast<double>(numbers[0]);
			this->anglesIncrement.psi = boost::lexical_cast<double>(numbers[1]);
			this->anglesIncrement.theta = boost::lexical_cast<double>(numbers[2]);
		}
		
		if( line.find("#angles end") != std::string::npos ){
			std::getline(fileIn,line);
			std::vector<std::string> numbers = splitString(line, ",");
			
			this->anglesEnd.phi = boost::lexical_cast<double>(numbers[0]);
			this->anglesEnd.psi = boost::lexical_cast<double>(numbers[1]);
			this->anglesEnd.theta = boost::lexical_cast<double>(numbers[2]);
		}

		if( line.find("#number of particles") != std::string::npos ){
			std::getline(fileIn,line);
			this->numberOfParticles = boost::lexical_cast<std::size_t>(line);
		}

		if( line.find("#save filter") != std::string::npos ){
			std::getline(fileIn,line);
			this->filterSave= boost::lexical_cast<bool>(line);
		}

		if( line.find("#save picklist") != std::string::npos ){
			std::getline(fileIn,line);
			this->listSave= boost::lexical_cast<bool>(line);
		}

		if( line.find("#save particle") != std::string::npos ){
			std::getline(fileIn,line);
			this->stackSave= boost::lexical_cast<bool>(line);
		}
	}
	
	fileIn.close();

}

/****************************************************************************//**
 * \brief Splits one job into multiple jobs.
 *
 * Splits one job into multiple jobs according to the read out settings.
 *******************************************************************************/
inline std::vector<tom::os3_job> tom::os3_options::applyReadoutOptions(tom::os3_job job){

	tom_io_em_header volumeHeader,patternHeader;
	std::vector<tom::os3_job> subJobs;
	
	if(tom_io_em_read_header(job.volumePath.c_str(), NULL, &volumeHeader,NULL) != TOM_ERR_OK)
			throw std::runtime_error("\n tom::os3_options::applyReadoutOptions() - Error reading volume file.");
	if(tom_io_em_read_header(job.patternPath.c_str(), NULL, &patternHeader,NULL))
			throw std::runtime_error("\n tom::os3_options::applyReadoutOptions() - Error reading pattern file.");

	//create a special job which creates the result contrainers on disk
	tom::os3_job containerJob;
	containerJob.jobType = __TOM_OS3_CONTAINER_JOB__;
	containerJob.volumePath = job.volumePath;
	containerJob.resultPath = job.resultPath;	
	
	tom::os3_job newJob;
	newJob.jobType = __TOM_OS3_CORRELATION_JOB__;

	newJob.volumePath = job.volumePath;
	newJob.patternPath = job.patternPath;
	newJob.resultPath = job.resultPath;
	newJob.maskPath   = job.maskPath;
	newJob.psfPath	  = job.psfPath;

	if(this->split.x >0 || this->split.y >0 || this->split.z >0){
		//list of subregion coordinates with extension (size of pattern volume)
		std::vector<tom::subVolumeCoordinates> coordinates;
		//list of subregion coordinates without extension 
		std::vector<tom::subVolumeCoordinates> writeBackCoordinates;
		std::pair<std::vector<tom::subVolumeCoordinates> ,std::vector<tom::subVolumeCoordinates> > p;
		if(volumeHeader.dims[2] == 1)
		//2d case 
			p= tom::split2d(&volumeHeader,&patternHeader,this->split);
		else
		//3d case
			p= tom::split3d(&volumeHeader,&patternHeader,this->split);
		
		coordinates = p.first;
		writeBackCoordinates = p.second;
		
		//append __TOM_OS3_CONTAINER_JOB__ to list if the settings cause a volume splitting
		if(coordinates.size() > 1)
			subJobs.push_back(containerJob);
		
		std::vector<tom::angleTupel> angles = tom::generateAngleList(this->anglesStart,this->anglesIncrement,this->anglesEnd);
		
		std::vector<tom::subVolumeCoordinates>::iterator coordinatesIterator;
		std::vector<tom::subVolumeCoordinates>::iterator writeBackIterator = writeBackCoordinates.begin();
		
		std::size_t counter = 0;
		std::size_t numberOfJobs = coordinates.size();

		uint32_t zeros[6] = {0,0,0};
		tom::setJobSampling(newJob,zeros);
		tom::setJobBinning(newJob,this->binning);

		for(coordinatesIterator = coordinates.begin();coordinatesIterator < coordinates.end();coordinatesIterator++){
		
			tom::setJobSubregion(newJob,*coordinatesIterator,*writeBackIterator++);
			newJob.id =counter++;
			newJob.numberJobs = numberOfJobs;
			newJob.angleList = angles;
			subJobs.push_back(newJob);
		}
		
	}else{
		uint32_t zeros[6] = {0,0,0,(uint32_t)volumeHeader.dims[0],(uint32_t)volumeHeader.dims[1],(uint32_t)volumeHeader.dims[2]};
		//tom::setJobSubregion(newJob,*coordinatesIterator,*writeBackIterator++);
		tom::setJobSampling(newJob,zeros);
		tom::setJobBinning(newJob,this->binning);	
		//job.setNumberJobs(numberOfJobs);
		job.jobType =__TOM_OS3_CORRELATION_JOB__ ;
		subJobs.push_back(newJob);
	}

	return subJobs;
}

}










