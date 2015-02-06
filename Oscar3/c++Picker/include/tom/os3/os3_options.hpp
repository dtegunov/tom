/****************************************************************************//**
 * \file tom_os3_options.hpp
 * \brief The header file for the class os3_options.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    08.01.2008
 *******************************************************************************/ 
#ifndef __INCLUDE_TOM_OS3_OPTIONS_HPP__
#define __INCLUDE_TOM_OS3_OPTIONS_HPP__

#include <stdio.h>
#include <dirent.h>

/*----------------------------*/
#include "tom/os3/os3_structures.hpp"
/*----------------------------*/
#include <vector>
#include <iostream>
#include <helper/filesystem.hpp>


namespace tom {

	class os3_options {
	
	public:
			
		os3_options();	

		~os3_options();

		void readDirectories();
	
		void setVolumeDir(std::string str);
		void setPatternDir(std::string str);
		void setMaskDir(std::string str);
		void setPsfDir(std::string str);
		void setAngles(tom::angleTupel start,tom::angleTupel inc,tom::angleTupel end);
		void setResultDir(std::string str);
		void setFilterWeight(double* weight);
		std::string getVolumeDir();
		std::string getPatternDir();
		std::string getMaskDir();
		std::string getPsfDir();
		std::string getResultDir();
		std::vector<tom::os3_job> createJobList();
		double* getFilterWeight();
		void setVolumeProperties(short split_x,short split_y,short split_z,uint32_t bin_x,uint32_t bin_y,uint32_t bin_z);

		void readFromFile(std::string path);

		static void createOptionsFile(boost::filesystem::path path);
		static void createDirectories(std::string path);
	protected:
	
	private:	

		std::vector<std::string> readFileDirectory(std::string dir);

		std::string jobType;	/* Define the job type - 2d or 3d */
		std::string jobName;	/* Set the jobName */

		std::string volumeDir; 		/* directory containing the images / volumes for PR */
		std::string resultDir;		/* target directory for the results */
		std::string patternDir;	/* path to the pattern  */
		std::string maskDir;		/* path to the mask used for XCF (FLCF - Roseman 2003) */
		std::string psfDir;			/* path to the psf used to manipulate the pattern */
		
		std::vector<std::string> volumeList;
		std::vector<std::string> patternList;
		std::vector<std::string> maskList;
		std::vector<std::string> psfList;


		tom::st_idx split;
		uint32_t binning[3];
		
		std::vector<tom::os3_job> applyReadoutOptions(tom::os3_job job);

		double weight[3];  	/* weightting of the three different correlation filters */
		tom::angleTupel anglesStart;		/* start of angles -phi psi theta */
		tom::angleTupel anglesIncrement; /* increment of angles */
		tom::angleTupel anglesEnd;		/* end of angles */

		std::size_t numberOfParticles;
	
		bool filterSave;
		bool listSave;
		bool stackSave;

	};


}
/** \brief  Sets the weighting of the tree corralation filter components.*/
inline void tom::os3_options::setFilterWeight(double* weight){
	this->weight[0] = weight[0];
	this->weight[1] = weight[1];
	this->weight[2] = weight[2];
}

inline void tom::os3_options::setVolumeDir(std::string str){
	this->volumeDir =str;
}
inline void tom::os3_options::setPatternDir(std::string str){
	this->patternDir =str;
}
inline void tom::os3_options::setMaskDir(std::string str){
	this->maskDir =str;
}
inline void tom::os3_options::setPsfDir(std::string str){
	this->psfDir =str;
}
inline void tom::os3_options::setResultDir(std::string str){
	this->resultDir=str;
}
inline void tom::os3_options::setAngles(tom::angleTupel start,tom::angleTupel inc,tom::angleTupel end){
	this->anglesStart = start;
	this->anglesIncrement = inc;
	this->anglesEnd = end;
}
/** \brief  Returns the 3 values defining the weighting of the tree corralation filter components.*/
inline	double* tom::os3_options::getFilterWeight(){
	return this->weight;
}
inline	std::string tom::os3_options::getVolumeDir(){
	return this->volumeDir;
}
inline	std::string tom::os3_options::getPatternDir(){
	return this->patternDir;
}
inline	std::string tom::os3_options::getMaskDir(){
	return this->maskDir;
}
inline	std::string tom::os3_options::getPsfDir(){
	return this->psfDir;
}
inline	std::string tom::os3_options::getResultDir(){
	return this->resultDir;
}
#endif










