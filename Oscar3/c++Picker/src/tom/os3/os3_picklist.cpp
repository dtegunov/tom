#include "tom/os3/os3_picklist.hpp"
/*----------*/
#include "tom/core/volume.hpp"
/*----------*/
#include "boost/lexical_cast.hpp"
#include "helper/filesystem.hpp"
/*----------*/
#include <iostream>
#include <fstream>

template<typename T>
tom::os3_picklist<T>::os3_picklist(){
	this->particleStack = NULL;
}

template<typename T>
tom::os3_picklist<T>::os3_picklist(const os3_picklist& original){
	this->particleStack = NULL;
	this->picklist = original.picklist;
}

template<typename T>
tom::os3_picklist<T>::os3_picklist(const std::string path){
	this->particleStack = NULL;
}



/****************************************************************************//**
 * \brief Creates a particle stack from the picklist stored.
 *
 * The picklist may contain picks from different volumes and different templates, \n but the templates must have the same size!
 *******************************************************************************/
template<typename T>
void tom::os3_picklist<T>::createParticleStack(){
	
	tom::Volume<T>* volume = NULL;
	std::auto_ptr<tom::Volume<T> > autoVolume(volume); 
	
	//read the header of the tempalte
	tom_io_em_header patternHeader;		
	if(tom_io_em_read_header((*this->picklist.begin()).patternPath.c_str(), NULL, &patternHeader,NULL))
		throw std::runtime_error("\n tom::os3_picker::createParticleStack() - Error reading pattern file.");

	//create the particlestack in mem
	std::auto_ptr<tom::Volume<T> > autoStack(this->particleStack);	
	this->particleStack = new tom::Volume<T>(patternHeader.dims[0],patternHeader.dims[1],patternHeader.dims[2] * this->picklist.size(),NULL,NULL);
	
	typename std::vector<tom::os3_pick<T> >::iterator iterator;
	std::size_t counter = 0;
	std::string volumeName = "NONE";

	//for each entry in the picklist, select the volume, locate the particle, cut it out, rotate and place it into the stack
	for(iterator = this->picklist.begin();iterator != this->picklist.end();iterator++){
		
		if(volumeName.compare(std::string("NONE")) == 0 || !(*iterator).volumePath.compare(volumeName) == 0){
			volumeName = (*iterator).volumePath;
			//read the volume
			tom::read_from_em<T>(volume, volumeName.c_str(), NULL,NULL,NULL, NULL, NULL,NULL);
		}
		
		tom::Volume<T> volumeWindow(*volume,  &volume->get((*iterator).position.x,(*iterator).position.y,(*iterator).position.z),
									patternHeader.dims[0],patternHeader.dims[1],patternHeader.dims[2],
									volume->getStrideX(), volume->getStrideY(), volume->getStrideZ());

		//rotate particle
		tom::Volume<T> particle(volumeWindow.getSizeX(),volumeWindow.getSizeY(),volumeWindow.getSizeZ(),NULL,NULL);
		tom::rotate(volumeWindow,particle,-(*iterator).angles.phi,-(*iterator).angles.theta,-(*iterator).angles.psi);
		

		tom::Volume<T> stackWindow(*this->particleStack,&this->particleStack->get(0,0,patternHeader.dims[2]*(counter)),
								patternHeader.dims[0],patternHeader.dims[1],patternHeader.dims[2],
								particleStack->getStrideX(),particleStack->getStrideY(),particleStack->getStrideZ());
				
		//copy values into the particleStack
		stackWindow.setValues(particle);
		counter++;
	}
	autoStack.release();
}


/****************************************************************************//**
 * \brief Saves the picklist to path.
 *
 * \param[in] path Path to picklist file.
 *
 * Saves the picklist to ascii file. 
 *******************************************************************************/
template<typename T>
void tom::os3_picklist<T>::savePicklist(std::string path){

	std::ofstream fileOut;

	fileOut.open(path.c_str());
	
	if(fileOut.bad())
		throw std::runtime_error("\n tom::os3_picker::os3_savePicklist() - Error saving picklist! Wrong path?");
		
	typename std::vector<tom::os3_pick<T> >::iterator iterator;
					
	for(iterator= this->picklist.begin();iterator != this->picklist.end();iterator++){
				
		std::string x,y,z,phi,psi,theta,value;
				
		x = boost::lexical_cast<std::string>((*iterator).position.x);
		y = boost::lexical_cast<std::string>((*iterator).position.y);
		z = boost::lexical_cast<std::string>((*iterator).position.z);
							
		phi = boost::lexical_cast<std::string>((*iterator).angles.phi);
		psi = boost::lexical_cast<std::string>((*iterator).angles.psi);
		theta = boost::lexical_cast<std::string>((*iterator).angles.theta);
		
		value = boost::lexical_cast<std::string>((*iterator).value);

		std::string outString = (*iterator).volumePath + " | " + (*iterator).patternPath + " | " 
								+ x + " | " + y + " | " + z + " | " 
								+ phi+ " | " + psi+ " | " + theta + " | " 
								+ value; 
		fileOut << outString << std::endl;
	}
	fileOut.close();
}

/** \brief Helper function.*/
namespace{
	template <typename T>
	inline tom::os3_pick<T> parsePicklistLine(std::string line,std::string splitChar){
	
		tom::os3_pick<T> pick;
		// /fs/home/hrabe/develop/images/20S/volumes/20S_3_2bin.em | /fs/home/hrabe/develop/images/20S/pattern/20S_2bin.em | 904 | 577 | 0 | 5.759586531581287 | 0 | 0 | 0.57
		std::vector<std::string> stringV;
		std::string::size_type loc;
		while( (loc = line.find_first_of(splitChar)) != std::string::npos){
			
			stringV.push_back(line.substr(0,loc));
			line = line.substr(loc+1);
		}
		stringV.push_back(line);

		pick.volumePath = stringV[0];
		pick.patternPath = stringV[1];
		
		
		pick.position.x = boost::lexical_cast<std::size_t> (stringV[2].substr(1,stringV[2].size()-2));
		pick.position.y = boost::lexical_cast<std::size_t>(stringV[3].substr(1,stringV[3].size()-2));
		pick.position.z = boost::lexical_cast<std::size_t>(stringV[4].substr(1,stringV[4].size()-2));
	
		pick.angles.phi = boost::lexical_cast<double>(stringV[5].substr(1,stringV[5].size()-2));
		pick.angles.psi = boost::lexical_cast<double>(stringV[6].substr(1,stringV[6].size()-2));
		pick.angles.theta = boost::lexical_cast<double>(stringV[7].substr(1,stringV[7].size()-2));

		pick.value = boost::lexical_cast<T>(stringV[8].substr(1,stringV[8].size()-2));
	
		return pick;
	}
}

/****************************************************************************//**
 * \brief Loads a picklist.
 *
 * \param[in] path Path to picklist file.
 *
 * Loads the picklist from ascii file. 
 *******************************************************************************/
template<typename T>
void tom::os3_picklist<T>::loadPicklist(std::string path,bool clear){

	std::ifstream fileIn;

	fileIn.open(path.c_str());

	if(fileIn.bad())
		throw std::runtime_error("\n tom::os3_picker::os3_savePicklist() - Error saving picklist! Wrong path?");


	if(clear)
		this->picklist.clear();

	std::string line;
	std::getline(fileIn,line);
	while(! fileIn.eof()){
		
		//split line, parse values into this->picklist
		tom::os3_pick<T> pick= parsePicklistLine<T>(line,std::string("|"));
		
		this->picklist.push_back(pick);
		std::getline(fileIn,line);	
	}
	
	fileIn.close();
}

template class tom::os3_picklist<TFLOAT>;
template class tom::os3_picklist<TDOUBLE>;






