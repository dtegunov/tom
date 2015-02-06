 /****************************************************************************//**
 * \file tom_os3_results.cpp
 * \brief The source file for the class os3_results.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    23.01.2008
 *******************************************************************************/ 


#include <iostream>
#include <memory>
#include <sstream>

//os3 include files
#include "tom/os3/os3_results.hpp"
#include "tom/os3/os3_functions.hpp"
#include "tom/os3/os3_volume.hpp"
#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_picklist.hpp"
/*include files from tomc*/
#include "tom/core/fftw_plan.hpp"
#include "tom/core/volume.hpp"
#include "tom/core/volume_fcn.hpp"



/****************************************************************************//**
 * \brief Constructor.
 *
 * Each pointer is set to NULL. 
 *******************************************************************************/
template<typename T> tom::os3_results<T>::os3_results(){

	this->xcfResult = NULL;
	this->psrResult = NULL;	
	this->socResult = NULL;
	this->angResult = NULL;
	
}

/****************************************************************************//**
 * \brief Constructor.
 *
 * Each pointer is initialised. 
 *******************************************************************************/
template<typename T> tom::os3_results<T>::os3_results(std::size_t x,std::size_t y,std::size_t z){

	this->xcfResult = new tom::Volume<T>(x,y,z,NULL,NULL);
	this->psrResult = new tom::Volume<T>(x,y,z,NULL,NULL);	
	this->socResult = new tom::Volume<T>(x,y,z,NULL,NULL);
	this->angResult = new tom::Volume<int>(x,y,z,NULL,NULL);
}

/****************************************************************************//**
 * \brief Constructor.
 *
 * Each pointer is initialised. 
 *******************************************************************************/
template<typename T> tom::os3_results<T>::os3_results(tom::os3_job job){

	tom::read_from_em<T>(this->xcfResult,tom::appendStringToFilename(job.resultPath,"XCF.em").c_str(), NULL,NULL,NULL, NULL, NULL,NULL);
	tom::read_from_em<T>(this->psrResult,tom::appendStringToFilename(job.resultPath,"PSR.em").c_str(), NULL,NULL,NULL, NULL, NULL,NULL);
	tom::read_from_em<T>(this->socResult,tom::appendStringToFilename(job.resultPath,"SOC.em").c_str(), NULL,NULL,NULL, NULL, NULL,NULL);
	tom::read_from_em<int>(this->angResult,tom::appendStringToFilename(job.resultPath,"ANG.em").c_str(), NULL,NULL,NULL, NULL, NULL,NULL);

}


/****************************************************************************//**
 * \brief Destructor.
 *
 * Each volume is freed.
 *******************************************************************************/
template<typename T> tom::os3_results<T>::~os3_results(){

	delete this->xcfResult;
	delete this->psrResult;
	delete this->socResult;
	delete this->angResult;

}


/****************************************************************************//**
 * \brief Sets the values of the xcf result member.
 *
 * Updates the values of the xcf result member.
 *******************************************************************************/
template<typename T> 
void tom::os3_results<T>::update(tom::Volume<T>* xcfVolume,tom::Volume<T>* psrVolume,tom::Volume<T>* socVolume,int angleIndex){

	if(this->xcfResult == NULL){
		std::auto_ptr<tom::Volume<T> > autoXCF(this->xcfResult); 
		this->xcfResult = new tom::Volume<T>(xcfVolume->getSizeX(),xcfVolume->getSizeY(),xcfVolume->getSizeZ(),NULL,NULL);
		this->xcfResult->setValues(*xcfVolume);

		std::auto_ptr<tom::Volume<T> > autoPSR(this->psrResult); 
		this->psrResult = new tom::Volume<T>(xcfVolume->getSizeX(),xcfVolume->getSizeY(),xcfVolume->getSizeZ(),NULL,NULL);
		this->psrResult->setValues(*psrVolume);

		std::auto_ptr<tom::Volume<T> > autoSOC(this->socResult); 
		this->socResult = new tom::Volume<T>(xcfVolume->getSizeX(),xcfVolume->getSizeY(),xcfVolume->getSizeZ(),NULL,NULL);
		this->socResult->setValues(*socVolume);

		std::auto_ptr<tom::Volume<int> > autoANG(this->angResult); 
		this->angResult = new tom::Volume<int>(xcfVolume->getSizeX(),xcfVolume->getSizeY(),xcfVolume->getSizeZ(),NULL,NULL);
		this->angResult->setValues(angleIndex);

		autoXCF.release();
		autoPSR.release();
		autoSOC.release();
		autoANG.release();
	}
	else{
		
		T* xcfOld = NULL;
		T* xcfNew = NULL;
		T* psrOld = NULL;
		T* psrNew = NULL;
		T* socOld = NULL;
		T* socNew = NULL;
		
		bool xcfMax;
		bool psrMax;
		bool socMax;
		for(std::size_t z = 0;z < (T)this->xcfResult->getSizeZ();z++)
			for(std::size_t y = 0;y <(T)this->xcfResult->getSizeY();y++)
				for(std::size_t x = 0;x <(T)this->xcfResult->getSizeX();x++){
					
					xcfOld = &this->xcfResult->get(x,y,z);
					xcfNew = &xcfVolume->get(x,y,z);

					psrOld = &this->psrResult->get(x,y,z);
					psrNew = &psrVolume->get(x,y,z);

					socOld = &this->socResult->get(x,y,z);
					socNew = &socVolume->get(x,y,z);

					xcfMax = *xcfOld < *xcfNew;
					psrMax = *psrOld < *psrNew;
					socMax = *socOld < *socNew;

					if(xcfMax && psrMax && socMax){
						*xcfOld = *xcfNew;
						*psrOld = *psrNew;
						*socOld = *socNew;
						this->angResult->get(x,y,z) = angleIndex;
					}
		}	
	}
}

/****************************************************************************//**
 * \brief Writes the results to file.
 *
 * Writes the results to file.
 *******************************************************************************/
template<typename T> 
inline void tom::os3_results<T>::writeToEm(std::string filename){
	
	this->xcfResult->write_to_em(tom::appendStringToFilename(filename,"XCF.em").c_str(),NULL);
	
	this->psrResult->write_to_em(tom::appendStringToFilename(filename,"PSR.em").c_str(),NULL);
	
	this->socResult->write_to_em(tom::appendStringToFilename(filename,"SOC.em").c_str(),NULL);
	
	this->angResult->write_to_em(tom::appendStringToFilename(filename,"ANG.em").c_str(),NULL);
}

/****************************************************************************//**
 * \brief Writes the results to file.
 *
 * Writes the results to file.
 *******************************************************************************/
template<typename T> 
void tom::os3_results<T>::setPicklist(tom::os3_job& job,bool printStatus){

	if(printStatus)
			jobInfo(job);

		//read the header of the tempalte
		tom_io_em_header patternHeader;
		
		if(tom_io_em_read_header(job.patternPath.c_str(), NULL, &patternHeader,NULL))
			throw std::runtime_error("\n tom::os3_picker::createParticleStack() - Error reading pattern file.");

		if(printStatus)
			std::cout << "createParticleStack : results read." << std::endl;
		//sum the results up a*xcf+b*psr+c*soc

		tom::Volume<T>* xcf;
		std::auto_ptr<tom::Volume<T> > xcfP(xcf);
		
		//copy, keep the orignal results
		xcf= new tom::Volume<T>(*this->xcfResult);
		xcf->shift_scale(0.,job.weight[0]);	
		
		tom::Volume<T>* psr;
		std::auto_ptr<tom::Volume<T> > psrP(psr);
		psr= new tom::Volume<T>(*this->psrResult);
		psr->shift_scale(0.,job.weight[1]);	
		
		tom::Volume<T>* soc;
		std::auto_ptr<tom::Volume<T> > socP(soc);
		soc= new tom::Volume<T>(*this->socResult);
		soc->shift_scale(0.,job.weight[2]);	
		
		tom::element_wise_add(*xcf,*psr);
		tom::element_wise_add(*xcf,*soc);
	
		if(printStatus)
			std::cout << "createParticleStack : weighting applied and sum created." << std::endl;

		bool twoDimensions = patternHeader.dims[2] <= 1;		
		std::size_t offsetX = patternHeader.dims[0]/2+1;
		std::size_t offsetY = patternHeader.dims[1]/2+1;
		std::size_t offsetZ;
		if(twoDimensions)
			offsetZ = 1;
		else
			offsetZ = patternHeader.dims[2]/2+1;


		//set values at the volume boundaries to -1000
		tom::Volume<T> zeros(*xcf);
		zeros.setValues(0);
		tom::Volume<T>* zerosWindow;
		std::auto_ptr<tom::Volume<T> > zerosPrt(zerosWindow);

		if(twoDimensions)
			zerosWindow = new tom::Volume<T>(zeros,  &zeros.get(offsetX,offsetY,0),
											zeros.getSizeX()-(offsetX*2),zeros.getSizeY()-(offsetY*2),1, 
											zeros.getStrideX(), zeros.getStrideY(), zeros.getStrideZ());
		else
			zerosWindow = new tom::Volume<T>(zeros,  &zeros.get(offsetX,offsetY,offsetZ),
											zeros.getSizeX()-(offsetX*2),zeros.getSizeY()-(offsetY*2),zeros.getSizeZ()-(offsetZ*2), 
											zeros.getStrideX(), zeros.getStrideY(), zeros.getStrideZ());
		
		zerosWindow->setValues(-1000);
		tom::element_wise_sub(*xcf,zeros);

		
		//create a window to the peaks - erase around the edges
		tom::Volume<T>* xcfWindow = NULL;
		std::auto_ptr<tom::Volume<T> > windowPtr(xcfWindow);

		if(twoDimensions)
			xcfWindow = new tom::Volume<T>(*xcf,  &xcf->get(offsetX,offsetY,0),
									xcf->getSizeX()-(offsetX*2),xcf->getSizeY()-(offsetY*2),1, 
									xcf->getStrideX(), xcf->getStrideY(), xcf->getStrideZ());
		else
			xcfWindow = new tom::Volume<T>(*xcf,  &xcf->get(offsetX,offsetY,offsetZ),
									xcf->getSizeX()-(offsetX*2),xcf->getSizeY()-(offsetY*2),xcf->getSizeZ()-(offsetZ*2), 
									xcf->getStrideX(), xcf->getStrideY(), xcf->getStrideZ());
		
		tom::os3_pick<T> pick;

		pick.volumePath = job.volumePath;
		pick.patternPath = job.patternPath;
		//get the first n particles
		std::size_t counter=0;
		
		while(counter < job.numberParticles){
			
			if(printStatus)
				std::cout << "Peak no: " << counter << " " << std::endl;
			
			std::vector<st_idx> peaks = tom::peak(*xcfWindow);
			st_idx peak = *(peaks.begin());

			//delete pattern area			
			tom::Volume<T> xcfSmallWindow(*xcf,&xcf->get(peak.x,peak.y,peak.z),
										  patternHeader.dims[0],patternHeader.dims[1],patternHeader.dims[2],
										  xcf->getStrideX(), xcf->getStrideY(), xcf->getStrideZ());
			
			xcfSmallWindow.setValues(-1000);
			
			int angleIndex = this->angResult->get(peak.x,peak.y,peak.z);
			pick.angles = job.angleList[angleIndex];
			
			peak.x += offsetX + job.subregion[0];
			peak.y += offsetY + job.subregion[1];

			(twoDimensions)?peak.z = 1 : peak.z = offsetZ + job.subregion[2];

			pick.position = peak;
			this->picklist.push_back(pick);
			counter++;
		}	
}

/****************************************************************************//**
 * \brief Writes the results to file.
 *
 * Writes the results to file.
 *******************************************************************************/
template<typename T> 
void tom::os3_results<T>::writeToEm(tom::os3_job job){
	
	if(job.numberJobs >1){
		
		std::string filename = job.resultPath;

		const uint32_t* subregion = job.subregion;
		const uint32_t* writeBackSubregion = job.writeBackSubregion;
		std::size_t posX =(std::size_t) writeBackSubregion[0];
		std::size_t posY =(std::size_t) writeBackSubregion[1];
		std::size_t posZ =(std::size_t) writeBackSubregion[2];
		
		std::size_t numelX=(std::size_t) writeBackSubregion[3];
		std::size_t numelY=(std::size_t) writeBackSubregion[4];
		std::size_t numelZ=(std::size_t) writeBackSubregion[5];
		
		uint32_t pos[3] = {(uint32_t)subregion[0]+(uint32_t)posX,
						   (uint32_t)subregion[1]+(uint32_t)posY,
						   (uint32_t)subregion[2]+(uint32_t)posZ};

		tom::Volume<T> xcfWindow(*this->xcfResult,  &this->xcfResult->get(posX,posY,posZ),
								numelX,numelY,numelZ, 
								this->xcfResult->getStrideX(), this->xcfResult->getStrideY(), this->xcfResult->getStrideZ());
		xcfWindow.write_to_em(tom::appendStringToFilename(filename,"XCF.em").c_str(),NULL,&pos[0]);
		
		tom::Volume<T> psrWindow(*this->psrResult,  &this->psrResult->get(posX,posY,posZ),
								numelX,numelY,numelZ, 
								this->psrResult->getStrideX(), this->psrResult->getStrideY(), this->psrResult->getStrideZ());
		psrWindow.write_to_em(tom::appendStringToFilename(filename,"PSR.em").c_str(),NULL,&pos[0]);

		tom::Volume<T> socWindow(*this->socResult,  &this->socResult->get(posX,posY,posZ),
								numelX,numelY,numelZ, 
								this->socResult->getStrideX(), this->socResult->getStrideY(), this->socResult->getStrideZ());
		socWindow.write_to_em(tom::appendStringToFilename(filename,"SOC.em").c_str(),NULL,&pos[0]);

		tom::Volume<int> angWindow(*this->angResult,  &this->angResult->get(posX,posY,posZ),
									numelX,numelY,numelZ, 
									this->angResult->getStrideX(), this->angResult->getStrideY(), this->angResult->getStrideZ());
		angWindow.write_to_em(tom::appendStringToFilename(filename,"ANG.em").c_str(),NULL,&pos[0]);
	}else
		this->writeToEm(job.resultPath);
}








template class tom::os3_results<TFLOAT>;
template class tom::os3_results<TDOUBLE>;




































