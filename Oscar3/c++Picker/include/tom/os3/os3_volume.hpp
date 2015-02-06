/****************************************************************************//**
 * \file tom_os3_volume.h
 * \brief The header file for the class os3_volume.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    08.01.2008
 *******************************************************************************/ 
#ifndef __INCLUDE_TOM_OS3_VOLUME_H__
#define __INCLUDE_TOM_OS3_VOLUME_H__

#include "tom/os3/os3_types.hpp"
#include "tom/core/volume.hpp"
#include "tom/core/volume_fcn.hpp"
#include "tom/core/fftw_plan.hpp"

#include <iostream>

namespace tom{

/****************************************************************************//**
 * \brief The os3_volume class is used to store the volume and its statistics.
 * 
 * 
 * The os3_volume class used to store the volume and its statistics. \n 
 *******************************************************************************/

	template<typename T>
	class os3_volume{
	
	
		public:
			/*constructors*/		
			os3_volume();
			os3_volume(std::string fileName);
			os3_volume(const std::string fileName,const uint32_t* subregion,const uint32_t* resampling,const uint32_t* binning);
			os3_volume(std::string fileName,std::size_t x,std::size_t y,std::size_t z,std::size_t posX,std::size_t posY,std::size_t posZ);
			os3_volume(tom::os3_volume<T>* original);
			os3_volume(tom::Volume<T>* original);
			~os3_volume();

			void applyMask(tom::os3_volume<T>* maskV);
			void calculateMeanVolume(tom::os3_volume<T>* mask);
			void calculateStdVolume(tom::os3_volume<T>* mask);
			void calculateStatisticVolumes(tom::os3_volume<T>* mask);
			void normaliseVolume(tom::os3_volume<T>* exterStatistics);
			void resetStatistics(tom::os3_volume<T>* mask);
			void resize(tom::os3_volume<T>* biggerVolume);
			void resize(std::size_t x,std::size_t y,std::size_t z,tom::os3_volume<T>* biggerVolume);
			void rotate(double* angles);
			void replaceCurrentPlan(tom::os3_volume<T>* original);
			void setVolume(tom::Volume<T>* original);
			void setVolume(tom::os3_volume<T>* original);
			void transform();
			void transformBack(bool normalise);	
			void writeToEm(std::string fileName);
			
			bool equalSize(const tom::os3_volume<T>* volume);

			double numel();
	
			//some inline methods - code below
			std::size_t getSizeX() const;
			std::size_t getSizeY() const;
			std::size_t getSizeZ() const;

			tom::Volume<T>* getVolume() const;
			tom::Volume<T>* getMeanVolume() const;
			tom::Volume<T>* getStdVolume() const;
			tom::Volume<T>* getTransVolume();
			tom::Volume<std::complex<T> >* getTransComplVolume();
			tom::Volume<std::complex<T> >* getFourierVolume() ; 	

			tom::fftw::Plan<T>* getForwardPlan();
			tom::fftw::Plan<T>* getBackwardPlan();

			const T getCentralMean() const;
			const T getCentralStd() const;
	
			const T getInnerSum();
		
		private:

			/*the image / volume */
			tom::Volume<T>* volume;
		
			/*the fouriertransformed volume*/
			tom::Volume<std::complex<T> >* fourierVolume;
	
			/*statistics of the volume used for flcf*/
			tom::Volume<T>* meanVolume;
			tom::Volume<T>* stdVolume;
			
			/*points to the object containing the mask used for normalisation*/
			tom::os3_volume<T>* mask;
			

			tom::Volume<T>* transVolume; //pointer to memory of same size as volume used for fftw only
			tom::Volume<std::complex<T> >* transComplVolume; //pointer to memory of same size as volume used for fftw only
			tom::fftw::Plan<T>* forwardFtPlan; //the volume specific forward FT plan
			tom::fftw::Plan<T>* backwardFtPlan; //the volume specific backwardFtPlan FT plan

			bool planCalculated;
			void initialise();
			std::string fileName;
			bool recursion;
			bool statisticsAvailable;
	
			T innerSum;
	};
}

 /** \brief Returns the size of the volume along X. */
template<typename T>
inline std::size_t tom::os3_volume<T>::getSizeX() const {
    return this->volume->getSizeX();
}

 /** \brief Returns the size of the volume along Y. */
template<typename T>
inline std::size_t tom::os3_volume<T>::getSizeY() const {
    return this->volume->getSizeY();
}

 /** \brief Returns the size of the volume along Z. */
template<typename T>
inline std::size_t tom::os3_volume<T>::getSizeZ() const {
    return this->volume->getSizeZ();
}

 /** \brief Returns a pointer to the volume data. */
template<typename T>
inline tom::Volume<T>* tom::os3_volume<T>::getVolume() const {
    return this->volume;
}

 /** \brief Returns a pointer to the mean volume data. */
template<typename T>
inline tom::Volume<T>* tom::os3_volume<T>::getMeanVolume() const {
    return this->meanVolume;
}

 /** \brief Returns a pointer to the std volume data. */
template<typename T>
inline tom::Volume<T>* tom::os3_volume<T>::getStdVolume() const {
    return this->stdVolume;
}

 /** \brief Returns a pointer to the std volume data. */
template<typename T>
inline tom::Volume<std::complex<T>  >* tom::os3_volume<T>::getFourierVolume() {
	this->transform();
    return this->fourierVolume;
}

 /** \brief Returns a pointer to the fftw forward plan. */
template<typename T>
inline tom::fftw::Plan<T>* tom::os3_volume<T>::getForwardPlan(){
	return this->forwardFtPlan;
}

 /** \brief Returns a pointer to the fftw backward plan. */
template<typename T>
inline tom::fftw::Plan<T>* tom::os3_volume<T>::getBackwardPlan(){
	return this->backwardFtPlan;
}

 /** \brief Returns a pointer to the transformation volume. */
template<typename T>
inline tom::Volume<T>* tom::os3_volume<T>::getTransVolume(){
	return this->transVolume;
}

 /** \brief Returns a pointer to the complex transformation volume. */
template<typename T>
inline tom::Volume<std::complex<T> >* tom::os3_volume<T>::getTransComplVolume(){
	return this->transComplVolume;
}

/** \brief Returns the central mean value of the volumes meanVolume. */
template<typename T>
inline T const tom::os3_volume<T>::getCentralMean() const {
	std::size_t a = this->getSizeX()/2 ;
	std::size_t b = this->getSizeY()/2 ;
	std::size_t c = this->getSizeZ()/2 ;
	if(this->getSizeZ() == 1 ) c = 0;
	//std::cout << "(" << a << "," << b << "," << c << ")" << std::endl;
    return this->meanVolume->get(a,b,c);
}

/** \brief Returns the central std value of the volumes stdVolume. */
template<typename T>
inline T const tom::os3_volume<T>::getCentralStd() const {
	std::size_t a = this->getSizeX()/2 ;
	std::size_t b = this->getSizeY()/2 ;
	std::size_t c = this->getSizeZ()/2 ;

	if(this->getSizeZ() == 1 ) c = 0;
	//std::cout << "(" << a << "," << b << "," << c << ")" << std::endl;

    return this->stdVolume->get(a,b,c);
}

/** \brief Initialises the volume. */
template<typename T> 
void tom::os3_volume<T>::initialise(){

	//allocate memory for the real space volume where all the transformations take place
	this->transVolume = new tom::Volume<T>(this->volume->getSizeX(), this->volume->getSizeY(), this->volume->getSizeZ(), NULL,NULL);
	//allocate memory, depends on the rank of the volume
	if(volume->getSizeZ() == 1){
		this->fourierVolume = new tom::Volume<std::complex<T> >(this->volume->getSizeX(), this->volume->getSizeY()/2+1,1, NULL,NULL);
		this->transComplVolume = new tom::Volume<std::complex<T> >(this->volume->getSizeX(), this->volume->getSizeY()/2+1,1, NULL,NULL);
	}
	else{
		this->fourierVolume = new tom::Volume<std::complex<T> >(this->volume->getSizeX(), this->volume->getSizeY(), this->volume->getSizeZ()/2+1, NULL,NULL);
		this->transComplVolume = new tom::Volume<std::complex<T> >(this->volume->getSizeX(), this->volume->getSizeY(),this->volume->getSizeZ()/2+1, NULL,NULL);
	}

	//allocate memory for the statistics volumes
	this->meanVolume = new tom::Volume<T>(this->volume->getSizeX(), this->volume->getSizeY(), this->volume->getSizeZ(), NULL,NULL);
	this->stdVolume = new tom::Volume<T>(this->volume->getSizeX(), this->volume->getSizeY(), this->volume->getSizeZ(), NULL,NULL);

	this->mask = NULL;

	this->forwardFtPlan = NULL;
	this->backwardFtPlan= NULL;
	this->planCalculated = false;
	this->recursion = true;
	this->statisticsAvailable = false;
	this->innerSum = 0;
}

/** \brief Initialises the volume. */
template<typename T> 
double tom::os3_volume<T>::numel(){
	return (double) this->getSizeX() *this->getSizeY()*this->getSizeZ();
}

/** \brief Returns the inner sum of the volume. */
template<typename T>
inline T const tom::os3_volume<T>::getInnerSum(){

	
	if(innerSum == (T)0){
		
		for(std::size_t z = 0;z < (T)this->getSizeZ();z++)
			for(std::size_t y = 0;y <(T)this->getSizeY();y++)
				for(std::size_t x = 0;x <(T)this->getSizeX();x++){
					this->innerSum += this->volume->get(x,y,z);
				}
	}

	return innerSum;
	
}

#endif































