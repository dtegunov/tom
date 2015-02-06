 /****************************************************************************//**
 * \file tom_os3_volume.cpp
 * \brief The source file for the class os3_volume.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    09.01.2008
 *******************************************************************************/ 

#include <iostream>
#include <sstream>

//os3 include files
#include "tom/os3/os3_volume.hpp"
#include "tom/os3/os3_types.hpp"
/*include files from tomc*/
#include "tom/core/fftw_plan.hpp"
#include "tom/core/volume.hpp"
#include "tom/core/volume_fcn.hpp"

/****************************************************************************//**
 * \brief Constructor of the os3_volume class.
 *
 * Each pointer is set to NULL. 
 *******************************************************************************/
template<typename T> tom::os3_volume<T>::os3_volume(){

		this->volume = NULL;
		this->fourierVolume = NULL;
		
		/*statistics of the volume used for flcf*/
		this->meanVolume = NULL;
		this->stdVolume = NULL;

		/* points to the object containing the mask used for normalisation */
		this->mask = NULL;
	
		/*memory allocated for fourier transformation*/
		this->transVolume = NULL;
		this->transComplVolume= NULL;

		/*plans*/
		this->forwardFtPlan = NULL;
		this->backwardFtPlan = NULL;
		this->recursion = false;
	}


template<typename T> tom::os3_volume<T>::os3_volume(tom::Volume<T>* original){

		this->volume = new tom::Volume<T>(*original);

		this->initialise();
	}

/****************************************************************************//**
 * \brief Constructor of the os3_volume class.
 *
 * 
 *******************************************************************************/
template<typename T> tom::os3_volume<T>::os3_volume(std::string fileName){

	try {
		//read the volume from file
		tom::read_from_em<T>(this->volume, fileName.c_str(), NULL,NULL,NULL, NULL, NULL,NULL);
	}catch(std::exception &e){
		std::cerr << e.what();
	}
	this->fileName = fileName;
	this->initialise();
}
/****************************************************************************//**
 * \brief Constructor of the os3_volume class.
 *
 * 
 *******************************************************************************/
template<typename T> tom::os3_volume<T>::os3_volume(const std::string fileName,const uint32_t* subregion,const uint32_t* resampling,const uint32_t* binning){
	try {
		//read the volume from file
		tom::read_from_em<T>(this->volume, fileName.c_str(), subregion,resampling,binning, NULL, NULL,NULL);
		//std::cout << subregion[0] << "," << subregion[1]<< "," << subregion[2]<< "," << subregion[3] << "," << subregion[4]<< "," << subregion[5]<< std::endl;
	}catch(std::exception &e){
		std::cout << e.what();
	}
	this->fileName = fileName;
	this->initialise();

}
/****************************************************************************//**
 * \brief A constructor of the os3_volume class.
 *
 * \param[in] fileName The absolute path to the em file.
 * \param[in] x X size of the bigger volume.
 * \param[in] y Y size of the bigger volume.
 * \param[in] z Z size of the bigger volume.
 * \param[in] posX The X position of the first voxel which is beeing pasted.
 * \param[in] posY The Y position of the first voxel which is beeing pasted.
 * \param[in] posZ The Z position of the first voxel which is beeing pasted.
 *
 * Reads the given em file to memory and pastes the volume read into another volume of size x,y,z.
 *******************************************************************************/
template<typename T> tom::os3_volume<T>::os3_volume(std::string fileName,std::size_t x,std::size_t y,std::size_t z,std::size_t posX,std::size_t posY,std::size_t posZ){

	tom::Volume<T>* smallerVolume;
	try {
		//read the volume from file
		tom::read_from_em<T>(smallerVolume, fileName.c_str(), NULL,NULL,NULL, NULL, NULL,NULL);
		
	}catch(std::exception &e){
		std::cout << e.what();
	}

	this->volume = new tom::Volume<T>(x,y,z,NULL,NULL);
	this->volume->setValues(0);
	tom::Volume<T> windowVolume(*this->volume,  &this->volume->get(posX,posY,posZ),
								smallerVolume->getSizeX(),smallerVolume->getSizeY(),smallerVolume->getSizeZ(), 
								this->volume->getStrideX(), this->volume->getStrideY(), this->volume->getStrideZ());

	windowVolume.setValues(*smallerVolume);
	this->fileName = fileName;
	this->initialise();
}

/****************************************************************************//**
 * \brief Copy constructor of the class.
 *
 * Copy constructor of the class.
 *******************************************************************************/
template<typename T> tom::os3_volume<T>::os3_volume(tom::os3_volume<T>* original){
		
		this->volume = new tom::Volume<T>(*original->getVolume());
		this->fourierVolume = new tom::Volume<std::complex<T> >(*original->getFourierVolume());
		
		this->meanVolume = new tom::Volume<T>(*original->getMeanVolume());
		this->stdVolume= new tom::Volume<T>(*original->getStdVolume());
		this->statisticsAvailable = true;
			
		this->transVolume = NULL;
		this->transComplVolume = NULL;
		this->forwardFtPlan = NULL;
		this->backwardFtPlan = NULL;

		if(original->getForwardPlan() != NULL){
			this->replaceCurrentPlan(original);	
			this->recursion = true;
			this->fileName = original->fileName;
		}	
	}


/****************************************************************************//**
 * \brief Destructor of the os3_volume class.
 *
 * Frees all the memory used by the object.
 *******************************************************************************/
template<typename T> tom::os3_volume<T>::~os3_volume(){
		
		//std::cout << "Deleting volume : " << this->fileName << std::endl;

		delete this->volume;
		delete this->fourierVolume;
		delete this->meanVolume;
		delete this->stdVolume;
		delete this->transVolume;
		delete this->transComplVolume;
		delete this->forwardFtPlan;		
		delete this->backwardFtPlan;
		
	}

/****************************************************************************//**
 * \brief Replaces the volume stored to the one determined by extVolume
 *
 * \param[in] extVolume 
 *
 * Replaces the volume stored to the one determined by volume.
 *******************************************************************************/
template<typename T> 
void tom::os3_volume<T>::setVolume(tom::Volume<T>* extVolume){
		/*clear the old values from memory*/
		delete this->volume;
		delete this->fourierVolume;
		delete this->meanVolume;
		delete this->stdVolume;
		delete this->transVolume;
		delete this->transComplVolume;
		delete this->forwardFtPlan;
		delete this->backwardFtPlan;

		this->volume = NULL;
		this->fourierVolume = NULL;
		this->meanVolume = NULL;
		this->stdVolume = NULL;
		this->transVolume = NULL;
		this->transComplVolume =NULL;
		this->forwardFtPlan = NULL;
		this->backwardFtPlan = NULL;

		this->volume = new tom::Volume<T> (extVolume->getSizeX(), extVolume->getSizeY(), extVolume->getSizeZ(), NULL,NULL);
		this->volume->setValues(*extVolume);
		this->initialise();
}

/****************************************************************************//**
 * \brief Replaces the plan and the transform volumes by the ones of the original volume.
 *
 * \param[in] original The original storing the new plan.
 *
 * Replaces the plan and the transform volumes by the ones of the original volume.
 *******************************************************************************/
template<typename T> 
inline void tom::os3_volume<T>::replaceCurrentPlan(tom::os3_volume<T>* original){

	delete this->transVolume;
	delete this->transComplVolume;
	delete this->forwardFtPlan;
	delete this->backwardFtPlan;

	std::auto_ptr<tom::Volume<T> > transV(transVolume);
	std::auto_ptr<tom::Volume<std::complex<T> > > transCV(transComplVolume);
	
	this->transVolume = new tom::Volume<T>(*original->getTransVolume(),&(original->getTransVolume()->get()),
											   original->getTransVolume()->getSizeX(),original->getTransVolume()->getSizeY(),original->getTransVolume()->getSizeZ(),
											   original->getTransVolume()->getStrideX(),original->getTransVolume()->getStrideY(),original->getTransVolume()->getStrideZ());
	
	this->transComplVolume = new tom::Volume<std::complex<T> >(*original->getTransComplVolume(),&(original->getTransComplVolume()->get()),
											   				    original->getTransComplVolume()->getSizeX(),original->getTransComplVolume()->getSizeY(),original->getTransComplVolume()->getSizeZ(),
											   				    original->getTransComplVolume()->getStrideX(),original->getTransComplVolume()->getStrideY(),original->getTransComplVolume()->getStrideZ());
	
	this->forwardFtPlan = new tom::fftw::Plan<T>(*original->getForwardPlan());
	this->backwardFtPlan= new tom::fftw::Plan<T>(*original->getBackwardPlan());

	this->planCalculated = true;

	transV.release();
	transCV.release();

}

/****************************************************************************//**
 * \brief Replaces the volume stored to the one determined by extVolume
 *
 * \param[in] extVolume 
 *
 * Replaces the volume stored to the one determined by volume.
 *******************************************************************************/
template<typename T> 
void tom::os3_volume<T>::setVolume(tom::os3_volume<T>* original){
		/*clear the old values from memory*/

		delete this->volume;
		delete this->fourierVolume;
		delete this->meanVolume;
		delete this->stdVolume;
		delete this->transVolume;
		delete this->transComplVolume;
		delete this->forwardFtPlan;
		delete this->backwardFtPlan;
	
		this->volume = NULL;
		this->fourierVolume = NULL;
		this->meanVolume = NULL;
		this->stdVolume = NULL;
		this->transVolume = NULL;
		this->transComplVolume =NULL;
		this->forwardFtPlan = NULL;
		this->backwardFtPlan = NULL;
	
		
		
		this->volume = new tom::Volume<T>(*original->getVolume());
		this->fourierVolume = new tom::Volume<std::complex<T> >(*original->getFourierVolume());
		
		this->meanVolume = new tom::Volume<T>(*original->getMeanVolume());
		this->stdVolume= new tom::Volume<T>(*original->getStdVolume());
		this->statisticsAvailable = true;

		this->replaceCurrentPlan(original);

		this->recursion = true;
		this->statisticsAvailable = false;
		this->innerSum = 0;
}

/****************************************************************************//**
 * \brief Resizes the volume stored.
 *
 * \param[in] x Size in x. 
 * \param[in] y Size in y. 
 * \param[in] z Size in z. 
 *
 * Resizes the volume stored. The original volume is pasted into the center of the new volume.
 *******************************************************************************/
template<typename T> 
void tom::os3_volume<T>::resize(tom::os3_volume<T>* biggerVolume){
	
	std::size_t x = biggerVolume->getSizeX();
	std::size_t y = biggerVolume->getSizeY();
	std::size_t z = biggerVolume->getSizeZ();
	
	this->resize(x,y,z,biggerVolume);

}
/****************************************************************************//**
 * \brief Resizes the volume stored.
 *
 * \param[in] x Size in x. 
 * \param[in] y Size in y. 
 * \param[in] z Size in z. 
 *
 * Resizes the volume stored. The original volume is pasted into the center of the new volume.
 *******************************************************************************/
template<typename T> 
void tom::os3_volume<T>::resize(std::size_t x,std::size_t y,std::size_t z,tom::os3_volume<T>* biggerVolume){
		/*clear the old values from memory*/

		std::size_t oldSizeX = this->getSizeX();
		std::size_t oldSizeY = this->getSizeY();
		std::size_t oldSizeZ = this->getSizeZ();

		delete this->fourierVolume;
		delete this->meanVolume;
		delete this->stdVolume;
		delete this->transVolume;
		delete this->transComplVolume;
		delete this->forwardFtPlan;
		delete this->backwardFtPlan;

		this->fourierVolume = NULL;
		this->meanVolume = NULL;
		this->stdVolume = NULL;
		this->transVolume = NULL;
		this->transComplVolume =NULL;
		this->forwardFtPlan = NULL;
		this->backwardFtPlan = NULL;

		std::size_t posX = x/2 - oldSizeX/2;
		std::size_t posY = y/2 - oldSizeY/2;
		std::size_t posZ = z/2 - oldSizeZ/2;
		
		if(z <= 1)
			posZ = 0;

		tom::Volume<T>* newVolume = new tom::Volume<T> (x,y,z, NULL,NULL);
		newVolume->setValues(0);

		tom::Volume<T> windowVolume(*newVolume,  &newVolume->get(posX,posY,posZ),
									oldSizeX,oldSizeY,oldSizeZ, 
									newVolume->getStrideX(), newVolume->getStrideY(), newVolume->getStrideZ());
		
		windowVolume.setValues(*this->volume);
		
		delete this->volume;
		this->volume = newVolume;
	
		//if the plan for the new size of the image is available, copy the plan from the bigger volume, 
		//generate it othervise
		if(biggerVolume != NULL && biggerVolume->getForwardPlan() != NULL){
			this->replaceCurrentPlan(biggerVolume);
			this->recursion = true;
			this->statisticsAvailable = false;
			this->innerSum = 0;
		}
		else
			this->initialise();

		
}


/****************************************************************************//**
 * \brief Rotates the volume stored.
 *
 * Rotates the volume stored.
 *******************************************************************************/
template<typename T> 
void tom::os3_volume<T>::rotate(double* angles){
	tom::Volume<T>* rotVol =new Volume<T>(this->getSizeX(),this->getSizeY(),this->getSizeZ(),NULL,NULL);
	tom::rotate(*this->volume,*rotVol,angles[0],angles[1],angles[2]);

	delete this->volume;
	this->volume = rotVol;

	this->transform();
	
//		throw std::runtime_error("tom::os3_volume::rotate - error in this->transform");
}


/****************************************************************************//**
 * \brief Applies a mask to the current volume.
 *
 * Applies a mask to the current volume.
 *******************************************************************************/
template<typename T> void tom::os3_volume<T>::applyMask(tom::os3_volume<T>* maskV){

		tom::Volume<T> mask(*(maskV->getVolume()));
		tom::element_wise_multiply(*this->volume,mask);
		
		this->planCalculated = true;
		this->recursion = true;
		this->statisticsAvailable = false;
		this->innerSum = 0;
		this->transform();
}

/****************************************************************************//**
 * \brief Calculates the statistics under the provided mask.
 *
 * \param[in] mask The mask used. 
 * \param[in] numberOfMaskVoxels The number of set voxels in the mask.
 * 
 * Calculates the statistics under the provided mask.
 * This is the C++ implementation of the matlab functions tom_os3_mean and tom_os3_std. 
 * These functions are a part of the fast local normalisation algorithm introduced by Roseman (2003).
 * See the matlab code and the Roseman paper for a deeper understanding of the algorithm.
 *******************************************************************************/
template<typename T> void tom::os3_volume<T>::calculateStatisticVolumes(tom::os3_volume<T>* mask){

		//skip if the statistics have been calculated already 
		if(this->statisticsAvailable)
			return;
		
		if(! this->equalSize(mask))
			throw std::runtime_error("tom::os3_volume<T>::calculateStatisticVolumes : The size of the volume and the mask do not match.\n");
		
		this->calculateMeanVolume(mask);
		
		this->calculateStdVolume(mask);

		//check deviation for values
		if(this->stdVolume->mean() <= 0.01 && this->recursion){
			//shift and scale the volume
// 			std::cout << std::endl<< "calculateStatisticVolumes : The deviation of the volume is very low. STD = " <<  stdVolume->mean() << ". The volume will be manipulated to ensure an accurate correlation." << std::endl;
			double m;
			double v;
			this->volume->stat(m,v,false);
			//std::cout << m << " " << v << std::endl;
			this->volume->shift_scale((T)-m,(T)(1.0/sqrt(v)));
			this->recursion = false;
			this->calculateStatisticVolumes(mask);
		}
		//save pointer to the mask
		this->mask = mask;
		this->statisticsAvailable = true;
	}

/****************************************************************************//**
 * \brief Resets the statistics of the volume.
 *
 * Resets the statistics of the volume. Must be called after update.
 *******************************************************************************/
template<typename T> void tom::os3_volume<T>::resetStatistics(tom::os3_volume<T>* mask){
	this->statisticsAvailable = false;
	this->recursion = true;
	delete this->meanVolume;
	meanVolume = NULL;
	delete this->stdVolume;
	stdVolume = NULL;

	if(!this->equalSize(mask)){
			tom::os3_volume<T> biggerMask(mask);
			biggerMask.resize(this->getSizeX(),this->getSizeY(),this->getSizeZ(),this);		
			this->calculateStatisticVolumes(&biggerMask);
	}
	else
		this->calculateStatisticVolumes(mask);

}

/****************************************************************************//**
 * \brief Calculates the standard deviation volume.
 *
 * 
 * Calculates the standard deviation volume.
 *******************************************************************************/
template<typename T>
inline void tom::os3_volume<T>::calculateStdVolume(tom::os3_volume<T>* mask){
	

	tom::Volume<std::complex<T> > maskF(*mask->getFourierVolume());

	/*square the volume for calculation*/
	tom::Volume<T> volumeSquared(*this->volume);
	tom::element_wise_multiply(volumeSquared,volumeSquared);
	
	/*transform the squared volume to fourierspace*/
	this->transVolume->setValues(volumeSquared);
	this->forwardFtPlan->execute(*this->transVolume,*this->transComplVolume);
	
	//convolve both volumes
	tom::element_wise_multiply(*this->transComplVolume,maskF);

	//use the plan of the volume to transform the volume back to real space		
	//process transformation to real space
	this->backwardFtPlan->execute(*this->transComplVolume,*this->transVolume);

	if(this->stdVolume == NULL)	
			this->stdVolume = new tom::Volume<T>(mask->getSizeX(), mask->getSizeY(), mask->getSizeZ(), NULL,NULL);

	
	//copy the result into the destination and shift it
	tom::fftshift(*transVolume,*this->stdVolume,true);

	//divide by the number of elements
	this->stdVolume->shift_scale(0.,1./((((T) this->numel()) * mask->getInnerSum())));

	tom::Volume<T> meanVolumeSquared(*this->meanVolume);
	tom::element_wise_multiply(meanVolumeSquared,meanVolumeSquared);

	//element wise sub stdVolume and meanVolume
	tom::element_wise_sub(*this->stdVolume,meanVolumeSquared);

	//replace each voxel smaller than 0 in stdVolume by 0
	tom::element_wise_max(*this->stdVolume,(T) 0.0);
		
	//element wise square root
	tom::element_wise_power(*this->stdVolume,(T) 0.5);
}

/****************************************************************************//**
 * \brief Calculates the mean volume.
 *
 * 
 * Calculates the mean volume.
 *******************************************************************************/
template<typename T>
inline void tom::os3_volume<T>::calculateMeanVolume(tom::os3_volume<T>* mask){
	
	tom::Volume<std::complex<T> > maskF(*mask->getFourierVolume());
	
	//copy the fourierVolume to the transformation volume	
	this->transComplVolume->setValues(*this->getFourierVolume());
	
	//convolve both volumes
	tom::element_wise_multiply(*this->transComplVolume,maskF);
	
	//use the plan of the volume to transform the volume back to real space		
	//process transformation to real space
	this->backwardFtPlan->execute(*this->transComplVolume,*this->transVolume);
	
	if(this->meanVolume == NULL)	
			this->meanVolume = new tom::Volume<T>(this->volume->getSizeX(), this->volume->getSizeY(), this->volume->getSizeZ(), NULL,NULL);
	
	//copy the result into the destination and shift it
	tom::fftshift(*this->transVolume,*this->meanVolume,true);

	//divide by the number of elements
	this->meanVolume->shift_scale(0.,1./((((T) this->numel()) * mask->getInnerSum())));
}

/****************************************************************************//**
 * \brief Writes the content of the os3_volume to file.
 *
 * 
 * Writes the content of the os3_volume to file. But only the content of the volume object.
 *******************************************************************************/
template<typename T> void tom::os3_volume<T>::writeToEm(std::string fileName){

	this->volume->write_to_em(fileName.c_str(),NULL);

}


/****************************************************************************//**
 * \brief Calculates the FFT of the volume.
 *
 * 
 * Calculates the FFT of the volume.
 *******************************************************************************/
template<typename T> void tom::os3_volume<T>::transform(){
	
	if(this->volume == NULL){
		throw std::runtime_error("transform : The volume has not been initialised before transformation!\n");
	}

	if(this->transVolume == NULL || this->transComplVolume == NULL){
		throw std::runtime_error("transform : The transform volumes have not been initialised before transformation!\n");
	}

	//create plan for volume if needed
	if(this->forwardFtPlan == NULL){
		this->forwardFtPlan = new tom::fftw::Plan<T>(*transVolume, *transComplVolume, FFTW_MEASURE);
		this->backwardFtPlan= new tom::fftw::Plan<T>(*transComplVolume, *transVolume, FFTW_MEASURE);
	}

	if(this->fourierVolume == NULL){
		if(volume->getSizeZ() == 1){
			this->fourierVolume = new tom::Volume<std::complex<T> >(this->volume->getSizeX(), this->volume->getSizeY()/2+1,1, NULL,NULL);
		}
		else{
			this->fourierVolume = new tom::Volume<std::complex<T> >(this->volume->getSizeX(), this->volume->getSizeY(), this->volume->getSizeZ()/2+1, NULL,NULL);
		}
	}
		
	//execute the plan
	this->transVolume->setValues(*volume);
	
	this->forwardFtPlan->execute(*transVolume,*transComplVolume);

	//the result is stored in fourierVolume and may be used for later correlations
	this->fourierVolume->setValues(*transComplVolume);

}


/****************************************************************************//**
 * \brief Transforms the volume object back to real space.
 *
 * 
 * Calculates the real space representation of the fourier volume stored.
 *******************************************************************************/
template<typename T> void tom::os3_volume<T>::transformBack(bool normalise){

	//check whether the object has already been transformed
	if(this->volume == NULL || this->transVolume == NULL || this->transComplVolume == NULL){
		throw std::runtime_error("transformBack : The object can not be transformed back because it has not been initialised yet!\n");
	}

	

	//copy the fourierVolume to the transformation volume
	this->transComplVolume->setValues(*this->fourierVolume);

	//process transformation to real space
	this->backwardFtPlan->execute(*transComplVolume,*transVolume);
	//process ifftshift and copy the values to the volume
	tom::fftshift(*transVolume,*this->volume,true);

	if(normalise){
		double numberOfVoxels =((double)this->getSizeX())*((double)this->getSizeY())*((double)this->getSizeZ());
		//divide by the number of voxels
		this->volume->shift_scale(0.,1./((T) numberOfVoxels));
	}
}


/****************************************************************************//**
 * \brief Checks whether two objects have the same size
 *
 * 
 * Checks whether two objects have the same size.
 *******************************************************************************/
template<typename T> bool tom::os3_volume<T>::equalSize(const tom::os3_volume<T>* volume){

	return (this->volume->is_equal_size(*volume->getVolume()));

}







template class tom::os3_volume<TFLOAT>;
template class tom::os3_volume<TDOUBLE>;

