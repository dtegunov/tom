/***********************************************************************//**
 * \file tom_os3_correlation.cpp
 * \brief Different correlation functions
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    08.01.2008
 **************************************************************************/
#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_volume.hpp"
//----------------------------------
#include "tom/core/fftw_plan.hpp"
#include "tom/core/volume.hpp"
#include "tom/core/volume_fcn.hpp"

#include <iostream>

namespace tom{
/****************************************************************************//**
 * \brief Processes the fast localised correlation function (Roseman 2003)
 *
 * \param[in] searchV The em image / volume containing the searched patterns
 * \param[in] patterV The pattern searched in searchV
 *
 * A correlation of two volumes.
 *******************************************************************************/
	template<typename T> 
	void os3_correlate(tom::os3_volume<T>* searchV, tom::os3_volume<T>* patternV,const bool normalise){

		tom::Volume<std::complex<T> >* searchF;
		tom::Volume<std::complex<T> >* patternF;

		if(! searchV->equalSize(patternV)){
			throw std::runtime_error("os3_correlate : The size of both volumes does not match.\n");
		}
		
		
		searchF = searchV->getFourierVolume();
		
		patternF = patternV->getFourierVolume();
		
		//correlation, the result is written to the patterF
		tom::element_wise_conj_multiply(*patternF,*searchF);
		
		patternV->transformBack(normalise);

	}


/****************************************************************************//**
 * \brief Processes the fast localised correlation function (Roseman 2003)
 *
 * \param[in] searchV The em image / volume containing the searched patterns
 * \param[in,out] patternV The pattern searched in searchV
 * \param[in] maskV   The mask used for the local correlation
 *
 * Please read the Roseman 2003 paper for a detailed description of the FLCF.
 * Comments in the code refer to the program only and not to the algorithm.
 * The result is written into the patterV object. No new memory is allocated for the calculation.	
 *
 * Size of the parameters: pattern and mask must be of same size, the searchV may differ. Pattern and mask are extended to searchV size if they differ.
 *******************************************************************************/
	template<typename T> 
	void os3_flcf(tom::os3_volume<T>* searchV,tom::os3_volume<T>* patternV, tom::os3_volume<T>* maskV){

		patternV->applyMask(maskV);

		patternV->calculateStatisticVolumes(maskV);
		
		tom::Volume<T> tm(*patternV->getMeanVolume());
		tom::Volume<T> ts(*patternV->getStdVolume());
		
		T patternMean = patternV->getCentralMean();
		T patternStd = patternV->getCentralStd();

		
		if(!searchV->equalSize(maskV)){
			
			tom::os3_volume<T> biggerMask(maskV);
			
			biggerMask.resize(searchV);	
			
			searchV->calculateStatisticVolumes(&biggerMask);
			
		}
		else
			searchV->calculateStatisticVolumes(maskV);

		if(! searchV->equalSize(patternV))
			patternV->resize(searchV);

		tom::Volume<T> meanVolume(*searchV->getMeanVolume());
		tom::Volume<T> stdVolume(*searchV->getStdVolume());

		os3_correlate(searchV,patternV,false);
		tom::Volume<T>* correlationResult = patternV->getVolume();
		correlationResult->shift_scale((T)0.0,(T)(1.0/(searchV->numel() * (double)maskV->getInnerSum())));

		correlationResult->shift_scale((T)0.0,(T)10000.0);

		/* mu_seachV * mu_patternV */
		meanVolume.shift_scale((T)0.0,patternMean);
		meanVolume.shift_scale((T)0.0,(T)10000.0);

		/*std_searchV * std_patternV*/
		stdVolume.shift_scale((T)0.,patternStd); // patternStd
		stdVolume.shift_scale((T)0.0,(T)10000.0);

		tom::element_wise_sub(*correlationResult,meanVolume);

		tom::element_wise_div(*correlationResult,stdVolume,(T)0.0);

		patternV->replaceCurrentPlan(searchV);
		
	}



/****************************************************************************//**
 * \brief Processes the fast localised correlation function (Roseman 2003)
 *
 * \param[in] searchV The em image / volume containing the searched patterns
 * \param[in] patterV The pattern searched in searchV
 * \param[in] maskV   The mask used for the local correlation
 *
 * Please read the vanHeel 1993 paper for a detailed description of the MCF.
 * Comments in the code refer to the program only and not to the algorithm.	 
 *******************************************************************************/
	template<typename T> 
	void os3_mcf(tom::os3_volume<T>* searchV,tom::os3_volume<T>* patternV,tom::os3_volume<T>* maskV){

		

	}

/****************************************************************************//**
 * \brief Processes the phase only correlation function (Roseman 2003)
 *
 * \param[in] searchV The em image / volume containing the searched patterns
 * \param[in] patterV The pattern searched in searchV
 * \param[in] maskV   The mask used for the local correlation
 *
 * Please read the Correlation Pattern Recognition book for a detailed description of the POF.
 * Comments in the code refer to the program only and not to the algorithm.	 
 *******************************************************************************/
	template<typename T> 
	void os3_psf(tom::os3_volume<T>* searchV,tom::os3_volume<T>* patternV,tom::os3_volume<T>* maskV){

		

	}


/****************************************************************************//**
 * \brief Processes the second order correlation function 
 *
 * \param[in] corrV The correlation volume storing the flcf peaks.
 * \param[in] autoc The autocorrelation of the pattern searched.
 * \param[in] maskV   The mask used for the local correlation.
 *
 * Please read the Correlation Pattern Recognition book for a detailed description of the POF.
 * Comments in the code refer to the program only and not to the algorithm.	 
 *******************************************************************************/
	template<typename T> 
	void os3_soc(tom::os3_volume<T>* corrV,tom::os3_volume<T>* autoV,tom::os3_volume<T>* maskV){

		corrV->resetStatistics(maskV);
		autoV->resetStatistics(maskV);

		tom::os3_flcf(corrV,autoV,maskV);

	}


/****************************************************************************//**
 * \brief Processes the peak to sidelobe ratio (Kumar et al. 2005 [Correlation Pattern Recognition])
 *
 * \param[in] xcfV The correalation image / volume containing the calculated normalised correlation coefficients
 * \param[in] patterV The pattern searched in searchV
 *
 * Please read the Correlation Pattern Recognition book for a detailed description of the PSR.
 * Comments in the code refer to the program only and not to the algorithm.	 
 *******************************************************************************/
	template<typename T> 
	void os3_psr(tom::os3_volume<T>* xcfV,tom::os3_volume<T>* patternV){

		tom::Volume<T> psrMask(*xcfV->getVolume());
		psrMask.setValues((T)0);
	
		//psrMask.printInfo("psrMask");
		std::size_t sx = xcfV->getSizeX();
		std::size_t sy = xcfV->getSizeY();		
		std::size_t sz = xcfV->getSizeZ();

		std::size_t px = patternV->getSizeX();
		std::size_t py = patternV->getSizeY();		
		std::size_t pz = patternV->getSizeZ();

		std::size_t posX = sx/2-px/2;
		std::size_t posY = sy/2-py/2;
		std::size_t posZ = sz/2-pz/2;
		if(sz <= 1)
			posZ = 0;
		
		//create window of same size as the pattern at the center of psrMask
		tom::Volume<T> windowVolume(psrMask,&psrMask.get(posX,posY,posZ),
									px,py,pz, 
									psrMask.getStrideX(), psrMask.getStrideY(), psrMask.getStrideZ());

		//windowVolume.printInfo("windowVolume");
		float outerRadius = patternV->getSizeX()/2;
		//create sphere with outer radius size
		tom::init_spheremask(windowVolume,outerRadius,(float).0,(float).0);

		float innerRadius = patternV->getSizeX() /4;
		if(innerRadius <5){
			innerRadius = 5;
		}

		sx = windowVolume.getSizeX()/2;
		sy = windowVolume.getSizeY()/2;		
		sz = windowVolume.getSizeZ()/2;

		posX = sx - (std::size_t)innerRadius;
		posY = sy - (std::size_t)innerRadius;
		posZ = sz - (std::size_t)innerRadius;
		if(sz <= 1)
			posZ = 0;

		//create an volume of same size as the pattern
		tom::Volume<T> window2(windowVolume);
		window2.setValues((T)0);
		//leap into the new volume
	
		tom::Volume<T> smallWindow(window2,&window2.get(posX,posY,posZ),
								   (std::size_t)innerRadius*2,(std::size_t)innerRadius*2,((sz<=0) ? (1):((std::size_t)innerRadius*2)),
								   window2.getStrideX(), window2.getStrideY(), window2.getStrideZ());
		//create a small sphere there
	
		tom::init_spheremask(smallWindow,innerRadius,(float).0,(float).0);
		
		//create a ring sphere in the big volume (psrMask)
		tom::element_wise_sub(windowVolume,window2);
		//psrMask.write_to_em("/fs/home/hrabe/develop/images/psr.em",NULL);
		//psr mask is finished
		
		//patternV now carries the psrMask
		patternV->setVolume(&psrMask);
		patternV->replaceCurrentPlan(xcfV);
		//std::cout << std::endl << __FILE__ << " " << __LINE__ << std::endl;
		tom::os3_volume<T> xcfVCopy(xcfV);
		//std::cout << std::endl << __FILE__ << " " << __LINE__ << std::endl;
		xcfVCopy.resetStatistics(patternV);
		tom::Volume<T>* meanVolume = xcfVCopy.getMeanVolume();
		tom::Volume<T>* stdVolume = xcfVCopy.getStdVolume();
		
		tom::Volume<T>* volume = xcfVCopy.getVolume();
	
		tom::element_wise_sub(*volume,*meanVolume);
		tom::element_wise_div(*volume,*stdVolume,(T)0.0);

		patternV->setVolume(&xcfVCopy);
	}

/****************************************************************************//**
 * \brief Applies a point spread function to an volume.
 *
 * \param[in] volume The volume modified. (Results overwrite the volume)
 * \param[in] psf The point spread function.
 *
 * Applies a point spread function to an volume. Same as MATLAB tom_apply_weight_function
 *******************************************************************************/
	template<typename T> 
	void os3_apply_weight_function(tom::os3_volume<T>& volume,const tom::os3_volume<T>& psf){

		if(!volume.getVolume()->is_equal_size(*psf.getVolume()))
			throw std::runtime_error("os3_apply_weight_function() : volume and psf do not have the same size!\n");

		//get the transformed fourier volume
		tom::Volume<std::complex<T> >* fVolume = volume.getFourierVolume();
		
		//create volume for shifted psf
		std::auto_ptr<tom::Volume<T> > psf_shifted(new tom::Volume<T>(psf.getSizeX(), psf.getSizeY(), psf.getSizeZ(), NULL,NULL));	
		std::auto_ptr<tom::Volume<T> > psf_shifted_;
		if(psf.getSizeZ() > 1)
			psf_shifted_.reset(new tom::Volume<T>(*psf_shifted,&(*psf_shifted).get(0,0,0), (*psf_shifted).getSizeX(), (*psf_shifted).getSizeY(), (*psf_shifted).getSizeZ()/2+1, 
												  (*psf_shifted).getStrideX(),(*psf_shifted).getStrideY(),(*psf_shifted).getStrideZ()));	
		else
			psf_shifted_.reset(new tom::Volume<T>(*psf_shifted, &(*psf_shifted).get(0,0,0),(*psf_shifted).getSizeX(), (*psf_shifted).getSizeY()/2+1, 1,
												  (*psf_shifted).getStrideX(),(*psf_shifted).getStrideY(),(*psf_shifted).getStrideZ()));	
		
		tom::fftshift(*psf.getVolume(),*psf_shifted,true);

		
		//apply weighting - element_wise_multiply 
		tom::element_wise_multiply(*fVolume,*psf_shifted_);
		
		
		tom::Volume<std::complex<T> >* transComplexVolume = volume.getTransComplVolume();
		tom::Volume<T>* transVolume = volume.getTransVolume();
		
		//fVolume contains the weighted and transformed volume
		transComplexVolume->setValues(*fVolume);
		tom::fftw::Plan<T>* backwardPlan = volume.getBackwardPlan();
		
		backwardPlan->execute(*transComplexVolume,*transVolume);
		
		//reset the volume object with the new results
		
		volume.getVolume()->setValues(*transVolume);
	}

}





template void tom::os3_correlate(tom::os3_volume<TFLOAT>* searchV, tom::os3_volume<TFLOAT>* patternV, const bool normalise);
template void tom::os3_flcf(tom::os3_volume<TFLOAT>* searchV,tom::os3_volume<TFLOAT>* patternV, tom::os3_volume<TFLOAT>* maskV);
template void tom::os3_psr(tom::os3_volume<TFLOAT>* xcfV,tom::os3_volume<TFLOAT>* patternV);
template void tom::os3_soc(tom::os3_volume<TFLOAT>* corrV,tom::os3_volume<TFLOAT>* autoV,tom::os3_volume<TFLOAT>* maskV);
template void tom::os3_apply_weight_function(tom::os3_volume<TFLOAT>& volume,const tom::os3_volume<TFLOAT>& psf);

template void tom::os3_correlate(tom::os3_volume<TDOUBLE>* searchV, tom::os3_volume<TDOUBLE>* patternV, const bool normalise);
template void tom::os3_flcf(tom::os3_volume<TDOUBLE>* searchV,tom::os3_volume<TDOUBLE>* patternV, tom::os3_volume<TDOUBLE>* maskV);
template void tom::os3_psr(tom::os3_volume<TDOUBLE>* xcfV,tom::os3_volume<TDOUBLE>* patternV);
template void tom::os3_soc(tom::os3_volume<TDOUBLE>* corrV,tom::os3_volume<TDOUBLE>* autoV,tom::os3_volume<TDOUBLE>* maskV);
template void tom::os3_apply_weight_function(tom::os3_volume<TDOUBLE>& volume,const tom::os3_volume<TDOUBLE>& psf);

