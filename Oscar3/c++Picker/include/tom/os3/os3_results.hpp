/****************************************************************************//**
 * \file tom_os3_volume.h
 * \brief The header file for the class os3_results.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    23.01.2008
 *******************************************************************************/ 


#ifndef __INCLUDE_TOM_OS3_RESULTS_HPP__
#define __INCLUDE_TOM_OS3_RESULTS_HPP__

#include "tom/core/volume.hpp"
#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_picklist.hpp"
namespace tom{

/****************************************************************************//**
 * \brief The os3_results class stores the processing results of the oscar system.
 * 
 * 
 * The os3_results class stores the processing results of the oscar system. \n 
 *******************************************************************************/

	template<typename T>
	class os3_results{

		public:
			os3_results();
			os3_results(std::size_t x,std::size_t y,std::size_t z);
			os3_results(tom::os3_job job);
			~os3_results();	

			void update(tom::Volume<T>* xcfVolume,tom::Volume<T>* psrVolume,tom::Volume<T>* socVolume,int angleIndex);
			void writeToEm(tom::os3_job job);
			void writeToEm(std::string str);
			void setPicklist(tom::os3_job& job,bool printStatus);
			
			tom::os3_picklist<T> getPicklist() const;
			
		private:
			tom::Volume<T>* xcfResult;
			tom::Volume<T>* psrResult;
			tom::Volume<T>* socResult;
			tom::Volume<int>* angResult;
			tom::os3_picklist<T> picklist;
	};
}


#endif



/** \brief Returns a copy of the current picklist. */
template<typename T>
inline tom::os3_picklist<T> tom::os3_results<T>::getPicklist() const {
    return this->picklist;
}








