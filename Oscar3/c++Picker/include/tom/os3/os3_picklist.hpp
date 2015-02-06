/****************************************************************************//**
 * \file tom_os3_volume.h
 * \brief The header file for the class os3_volume.
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    08.01.2008
 *******************************************************************************/ 
#ifndef __INCLUDE_OS3_PICKLIST_HPP__
#define __INCLUDE_OS3_PICKLIST_HPP__

#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_structures.hpp"
#include "tom/core/volume.hpp"
#include "tom/core/volume_fcn.hpp"
#include "tom/core/fftw_plan.hpp"

#include <iostream>
#include <algorithm>
namespace tom{

/****************************************************************************//**
 * \brief Stores a picklist.
 * 
 * 
 * Stores a picklist. Picklists can be merged, sorted, returned. 
 * Make sure piclists come from the same dataset and the patterns searched are same. 
 *******************************************************************************/

	template<typename T>
	class os3_picklist{
	
	
		public:
			
			os3_picklist();
			os3_picklist(const os3_picklist& original);
			os3_picklist(const std::string path);
			void loadPicklist(std::string path,bool clear);
			void push_back(tom::os3_pick<T>& pick);
			void push_back(std::vector<tom::os3_pick<T> >& picklist);
			void push_back(tom::os3_picklist<T>& other);
			void sortPicklist();
			
			std::vector<tom::os3_pick<T> > getPicklist() const;
			
			void createParticleStack();
			tom::Volume<T>* getParticleStack() const;
			
			void savePicklist(std::string path);
			
		private:
			
			std::vector<tom::os3_pick<T> > picklist;
			
			tom::Volume<T>* particleStack;
	};
}


/** \brief Returns the picklist stored*/
template<typename T>
inline std::vector<tom::os3_pick<T> > tom::os3_picklist<T>::getPicklist() const{
	return this->picklist;
}
/** \brief Appends a pick to the picklist*/
template<typename T>
inline void tom::os3_picklist<T>::push_back(tom::os3_pick<T> & pick){
	this->picklist.push_back(pick);
}

/** \brief Appends a pick to the picklist*/
template<typename T>
inline void tom::os3_picklist<T>::push_back(tom::os3_picklist<T>& other){
	
	std::vector<tom::os3_pick<T> > otherList = other.getPicklist();
	
	this->push_back(otherList);
}

/** \brief Concatenates a vector with the picklist*/
template<typename T>
inline void tom::os3_picklist<T>::push_back(std::vector<tom::os3_pick<T> >& picklist){
	this->picklist.insert( this->picklist.end(), picklist.begin(), picklist.end());
}

template<typename T> 
inline tom::Volume<T>* tom::os3_picklist<T>::getParticleStack() const{
	return this->particleStack;
}

template<typename T>
inline void tom::os3_picklist<T>::sortPicklist(){

	std::sort( this->picklist.begin(), this->picklist.end());

}

#endif







