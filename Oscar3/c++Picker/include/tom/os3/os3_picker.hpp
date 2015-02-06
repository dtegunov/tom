#ifndef __TOM_OS3_PICKER_HPP__
#define __TOM_OS3_PICKER_HPP__


#include "tom/os3/os3_picklist.hpp"
#include "tom/os3/os3_structures.hpp"

namespace tom{
	
	template<typename T> void os3_correlationPicker(tom::os3_job &job,bool printStatus);
	template<typename T> void os3_createPicklist(tom::os3_job &job,bool printStatus);
	template<typename T> void os3_picker(tom::os3_job &job,bool printStatus);


}



#endif









