#pragma once
#ifndef _DATATYPES
	#include "DataTypes.h"
	#define _DATATYPES
#endif

template <class T>
class var_list
{
public:
	var_list(void);
	var_list(T*, int, int, int);
	virtual ~var_list(void);

	T*		vpData;
	int	lDimensions;
	int	lType;
	int*	lDim;

};
