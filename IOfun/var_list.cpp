//#include "stdafx.h"
#include "var_list.h"

template <class T>
var_list<T>::var_list(void)
{
}

template <class T>
var_list<T>::var_list(T* vpData, int lDimensions, int lType, int DimCount) {
	this->lDim			= new int[DimCount];
	this->lDimensions	= lDimensions;
	this->lType			= lType;
	this->vpData		= vpData;
}

template <class T>
var_list<T>::~var_list(void)
{
	
}

