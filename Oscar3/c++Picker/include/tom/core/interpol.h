/***********************************************************************//**
 * \file interpol.h
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    23.10.2007
 **************************************************************************/
#ifndef ___INCLUDE_CORE__INTERPOL_H__
#define ___INCLUDE_CORE__INTERPOL_H__



#ifdef __cplusplus
    extern "C" {
#endif


#define __TOM_INTERPOL_NEAREST          ((int)1)
#define __TOM_INTERPOL_LINEAR           ((int)2)
#define __TOM_INTERPOL_CUBIC            ((int)3)
#define __TOM_INTERPOL_CSPLINE          ((int)4)



float interpolate2D(const float* image, long sizeX, long sizeY,  float x,float y, short method);
float interpolate3D(const float* image, long sizeX, long sizeY, long sizeZ, float x,float y, float z, short method);


#ifdef __cplusplus
    }
#endif


#endif
