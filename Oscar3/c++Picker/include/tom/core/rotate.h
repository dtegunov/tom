/****************************************************************************//**
 * \file rotate.h
 * \author Thomas Hrabe
 * \date Oct. 26 2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__ROTATE_H__
#define ___INCLUDE_CORE__ROTATE_H__


#ifdef __cplusplus
    extern "C" {
#endif


void rot2dlin(float *image, float *rotimg, long sx, long sy, float *p_phi, float px, float py, int euler_dim);
void rot2dcubic(float *image, float *rotimg, long sx, long sy, float *p_phi, float px, float py, int euler_dim);
void rot2dcspline(float *image, float *rotimg, long sx, long sy, float *p_phi, float px, float py, int euler_dim, float *fact);
void rot3dcspline(float *image, float *rotimg, long sx, long sy, long sz, float *p_euler_A, float px, float py, float pz, int euler_dim, float *fact);
void rot3dcubic(float *image, float *rotimg, long sx, long sy, long sz, float *p_euler_A, float px, float py, float pz, int euler_dim);
void rot3dlin(float *image, float *rotimg, long sx, long sy, long sz, float *p_euler_A, float px, float py, float pz, int euler_dim);


#ifdef __cplusplus
    }
#endif


#endif


