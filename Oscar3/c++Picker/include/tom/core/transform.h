/****************************************************************************//**
 * \file transform.h
 * \author Thomas Haller
 * \date Oct. 26, 2007
 * \brief Transformation of volumes and images.
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__TRANSFORM_H__
#define ___INCLUDE_CORE__TRANSFORM_H__




#include <stddef.h>
#if defined __GNUC__ || defined __xlC__
    #include <stdint.h>
#else
    #include <system/pstdint.h>
#endif


#include "tom/core/interpol.h"







#ifdef __cplusplus
    extern "C" {
#endif



typedef struct tom_transf_st_interpol_param {
    const void *    src;
    size_t          sizex;
    size_t          sizey;
    size_t          sizez;
    size_t          stepy;
    size_t          stepz;
    double          x;
    double          y;
    double          z;
    void            *optparam;
} tom_transf_st_interpol_param;





int tom_transf_init_rotmat_3d(double *P, size_t stepy, size_t nrots, const int *axes, const double *angles);


int tom_transf_invert(const double *P, double *Pinv, int N);



/** \copydoc tom_transf_rotate3d_double */
int tom_transf_rotate3d_float(const float *src, float *dst, size_t sizex, size_t sizey, size_t sizez,
              double phi, double psi, double theta, double centerx, double centery, double centerz,
              float (*intfcn)(const tom_transf_st_interpol_param *param), void *optparam);
/** \copydoc tom_transf_rotate3d_double */
int tom_transf_rotate3d_double(const double *src, double *dst, size_t sizex, size_t sizey, size_t sizez,
              double phi, double psi, double theta, double centerx, double centery, double centerz,
              double (*intfcn)(const tom_transf_st_interpol_param *param), void *optparam);


/** \copydoc tom_transf_transf3d_double */
int tom_transf_transf3d_float(const float *src,
                   size_t srcsizex, size_t srcsizey, size_t srcsizez,
                   size_t srcstepy, size_t srcstepz,
                   float *dst,
                   size_t dstsizex, size_t dstsizey, size_t dstsizez,
                   size_t dststepy, size_t dststepz,
                   const double *P,
                   int is_affinity,
                   float valinf,
                   float (*intfcn)(const tom_transf_st_interpol_param *param),
                   void *optparam);


/** \copydoc tom_transf_transf3d_double */
int tom_transf_transf3d_double(const double *src,
                   size_t srcsizex, size_t srcsizey, size_t srcsizez,
                   size_t srcstepy, size_t srcstepz,
                   double *dst,
                   size_t dstsizex, size_t dstsizey, size_t dstsizez,
                   size_t dststepy, size_t dststepz,
                   const double *P,
                   int is_affinity,
                   double valinf,
                   double (*intfcn)(const tom_transf_st_interpol_param *param),
                   void *optparam);




/** \copydoc tom_transf_interpol3d_nearest_double */
float tom_transf_interpol3d_nearest_float(const tom_transf_st_interpol_param *param);
/** \copydoc tom_transf_interpol3d_nearest_double */
double tom_transf_interpol3d_nearest_double(const tom_transf_st_interpol_param *param);





typedef struct {
    float defaultval;
    short method;
} tom_transf_interpol3d_continous_float_param;
float tom_transf_interpol3d_continous_float(const tom_transf_st_interpol_param *param);


/** \copydoc tom_transf_interpol3d_trilinear_double */
float tom_transf_interpol3d_trilinear_float(const tom_transf_st_interpol_param *param);
/** \copydoc tom_transf_interpol3d_trilinear_double */
double tom_transf_interpol3d_trilinear_double(const tom_transf_st_interpol_param *param);



/** \copydoc tom_transf_proj3d_solid_zxz_double */
int tom_transf_proj3d_solid_zxz_float(const float *VOL, size_t NX, size_t NY, size_t NZ, double PHI, double THETA, float *_PROJ, int normalize, float defval);

/** \copydoc tom_transf_proj3d_solid_zxz_double */
int tom_transf_proj3d_solid_zxz_double(const double *VOL, size_t NX, size_t NY, size_t NZ, double PHI, double THETA, double *_PROJ, int normalize, double defval);

/** \copydoc tom_transf_proj3d_double */
int tom_transf_proj3d_float(const float *src, size_t srcsizex, size_t srcsizey, size_t srcsizez, size_t srcstepy, size_t srcstepz, const double *P, float *dst, size_t dstsizex, size_t dstsizey, size_t dststepy, double *weight);
/** \copydoc tom_transf_proj3d_double */
int tom_transf_proj3d_double(const double *src, size_t srcsizex, size_t srcsizey, size_t srcsizez, size_t srcstepy, size_t srcstepz, const double *P, double *dst, size_t dstsizex, size_t dstsizey, size_t dststepy, double *weight);


#ifdef __cplusplus
    }
#endif



#endif



