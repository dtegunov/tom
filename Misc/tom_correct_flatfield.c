/*=================================================================
 *
 * tom_correct_flatfield.c	Performs an interpolated flatfield correction
 *	                
 *
 * The syntax is:
 *
 *		tom_correct_flatfield(IMG,FLAT1,FLAT2,FACTOR,MEAN,FLAG)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 30.09.2008
 * By: Stephan Nickell
 * Revision: 1.00 by 
 *
 *=================================================================*/
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "mex.h"
#include "matrix.h"

/* Input Arguments */
#define	IMG     prhs[0]
#define	FLAT1	prhs[1]
#define	FLAT2	prhs[2]
#define	FACTOR	prhs[3]
#define	MEAN	prhs[4]
#define	FLAG	prhs[5]

static void tom_correct_flatfield( 
  short *img, 
  short *flat1, 
  short *flat2, 
  float factor, 
  float mean, 
  int dim_x,
  int dim_y
  ) 
{
int i;
float factor2;

factor2=1-factor;
printf("Flatfield correction applied\n"); 

 for (i=0;i<dim_x*dim_y;i++){
    img[i]=(float)img[i]/(((float)flat1[i]*factor+(float)flat2[i]*factor2)/mean);
 }

}
static void tom_correct_darkfield( 
  short *img, 
  short *flat1, 
  short *flat2, 
  float factor, 
  float mean, 
  int dim_x,
  int dim_y
  ) 
{
int i;
float factor2;

factor2=1-factor;
printf("Darfield correction applied\n"); 

 for (i=0;i<dim_x*dim_y;i++){
    img[i]=(float)img[i]-(((float)flat1[i]*factor+(float)flat2[i]*factor2));
 }

}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    short *img;
    short *flat1;
    short *flat2;
    double *p_offs;
    const mwSize *pp_dims;
    float factor; 
    float mean; 
    int flag; 
    size_t dim[3];
    size_t ndim; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 6) { 
	mexErrMsgTxt("Six input arguments required.\n Syntax: tom_correct_flatfield(IMG,FLAT1,FLAT2,FACTOR,MEAN,FLAG)\n IMG: uncorrected image\n FLAT1,FLAT2: flatfields which will be interpolated\n MEAN: mean value of IMG\n FLAG: 0 if dark-, 1 if flatfield correction is applied.\n"); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    ndim = mxGetNumberOfDimensions(prhs[1]);
    if (ndim != 2) { 
	mexErrMsgTxt("Two dimensional image as input requiered.\n"); 
    }
    if ( mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSingle(prhs[0]) )
                mexErrMsgTxt("Input image must be of type int16.\n");

    /* Check the dimensions of IN. */ 
    /* Get the number of dimensions in the input argument. */
    
    factor = mxGetScalar(FACTOR); 
    mean = mxGetScalar(MEAN); 
    flag = mxGetScalar(FLAG); 

    pp_dims=mxGetDimensions(prhs[0]);
    dim[0]=pp_dims[0];
    dim[1]=pp_dims[1];

    /* Assign pointers to the various parameters */ 
    img = mxGetData(IMG);
    flat1 = mxGetData(FLAT1);
    flat2 = mxGetData(FLAT2);
    /* Do the actual computations in a subroutine */
if (flag==1) tom_correct_flatfield(img,flat1,flat2,factor,mean,dim[0],dim[1]);
    else tom_correct_darkfield(img,flat1,flat2,factor,mean,dim[0],dim[1]);

    return;    
}
