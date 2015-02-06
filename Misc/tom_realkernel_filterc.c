/*=================================================================
 *
 * tom_realkernel_filterc.c	performs an in-place real space kernel filter
 *	                
 *
 * The syntax is:
 *
 *		tom_realkernel_filterc(IMG,SIZE)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 23.10.2008
 * By: Stephan Nickell
 *
 *=================================================================*/
#include "mex.h"
#include "math.h"

/* Input Arguments */
#define	IMG     prhs[0]
#define	SIZE	    prhs[1]


void realkernel_filter(
short *img,
int size,
int dim_x,
int dim_y
)
{
    int ix, iy, iix, iiy, idx, count, pos_x, pos_y,size_2;
    float val;
    
    size_2=size/2;
    
    for (ix=0;ix<dim_x-size;ix++){
        for (iy=0;iy<dim_y-size;iy++){
            /* interpolate */
                val=0; idx=0; count++;
                for (iix=-size_2;iix<size_2;iix++){
                    for (iiy=-size_2;iiy<size_2;iiy++){
                        pos_x = ix+iix;
                        pos_y = iy+iiy;
                        if (pos_x > -1 && pos_x < dim_x && pos_y > -1 && pos_y < dim_y) {
                            if ((float)img[ix+iy*dim_y+iix+iiy]>mm && (float)img[ix+iy*dim_y+iix+iiy]<mp){
                                val=(float)img[ix+iy*dim_y+iix+iiy]+val;
                                idx++;
                            }
                        }
                    }
                }
                val=val/(float)idx;
                img[ix+iy*dim_y]=(short)val;
        }
    }
    mexPrintf("Xray correction applied. %i value(s) outside range: %f to %f. Mean",count,mm,mp);
    
}


void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray*prhs[] )

{
    mxArray *img;
    const mwSize *pp_dims;
    int size;
    mwSize dim[2], ndim;
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments required.\n Syntax: tom_realkernel_filterc(IMG,SIZE)\n IMG: image\n SIZE: kernel size in pixel\n");
    } else if (nlhs > 0) {
        mexErrMsgTxt("No output arguments.");
    }
    ndim = mxGetNumberOfDimensions(IMG);
    if (ndim != 2) {
        mexErrMsgTxt("Two dimensional image as input required.\n");
    }
    if ( !mxIsInt16(IMG))
        mexErrMsgTxt("Input image must be of type int16.\n");
    
    /* Check the dimensions of IN. */
    /* Get the number of dimensions in the input argument. */
    
    size = mxGetScalar(SIZE);
    if (factor == 0) {
        mexErrMsgTxt("SIZE is 0.\n");
    }
    
    pp_dims=mxGetDimensions(IMG);
    dim[0]=pp_dims[0];
    dim[1]=pp_dims[1];
    
    /* Assign pointers to the various parameters */
    img = mxGetData(IMG);
    /* Do the actual computations in a subroutine */
    realkernel_filter(img,size,dim[0],dim[1]);
    return;
}
