/*=================================================================
 *
 * xraycorrect.c	Performs an in-place xray-correction
 *	                
 *
 * The syntax is:
 *
 *		tom_xraycorrectc(IMG,FACTOR)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 30.09.2008
 * By: Stephan Nickell
 * Revision: 1.00 by AK, added calculation of std and mean
 *
 *=================================================================*/
#include "mex.h"
#include "math.h"

/* Input Arguments */
#define	IMG     prhs[0]
#define	STD	    prhs[1]

float calc_std(
short *img,
int dim_x,
int dim_y,
float mean
)
{
    int i;
    float std;
    
    std = 0;
    
    for (i=0;i<dim_x*dim_y;i++){
        std = std + (((float)img[i] - mean)*((float)img[i] - mean));
    }
    
    std = sqrt(std / ((float)dim_x * (float) dim_y));
    return std;
}

float calc_mean(
short *img,
int dim_x,
int dim_y
)
{
    int i;
    float mean;
    
    mean = 0;    
    
    for (i=0;i<dim_x*dim_y;i++){
        mean = mean + (float)img[i];
    }
    
    mean = mean / ((float) dim_x * (float) dim_y);
    
    return mean;
}

void xraycorrect(
short *img,
float mean,
float std,
int dim_x,
int dim_y
)
{
    int ix, iy, iix, iiy, idx, count, pos_x, pos_y;
    float mm,mp, val;
    
    mm=mean-std;
    mp=mean+std;
    count=0;
    
    for (ix=0;ix<dim_x;ix++){
        for (iy=0;iy<dim_y;iy++){
            if((float)img[ix+iy*dim_y]<=mm || (float)img[ix+iy*dim_y]>=mp){
            /* interpolate */
                val=0; idx=0; count++;
                for (iix=-1;iix<2;iix++){
                    for (iiy=-1;iiy<2;iiy++){
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
    }
    mexPrintf("Xray correction applied. %i value(s) outside range: %f to %f. Mean",count,mm,mp);
    
}


void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray*prhs[] )

{
    mxArray *img;
    const mwSize *pp_dims;
    int factor;
    mwSize dim[2], ndim;
    float mean, std;
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments required.\n Syntax: xraycorrect(IMG,MEAN,FACTOR)\n IMG: uncorrected image\n FACTOR: standard deviation times factor\n");
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
    
    factor = mxGetScalar(STD);
    if (factor == 0) {
        mexErrMsgTxt("STD is 0. Image is constant ?!\n");
    }
    
    pp_dims=mxGetDimensions(IMG);
    dim[0]=pp_dims[0];
    dim[1]=pp_dims[1];
    
    /* Assign pointers to the various parameters */
    img = mxGetData(IMG);
    mean = calc_mean(img,dim[0],dim[1]);
    std = calc_std(img,dim[0],dim[1],mean);
    /* Do the actual computations in a subroutine */
    xraycorrect(img,mean,std * factor,dim[0],dim[1]);
    mexPrintf(" %f and Std: ",mean);
    mexPrintf("%f\n",std);    
    return;
}
