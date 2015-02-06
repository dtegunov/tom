/*=================================================================
 *
 * tom_bininc.c bins a input volume n times, the output is identical to tom_binc,
 * the time needed for binning is half.
 *
 * The calling syntax is:
 *
 *		[OUT] = tom_bininc(invol,[xbinning ybinning zbinning])
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 05/07/2006
 * By: Andreas Korinek
 *
 *=================================================================*/

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    
    unsigned int binning[3], ndims;
    int output_dims[3];
    const double *p_binning;
    mxArray *indata, *outdata;
    mxClassID datatype;
    const int  *input_dims;
    unsigned int idx=0, idy=0, idz=0, bidx, bidy, bidz, outpos=0; 
    float singlebinvalue = 0.0;
    float bin_vol;
    
    p_binning=mxGetData(prhs[1]);
    binning[0]=(unsigned int)p_binning[0];
    binning[1]=(unsigned int)p_binning[1];
    binning[2]=(unsigned int)p_binning[2];
   
            
    indata = mxGetData(prhs[0]);
    ndims=mxGetNumberOfDimensions(prhs[0]);
    datatype = mxGetClassID(prhs[0]);
    if(datatype != mxSINGLE_CLASS) { mexErrMsgTxt("Input volume must be single."); }
    
    if (ndims > 3 || ndims < 2) { mexErrMsgTxt("Input matrix must be 2D or 3D.\n"); }
    
    input_dims = mxGetDimensions(prhs[0]);
    output_dims[0] = (unsigned int) input_dims[0] / (unsigned int) binning[0]; 
    output_dims[1] = (unsigned int) input_dims[1] / (unsigned int) binning[1];
    
    if (ndims == 3) {
        output_dims[2] = (unsigned int) input_dims[2] / (unsigned int) binning[2];
    }
    else {
        output_dims[2] = 1;
        binning[2] = 1;
    }
    
    bin_vol = (float) (binning[0] * binning[1] * binning[2]);


    if ((plhs[0] = mxCreateNumericArray(3,output_dims,mxSINGLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_binc.\n"); }
    outdata = mxGetData(plhs[0]);
    /*3d data*/
    if (output_dims[2] > 1) {
        for(idz=0;idz<=input_dims[2]-binning[2];idz = idz + binning[2]) {
            for(idy=0;idy<=input_dims[1]-binning[1];idy = idy + binning[1]) {
                for(idx=0;idx<=input_dims[0]-binning[0];idx = idx + binning[0]) {

                    for(bidz=0;bidz<binning[2];bidz++) {
                        for(bidy=0;bidy<binning[1];bidy++) {
                            for(bidx=0;bidx<binning[0];bidx++) {
                                singlebinvalue = ((float*)indata)[idx+bidx + (idy+bidy)*input_dims[0] + (idz+bidz)*input_dims[1]*input_dims[0]] + singlebinvalue;
                            }
                        }
                    }
                    ((float*)outdata)[outpos] = singlebinvalue / bin_vol;
                    singlebinvalue = 0.0;
                    outpos++;
                }
            }
        }
    }
    /*2d data*/
    else {
        
        for(idy=0;idy<=input_dims[1]-binning[1];idy = idy + binning[1]) {
            for(idx=0;idx<=input_dims[0]-binning[0];idx = idx + binning[0]) {
                
                for(bidy=0;bidy<binning[1];bidy++) {
                    for(bidx=0;bidx<binning[0];bidx++) {
                        singlebinvalue = ((float*)indata)[idx+bidx + (idy+bidy)*input_dims[0]] + singlebinvalue;
                    }
                }
                
                ((float*)outdata)[outpos] = singlebinvalue / bin_vol;
                singlebinvalue = 0.0;
                outpos++;
            }
        }
    }
}
