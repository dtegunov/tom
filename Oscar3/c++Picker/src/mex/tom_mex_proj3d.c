/***********************************************************************//**
 * \file tom_mex_proj3d.c
 * \brief MEX-Function for proj3d_solid_XXX in tom_proj3d. Replaces tom_proj3d.m/tom_proj3d2c.m.
 * \author  Thomas Haller
 * \version 0.1
 * \date    04.11.2007
 **************************************************************************/



#include <stdio.h>
#ifdef _WIN32
    #define snprintf _snprintf
#endif


#include "mex.h"



#include "tom_mex_helpfcn.h"

#include "tom/core/transform.h"










/***********************************************************************//**
 * \brief Project a 3D-Volume to an image.
 *
 * This the a mex-interface to the function proj3d_solid_zxz_GENERIC
 * as defined in tom_transform.c.
 * It does the same as
 **************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double phi = 0., theta = 0.;
    const mxArray *PRHS_VOL = NULL;
    const mxArray *PRHS_ANGLES = NULL;
    const mxArray *PRHS_NORMALIZE = NULL;
    mxArray *PLHS_IMG = NULL;
    int numberOfDimensions;
    const mwSize *dimensions;
    mwSize NX, NY, NZ, NXYZ;
    int res;
    int normalize = 0;
    double normalize_defval;
    mxClassID classID;
    #define __STRBUFFER_MAX_LENGTH 1024
    char strbuffer[__STRBUFFER_MAX_LENGTH+1];

    #define __BUFFERLENGTH_SYNOPSIS__ 1024
    char synopsis[__BUFFERLENGTH_SYNOPSIS__];

    snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "img = %s(volume, [angles_deg, normalize])", mexFunctionName());
    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis);
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ "\n");
        return;
    }


    if (nrhs < 1 || nrhs > 3) {
        snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s requires up to 3 input arguments.\nSyntax: PROJ = %s(VOL, [ANGLES_DEG, NORMALIZE]);", mexFunctionName(), mexFunctionName());
        mexErrMsgTxt(strbuffer);
    }
    if (nlhs > 1) {
        snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s returns only the projected volume.\nSyntax: PROJ = %s(VOL, [ANGLES_DEG, NORMALIZE]);", mexFunctionName(), mexFunctionName());
        mexErrMsgTxt(strbuffer);
    }

    {
        /* check inputparameters/types. */
        mwSize i, j;
        mxArray *tmpArray;
        PRHS_VOL = prhs[0];
        numberOfDimensions = mxGetNumberOfDimensions(PRHS_VOL);
        dimensions = mxGetDimensions(PRHS_VOL);
        NX = dimensions[0];
        NY = dimensions[1];
        NZ = numberOfDimensions>=3 ? dimensions[2] : 1;
        NXYZ = NX*NY*NZ;
        classID = mxGetClassID(PRHS_VOL);

        if (!mxIsNumeric(PRHS_VOL) || mxIsComplex(PRHS_VOL) || !(mxIsDouble(PRHS_VOL) || mxIsSingle(PRHS_VOL)) || mxIsSparse(PRHS_VOL) ||
            (numberOfDimensions<2||numberOfDimensions>3) || NXYZ<1) {
            snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s: Input volume must be a non empty, real floating point volume.", mexFunctionName());
            mexErrMsgTxt(strbuffer);
        }

        if (nrhs >= 2) {

            PRHS_ANGLES = prhs[1];
            i = mxGetM(PRHS_ANGLES); j = mxGetN(PRHS_ANGLES);
            if (!mxIsNumeric(PRHS_ANGLES) || mxIsComplex(PRHS_ANGLES) || mxGetNumberOfDimensions(PRHS_ANGLES)!=2 ||
                !((i==0&&j==0) || ((i==1||j==1)&&(i*j<=2))) || (i!=0 && !(tmpArray=getDoubleArray(PRHS_ANGLES)))) {
                snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s: The second parameter must be the rotation angles phi and theta in degree.", mexFunctionName());
                mexErrMsgTxt(strbuffer);
            }
            if (i!=0) {
                phi = mxGetPr(tmpArray)[0];
                if (i*j==2) {
                    theta = mxGetPr(tmpArray)[1];
                }
                mxDestroyArray(tmpArray);
            }
        }

        if (nrhs >= 3) {
            PRHS_NORMALIZE = prhs[2];
            i = mxGetM(PRHS_NORMALIZE); j = mxGetN(PRHS_NORMALIZE);
            if (!mxIsNumeric(PRHS_NORMALIZE) || mxIsComplex(PRHS_NORMALIZE) || mxGetNumberOfDimensions(PRHS_NORMALIZE)!=2 ||
                !((i==0&&j==0) || (i==1&&j==1)) || (i!=0 && !(tmpArray=getDoubleArray(PRHS_NORMALIZE)))) {
                snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s: The third parameter specifies whether to normalize. No normalizing is [] (default), otherwise specify the default value for the untouched pixels.", mexFunctionName());
                mexErrMsgTxt(strbuffer);
            }
            if ((normalize = i!=0)) {
                normalize_defval = mxGetPr(tmpArray)[0];
                mxDestroyArray(tmpArray);
            }
        }
    }


    if (!(PLHS_IMG=mxCreateNumericMatrix(NX, NY, classID, mxREAL))) {
        snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s: Error creating numeric matrix for output (probably out of memory)", mexFunctionName());
        mexErrMsgTxt(strbuffer);
    }

    if (classID == mxDOUBLE_CLASS) {
        res = tom_transf_proj3d_solid_zxz_double((const double *)mxGetData(PRHS_VOL), NX, NY, NZ, phi, theta, (double *)mxGetData(PLHS_IMG), normalize, normalize_defval);
    } else {
        res = tom_transf_proj3d_solid_zxz_float((const float *)mxGetData(PRHS_VOL), NX, NY, NZ, phi, theta, (float *)mxGetData(PLHS_IMG), normalize, (float)normalize_defval);
    }

    if (!res) {
        snprintf(strbuffer, __STRBUFFER_MAX_LENGTH, "%s: Unexpected error during projection", mexFunctionName());
        mexErrMsgTxt(strbuffer);
    }

    if (nlhs>=1) {
        plhs[0] = PLHS_IMG;
    }
}



