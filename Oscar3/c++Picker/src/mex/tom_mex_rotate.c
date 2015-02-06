/****************************************************************************//**
 * \file tom_mex_rotate.c
 * \author Thomas Hrabe
 * \date Oct. 26 2007
 *******************************************************************************/



#include <string.h>

#ifdef _WIN32
    #define snprintf _snprintf
#endif

#include "tom/core/rotate.h"


#include "mex.h"

#define    INP    prhs[0]
#define    OUT    prhs[1]
#define    EULER  prhs[2]
#define    INT    prhs[3]
#define    CENT   prhs[4]
#define    FACT   prhs[5]




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
const mwSize *dims_i, *dims_o,*dimension_euler; /*dimension */
int dimension;
char *ip;
float *center, *euler_A;
float *fact;

#define __BUFFERLENGTH_SYNOPSIS__ 1024
char synopsis[__BUFFERLENGTH_SYNOPSIS__];


snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "%s(in, out, [angles], interpolation_type, [center])", mexFunctionName());
if (nrhs==0 && nlhs==0) {
    /* Print help */
    mexPrintf("SYNOPSIS: %s\n", synopsis);
    mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ ".\n");
    return;
}


/* Check for proper number of arguments */
if (nrhs != 6)
{

    mexErrMsgTxt("5 input arguments required.\n Syntax: rot3d(input_image,output_image,[Angle(s)],Interpolation_type,[center])");
}
else
{
    if (nlhs > 1)
    {
        mexErrMsgTxt("Too many output arguments.");
    }

}

/* Check data types */
if (!mxIsSingle(INP) || !mxIsSingle(OUT))
{
    mexErrMsgTxt("Input volumes must be single.\n");
}

if (mxGetNumberOfDimensions(INP)!= mxGetNumberOfDimensions(OUT))
{
    mexErrMsgTxt("Image volumes must have same dimensions.\n");
}

ip = mxArrayToString(INT);
if (strcmp ("linear",ip) != 0 && strcmp("cspline",ip) != 0 && strcmp("cubic",ip) != 0)
{
    mexErrMsgTxt("Unknown interpolation type ... only linear,cspline, cubic interpolation implemented\n");
}


dimension=mxGetNumberOfDimensions(INP);
dims_i=mxGetDimensions(INP);
dims_o=mxGetDimensions(OUT);
center=(float *)mxGetData(CENT);
fact = (float *)mxGetData(FACT);
dimension_euler=mxGetDimensions(EULER);
euler_A=(float *)mxGetData(EULER);


if (dimension==3)
{


    if (dims_o[0]!=dims_i[0] || dims_o[1]!=dims_i[1] || dims_o[2]!=dims_i[2])
    {
        mexErrMsgTxt("Volumes must have same size.\n");
    }
    /* Do the actual computations in a subroutine */
    if(strcmp("linear",ip) == 0){
        rot3dlin((float *)mxGetData(INP), (float *)mxGetData(OUT),dims_i[0],dims_i[1],dims_i[2],(float *)mxGetData(EULER),center[0],center[1],center[2],dimension_euler[0]);
    }else if(strcmp("cspline",ip) == 0){
        rot3dcspline((float *)mxGetData(INP),(float *)mxGetData(OUT),dims_i[0],dims_i[1],dims_i[2],(float *)mxGetData(EULER),center[0],center[1],center[2],dimension_euler[0],fact);
    }else if(strcmp("cubic",ip) == 0){
        rot3dcubic((float *)mxGetData(INP),(float *)mxGetData(OUT),dims_i[0],dims_i[1],dims_i[2],(float *)mxGetData(EULER),center[0],center[1],center[2],dimension_euler[0]);
    }else{
        mexErrMsgTxt("Unknown interpolation type for 3D... only linear, cubic and cspline interpolation implemented\n");
    }
}


if (dimension==2)
{
    if (dims_o[0]!=dims_i[0] || dims_o[1]!=dims_i[1])
    {
        mexErrMsgTxt("Images must have same size.\n");
    }
    /* Do the actual computations in a subroutine */
    /*rot2d(mxGetData(INP),mxGetData(OUT),dims_i[0],dims_i[1],mxGetData(EULER),ip[0],center[0],center[1],dimension_euler[0]);*/

    if(strcmp("linear",ip) == 0){
        rot2dlin ((float *)mxGetData(INP),(float *)mxGetData(OUT),dims_i[0],dims_i[1],(float *)mxGetData(EULER),center[0],center[1],dimension_euler[0]);
    }else if(strcmp("cspline",ip) == 0){
        rot2dcspline ((float *)mxGetData(INP),(float *)mxGetData(OUT),dims_i[0],dims_i[1],(float *)mxGetData(EULER),center[0],center[1],dimension_euler[0],fact);
    }else if(strcmp("cubic",ip) == 0){
        rot2dcubic ((float *)mxGetData(INP),(float *)mxGetData(OUT),dims_i[0],dims_i[1],(float *)mxGetData(EULER),center[0],center[1],dimension_euler[0]);
    }else {
        mexErrMsgTxt("Unknown interpolation type for 2D... only linear, cubic and cspline interpolation implemented\n");
    }
}

}
