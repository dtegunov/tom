#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#ifdef MATLAB_MEX_FILE
    #include "mex.h"
    #define __memfct_calloc (&mxCalloc)
    #define __memfct_free (&mxFree)
    #define __memfct_realloc (&mxRealloc)
    #define __memfct_malloc (&mxMalloc)
#else
    #define __memfct_calloc (&calloc)
    #define __memfct_free (&free)
    #define __memfct_realloc (&realloc)
    #define __memfct_malloc (&malloc)
#endif


#define __RET_PROJ_NOVOLUME -1
#define __RET_PROJ_INVALID_ARGS -2
#define __RET_PROJ_MALLOC_FAILS -2
#define __RET_PROJ_OK 0

#define __PI ((double)3.141592653589793238512808959406186204433)

#ifndef __TYPE
#define __TYPE float
#endif

/***************************************************************
 * Calculates the projection of an 2d or an 3d-array in arbitrary direction.
 * if NZ == 0 than an 2d-array has to be projected.
 * otherwise an 3d-volumse
 * NX = Size in x-direction of input volume
 * NY = Size in y-direction of input volume
 * NZ = Size in z-direction of input volume
 * PHI = ProjectionDirection (in deg!!!!)
 * THETA = Tiltangle (in deg!!!!)
 ***************************************************************/
int proj(int NX, int NY, int NZ, const double *VOL, double PHI, double THETA, float **_PROJ, int *NNX, int *NNY) {

    if (!_PROJ || !NNX || !NNY) {
        return __RET_PROJ_INVALID_ARGS;
    }
    *_PROJ = NULL;
    *NNX = 0;
    *NNY = 0;

    if (NX <= 0 || NY <= 0 || !VOL) {
        return __RET_PROJ_NOVOLUME;
    }

    int NI; /* Do a 2 or a 3-dimensional projection???? */
    if (NZ <= 0) {
        NI = 2;
        NZ = 1;
    } else {
        NI = 3;
    }

    /*
        VARIABLES IN proj.f:
        NA: an inputmatrix of ints. contains in the JX(1) colum the parameters of the inputmatrix VOL (=ARRAY) and
            in the JX(2)th column the parameters of the outputarray PROJ (=A);
        JX[2]: Where which colums of the matrix NA are the parameters for the inputarray and the outputarray?????
        JX1, JX2: takes values of JX(1), JX(2)
        KA: Maybe the number of elements in the Array PROJ (=A)
    */

    int I1, I2, I3, I4; /* variables for loops */



    /* const int NXYZ = NX * NY * NZ; // Number of voxels */

    const int NXP = NX; /*  Size in x-direction of output projection */
    const int NYP = NY; /*  Size in y-direction of output projection */

    /*  initialization output projection */
    float *PROJ = (float *)(*__memfct_calloc)(NXP*NYP, sizeof(*PROJ));
    if (!PROJ) {
        return __RET_PROJ_MALLOC_FAILS;
    }
    *_PROJ = PROJ;
    *NNX = NXP;
    *NNY = NYP;

    /*  Default initialization */
    const double RLX = NX; /* The volume is sampled as there are RLX voxls in X-direction.*/
    const double RLY = RLX * NY / NX;
    const int H = NZ;


    const double DX = RLX / NX;
    const double DY = RLY / NY;
    const double DZ = H / NZ;

    const double DP[3] = { (NX/2.)*DX, (NY/2.)*DY, (NZ/2.)*DZ };

    const double XFACT = 1.0;
    const double YFACT = 1.0;

    const double DXP = DX * XFACT;
    const double DYP = DY * YFACT;

    const double DPNEW[2] = { NXP/2.*DXP, NYP/2.*DYP };

    const double RLXNEW = RLX * NXP / NX;
    const double RLYNEW = RLXNEW * NYP / NXP;


    /*  Calculation projection direction */

    const double PHIR = PHI * (__PI / 180);
    const double THETAR = THETA * (__PI / 180);

    const double CPHI = cos(PHIR);
    const double SPHI = sin(PHIR);
    const double CTHE = cos(THETAR);
    const double STHE = sin(THETAR);

    double DMdata[9];
    const double *DM = DMdata;
    if (NI == 2) {
        DMdata[0] = CPHI;
        DMdata[3] = -SPHI;
        DMdata[1] = SPHI;
        DMdata[4] = CPHI;
    } else {
        DMdata[0] = CTHE*CPHI*CPHI + SPHI*SPHI;
        DMdata[3] = CTHE*CPHI*SPHI - CPHI*SPHI;
        DMdata[6] = -STHE*CPHI;

        DMdata[1] = DMdata[3];
        DMdata[4] = CTHE*SPHI*SPHI + CPHI*CPHI;
        DMdata[7] = -STHE*SPHI;

        DMdata[2] = -DMdata[6];
        DMdata[5] = -DMdata[7];
        DMdata[8] = CTHE;
    }


    /*  Calculation projection point */

    double X[3] = { 0., 0., -DZ };
    double XP[2];

    double ZWX, ZWY, DIPX, DIPY, W;
    int IPX, IPY; /* saves IFIX(ZWX). Integervalue of the calculated 2D-Position (1D) */
    int INDP;

    int INDX = 0;

    for (I1=0; I1<NZ; I1++) {

        X[2] += DZ;
        X[1] = -DY;

        for (I2=0; I2<NY; I2++) {

            X[1] += DY;
            X[0] = -DX;

            for (I3=0; I3<NX; I3++) {

                X[0] += DX;

                XP[0] = 0.;
                XP[1] = 0.;

                for (I4=0; I4<NI; I4++) {
                    const double tmp = X[I4]-DP[I4];
                    XP[0] += tmp * DM[I4*3];
                    XP[1] += tmp * DM[I4*3+1];
                }

                XP[0] += DPNEW[0];
                XP[1] += DPNEW[1];

                if (NI == 2) {
                    XP[1] = 0.;
                }

                if ((XP[0] > 0.)  && (XP[0] < RLXNEW-DXP) &&
                    (XP[1] > 0.)  && (XP[1] < RLYNEW-DYP)) {

                    /*  Bilinear interpolation */

                    ZWX = XP[0] / DXP;
                    ZWY = XP[1] / DYP;

                    /*  attention! assuming that ZWX and ZWY can only be positve!. otherwise use fix!!!! */
                    IPX = floor(ZWX);
                    IPY = floor(ZWY);
                    #if 1
                    if (ZWX < 0 || ZWY<0) {
                        mexPrintf("ZWX or ZWY is not 0.!!!!!!\n");
                    }
                    #endif

                    DIPX = ZWX - IPX;
                    DIPY = ZWY - IPY;

                    W = VOL[INDX];

                    INDP = IPX + 1 + IPY*NXP;

                    PROJ[INDP-1] += (1.-DIPX) * (1.-DIPY) * W;
                    PROJ[INDP] += DIPX * (1.-DIPY) * W;

                    if (NI == 3) {
                        INDP += NXP;

                        PROJ[INDP-1] += (1.-DIPX) * DIPY * W;
                        PROJ[INDP] += DIPX * DIPY * W;
                    }
                }
                INDX++;
            }
        }
    }


    return __RET_PROJ_OK;

}








/***************************************************************
 * NX = Size in x-direction of input volume
 * NY = Size in y-direction of input volume
 * NZ = Size in z-direction of input volume
 * PHI = ProjectionDirection (in deg!!!!)
 * THETA = Tiltangle (in deg!!!!)
 ***************************************************************/
int tom_proj3d2(int NX, int NY, int NZ, const double *VOL, double PHI, double THETA, float *PROJ) {

    if (NX <= 0 || NY <= 0 || NZ <=0) {
        return __RET_PROJ_NOVOLUME;
    }

    /*  Variables */
    int I1, I2, I3, I4;

    /* const int NXYZ = NX * NY * NZ; // Number of voxels */

    const int NXP = NX; /*  Size in x-direction of output projection */
    const int NYP = NY; /*  Size in y-direction of output projection */

    /*  initialization output projection */
    float *pPROJ = PROJ;
    for (I1=0; I1<NYP; I1++) {
        for (I2=0; I2<NXP; I2++) {
            *pPROJ++ = 0.;
        }
    }

    /*  Default initialization */
    const int RLX = NX;
    const int H = NZ;

    const int RLY = RLX * NY / NX;

    const int DX = RLX / NX;
    const int DY = RLY / NY;
    const int DZ = H / NZ;

    const double DP[3] = { (NX/2.)*DX, (NY/2.)*DY, (NZ/2.)*DZ };

    const double XFACT = 1.0;
    const double YFACT = 1.0;

    const double DXP = DX * XFACT;
    const double DYP = DY * YFACT;

    const double DPNEW[2] = { NXP/2.*DXP, NYP/2.*DYP };

    const double RLXNEW = RLX * NXP / NX;
    const double RLYNEW = RLXNEW * NYP / NXP;


    /*  Calculation projection direction */

    const double PHIR = PHI * (__PI / 180);
    const double THETAR = THETA * (__PI / 180);

    const double CPHI = cos(PHIR);
    const double SPHI = sin(PHIR);
    const double CTHE = cos(THETAR);
    const double STHE = sin(THETAR);

    double DMdata[9];
    const double *DM = DMdata;
    DMdata[0] = CTHE*CPHI*CPHI + SPHI*SPHI;
    DMdata[3] = CTHE*CPHI*SPHI - CPHI*SPHI;
    DMdata[6] = -STHE*CPHI;

    DMdata[1] = DMdata[3];
    DMdata[4] = CTHE*SPHI*SPHI + CPHI*CPHI;
    DMdata[7] = -STHE*SPHI;

    DMdata[2] = -DMdata[6];
    DMdata[5] = -DMdata[7];
    DMdata[8] = CTHE;


    /*  Calculation projection point */

    double X[3] = { 0., 0., -DZ };
    double XP[2];

    double ZWX, ZWY, DIPX, DIPY, W;
    int IPX, IPY;

    int INDP; /* Index for the resulting image (the projection) */
    int INDX = 0; /* Running index for the source-array (VOL) */

    for (I1=0; I1<NZ; I1++) {

        X[2] += DZ;
        X[1] = -DY;

        for (I2=0; I2<NY; I2++) {

            X[1] += DY;
            X[0] = -DX;

            for (I3=0; I3<NX; I3++) {

                X[0] += DX;

                XP[0] = 0.;
                XP[1] = 0.;

                for (I4=0; I4<3; I4++) {
                    const double tmp = X[I4]-DP[I4];
                    XP[0] += tmp * DM[I4*3];
                    XP[1] += tmp * DM[I4*3+1];
                }

                XP[0] += DPNEW[0];
                XP[1] += DPNEW[1];

                if ((XP[0] > 0.)  && (XP[0] < RLXNEW-DXP) &&
                    (XP[1] > 0.)  && (XP[1] < RLYNEW-DYP)) {

                    /*  Bilinear interpolation */

                    ZWX = XP[0] / DXP;
                    ZWY = XP[1] / DYP;

                    /*  attention! assuming that ZWX and ZWY can only be positve!. otherwise use fix!!!! */
                    IPX = floor(ZWX);
                    IPY = floor(ZWY);
                    #if 0
                    if (ZWX < 0 || ZWY<0) {
                        mexPrintf("ZWX or ZWY is not 0.!!!!!!\n");
                    }
                    #endif

                    DIPX = ZWX - IPX;
                    DIPY = ZWY - IPY;

                    W = VOL[INDX];

                    INDP = IPX + 1 + IPY*NXP;

                    PROJ[INDP-1] += (1.-DIPX) * (1.-DIPY) * W;
                    PROJ[INDP] += DIPX * (1.-DIPY) * W;


                    INDP += NXP;

                    PROJ[INDP-1] += (1.-DIPX) * DIPY * W;
                    PROJ[INDP] += DIPX * DIPY * W;


                }

                INDX++;
            }
        }
    }


    return __RET_PROJ_OK;

}



#ifdef MATLAB_MEX_FILE


#define __MATLABSYNTAXSTR "Syntax: PROJ = tom_proj3d2c(VOL, ANGLES);"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if (nrhs != 2) {
        mexErrMsgTxt("tom_proj3d2c requires 2 input arguments.\n" __MATLABSYNTAXSTR);
    }
    if (nlhs != 1) {
        return;
    }

    /* check inputparameters/types. */

    const mxArray *mxVOL = prhs[0];
    const mxArray *mxANGLES = prhs[1];

    if (!mxIsNumeric(mxVOL) || mxIsComplex(mxVOL) || !mxIsDouble(mxVOL)) {
        mexErrMsgTxt("Input volume must be real double.\n");
    }


    if (!mxIsNumeric(mxANGLES) || mxIsComplex(mxANGLES)) {
        mexErrMsgTxt("Angles must be real double");
    } else if (!mxIsDouble(mxANGLES)) {
        mexErrMsgTxt("Angles must be double");
    }

    if (mxGetNumberOfDimensions(mxANGLES) != 2 ||
        !((mxGetM(mxANGLES)==1 && mxGetN(mxANGLES)==2) || (mxGetM(mxANGLES)==2 && mxGetN(mxANGLES)==1))) {
        mexErrMsgTxt("Angles must be a two-vector (in deg) (ANGLES = [PHI THETA])");
    }
    const double *pANGLES = (const double *)mxGetData(mxANGLES);
    const double PHI = pANGLES[0];
    const double THETA = pANGLES[1];

    const int numberOfDimensions = mxGetNumberOfDimensions(mxVOL);
    const mwSize *dimensions = mxGetDimensions(mxVOL);
    const int NX = dimensions[0];
    const int NY = dimensions[1];
    const int NZ = numberOfDimensions==3 ? dimensions[2] : 0;

    if (numberOfDimensions<2 || numberOfDimensions>3) {
        mexErrMsgTxt("Input volume must be 2 or 3 dimensional.");
    }

    if (NX <= 0 || NY <= 0 || (numberOfDimensions==3 && NZ <= 0)) {
        mexErrMsgTxt("Input volume is empty.");
    }

    const double *pVOL = (const double *)mxGetData(mxVOL);

    mxArray *mxPROJ = NULL;
    float *pPROJ = NULL;

    int proj_returncode;
    mwSize dims[2];
#ifdef __USE_MATLABVERSION

    dims[0] = NX;
    dims[1] = NY;
    if (!(mxPROJ=mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL))) {
        mexErrMsgTxt("Unable to create numeric matrix for output");
    }
    pPROJ = mxGetData(mxPROJ);;

    proj_returncode = tom_proj3d2(NX, NY, NZ, pVOL, PHI, THETA, pPROJ);

    if (proj_returncode != __RET_PROJ_OK) {
        mxDestroyArray(mxPROJ);
        mexErrMsgTxt("Unexpected error during projection!");
    }
#else
    int NNX=0, NNY=0;
    proj_returncode = proj(NX, NY, NZ, pVOL, PHI, THETA, &pPROJ, &NNX, &NNY);
    if (proj_returncode != __RET_PROJ_OK) {
        char s[100];
        sprintf(s, "Unexpected error during projection: errorcode %d", proj_returncode);
        mexErrMsgTxt(s);
    }

    dims[0] = NNX;
    dims[1] = NNY;
    if (!(mxPROJ=mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL))) {
        mxFree(pPROJ);
        mexErrMsgTxt("Unable to create numeric matrix for output");
    }
    mxFree(mxGetData(mxPROJ));
    mxSetData(mxPROJ, pPROJ);
#endif

    plhs[0] = mxPROJ;


}



#endif



