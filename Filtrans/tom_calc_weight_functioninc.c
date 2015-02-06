/*=================================================================
 *
 * tom_calc_weight_function calculates 2D and yes yes 3D Weighting function
 * The syntax is:
 *
 * w_func=tom_calc_weight_fucntion(DIMS,ALL_ANGLES,THICKNESS,ANGLE_PROJ)
 *
 * INPUT:
 *     DIMS:            Dimension of the weighting function ... works as a switch between 2d and 3d
 *     ALL_ANGLES:      All Angles used in the backproj phi and the
 *     THICKNESS        Sample Thickness
 *     ANGLE_PROJ       parameter for 2d angle must be an element of ALL Angles
 *
 *  OUTPUT
 *      out:             weighting function 2d or 3d as reduced complex ... see alse tom_complete
 *
 *
 *  Example:
 *   2d: w_func=tom_calc_weight_functioninc([30 30],[10 20; 30 40; 50 10]',30,[10 20]);
 *   3d: w_func=tom_calc_weight_functioninc([30 30 30],[10 20; 30 40; 50 10]',30);
 *
 *  NOTE:
 *  All input types are double. Returned w_func is single
 *  All_Angles matrix must be transposed ==> [10 20; 30 40; 10 10;]'
 *
 * see also Radermacher et. al., 1986, A new 3-D reconstruction scheme
 * applied to to the 50S ribosomal subunit of E.coli, J.Microssc.
 * 141:RP1-RP2
 * and Electron Tomography: Three-Dimensional Imaging with the TEM, edited
 * by J. Frank, M.Radermacher, p91-115
 *
 *
 * Last changes: NOV27 2005
 * SN and FB after a Fortran routine by RH
 *
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tom_calc_weight_functioninc.h"
#ifdef MATLAB
#include "mex.h"
#endif


void weight2d(float *BR, int *dim, double *angles, int dimangles, double *thickness, double *angles_projection) {
    
    int i, I;
    int ITLT, IND, IT;
    float VST, ZW, UST;
    float ARG, ADD, SUM;
    int IX, IY;
    float XST, YST, ZST;
    int dimx, dimy, NDIM;
    double SPHI, CPHI, STHE, CTHE, PHI, THETA, A1, A2;
    double THICK;
    float **TILT;
    
    
    SPHI = sin(angles_projection[0]*PI/180);
    CPHI = cos(angles_projection[0]*PI/180);
    STHE = sin(angles_projection[1]*PI/180);
    CTHE = cos(angles_projection[1]*PI/180);
    A1 = CTHE*CPHI;
    A2 = CTHE*SPHI;
    
    
    IND = 0;
    VST = 0.0;
    ITLT = 0;
    dimx=dim[0];
    dimy=dim[1];
    NDIM=dimx;
    THICK=thickness[0];
    IT=1;
    
    
    /* starting to count at 1, feels like yehaa fortran77*/
    TILT=calloc(dimangles+1, sizeof(*TILT));
    for (i=0;i < dimangles+1; i++) {
        TILT[i]=(float*)calloc(3, sizeof(*TILT));
    }
    
    for (I=1; I < dimangles+1 ; I++) {
        PHI=(float)angles[ITLT]*PI/180;
        ITLT=ITLT+1;
        THETA=(float)angles[ITLT]*PI/180;
        ZW=(float)sin(THETA)*PI*(float)THICK;
        TILT[I][1]=(float)ZW*(float)cos(PHI);
        TILT[I][2]=(float)ZW*(float)sin(PHI);
        TILT[I][3]=(float)cos(THETA)*PI*(float)THICK;
        ITLT=ITLT+1;
    }
    
    
    
    for (IY=1; IY < dimy+1; IY++) {
        UST=0.0;
        for(IX=1; IX < dimx+1 ; IX++) {
            
            BR[IND]=1;
            SUM=0;
            ITLT=1;
            XST = UST*A1 - VST*SPHI;
            YST = UST*A2 + VST*CPHI;
            ZST = -UST*(float)STHE;
            
            for (I=1;I < dimangles +1; I++) {
                
                ARG=XST*TILT[ITLT][1]+YST*TILT[ITLT][2]+ZST*TILT[ITLT][3];
                if (fabs(ARG) < (float) 1e-6) {
                    ADD=1;
                }
                else if ((fabs(ARG) > 1e-6) && (fabs(ARG) < PI)) {
                    ADD=sin(ARG)/ARG;
                }
                else {
                    ADD=0;
                }
                
                SUM=SUM+ADD;
                ITLT=ITLT+1;
                
            } /* for (I=1;I < dimangles +1; I++) */
            if (SUM>1.0)
            {BR[IND]=BR[IND]/(float)SUM;}
            else
            {BR[IND]=BR[IND]/1.0;}
            IND=IND+1;
            UST=UST+(float)1.0/(float)NDIM;
            if (IX==dimx/2) {
                UST=-UST;
            }            
        }
        VST=VST+(float)1.0/(float)NDIM;
        if (IY==dimy/2) {
            VST=-VST;
        }
    }
    
    /*free memory*/
    for (i=0;i < dimangles+1; i++) {
        free(TILT[i]);
    }
    free(TILT);
    
}



void weight3d(float *BR, int *dim, double *angles, int dimangles, double *thickness) {
    
    int size_stream;
    int dimx, dimy, dimz, NDIM;
    int THICK;
    int IT=0;
    int ITLT=0;
    int IX, IY, IZ;
    float XST, YST, ZST;
    float ARG, ADD, SUM;
    float ZW;
    int IND;
    int I, i;
    float WIDTH, pi;
    float PHI, THETA;
    float **TILT;
    
    
    WIDTH=(float)3.1416;
    pi=(float)3.1416;
    ZW=(float)0.0;
    THICK=(int)thickness[0];
    
    /* starting to count at 1, feels like yehaa fortran77*/
    TILT=calloc(dimangles+1, sizeof(*TILT));
    for (i=0;i < dimangles+1; i++) {
        TILT[i]=(float*)calloc(3, sizeof(*TILT));
    }
    
    for (I=1;I <= dimangles; I++) {
        PHI=(float)angles[ITLT]*pi/180;
        ITLT=ITLT+1;
        THETA=(float)angles[ITLT]*pi/180;
        ZW=(float)sin(THETA)*pi*THICK;
        TILT[I][1]=(float)ZW*(float)cos(PHI);
        TILT[I][2]=(float)ZW*(float)sin(PHI);
        TILT[I][3]=(float)cos(THETA)*pi*THICK;
        ITLT=ITLT+1;
    }
    
    
    dimx=dim[0];
    dimy=dim[1];
    dimz=dim[2];
    NDIM=dimx;
    XST=0;
    IND = 0;
    ZST = 0;
    IT=0;
    for (IZ=1;IZ < dimz+1; IZ++) {
        
        YST=0;
        for (IY=1; IY < dimy+1; IY++) {
            XST=0;
            for (IX=1;IX<dimx+1; IX++) {
                SUM=0;
                ITLT=IT+1;
                for (I=1; I < dimangles+1; I++) {
                    ARG=XST*TILT[ITLT][1]+YST*TILT[ITLT][2]+ZST*TILT[ITLT][3];
                    
                    if (fabs(ARG) < (float)1e-6) {
                        ADD=1;
                    }
                    else if (fabs(ARG) > 1e-6 && fabs(ARG) < WIDTH)
                        
                    {
                        ADD=(float)sin(ARG)/ARG;
                    }
                    else {
                        ADD=0;
                    }
                    SUM=SUM+ADD;
                    ITLT=ITLT+1;
                }
                if (SUM>1.0)
                {BR[IND]=1.0/SUM;}
                else
                {BR[IND]=1.0;}
                IND=IND+1;
                XST=XST+(float)1.0/dimx;
                if (IX==dimx/2) {
                    XST = -XST;
                }
                
            }
            YST = YST + (float)1.0/dimy;
            if (IY==dimy/2) {
                YST = -YST;
            }
            
        }
        ZST=ZST+(float)1.0/dimz;
        if (IZ==dimz/2) {
            ZST=-ZST;
        }
        
        
    }/*for (IZ=1;IZ < dimz; IZ++) */
    
    
}




#ifdef MATLAB

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    const int *dimsp, *dimensionp, *dim_anglesp;
    int dims[3];
    double *help;
    
    dimensionp=mxGetDimensions(DIM_WFUNK);
    dim_anglesp=mxGetDimensions(ALL_ANGLES);
    
    /* funny construction to feed mxCreateNumericArray*/
    /* const int * sucks */
    help=mxGetData(DIM_WFUNK);
    dims[0]=(int)help[0];
    dims[1]=(int)help[1];
    dims[2]=(int)help[2];
    dimsp=(const*)&dims;
    
    
    OUT =mxCreateNumericArray(dimensionp[1], dimsp, mxSINGLE_CLASS, mxREAL);
    
    if (dimensionp[1]==3) {
        weight3d(mxGetData(OUT), dims, mxGetData(ALL_ANGLES), dim_anglesp[1], mxGetData(THICKNESS));
    }
    
    else {
        weight2d(mxGetData(OUT), dims, mxGetData(ALL_ANGLES), dim_anglesp[1], mxGetData(THICKNESS), mxGetData(ANGLE_PROJ));
    }
    
    
}


#endif




