
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

#define PI 3.141592653589793238



void tom_proj2volc(float *VOL,const float *PROJ,int NX,int NY,int NZ,double PHI,double THETA,int NXP,int NYP,int x_offset,int y_offset,int z_offset) {

/*  Define variables */
int I1, I2, I3;
double T[4][4];
size_t index_v;
size_t index_p;
int iz, iy, ix;
float RADD, RXADD, RYADD, RX, RY;
int IPX, IPY;
float XD, YD;
float XIN1,XIN2,YIN;

/*  Calculate rotation matrix */
T[0][0] = 0.0;T[0][1] = 0.0;T[0][2] = 0.0;T[0][3] = 0.0;
T[1][0] = 0.0;T[1][3] = 0.0;T[2][3] = 0.0;T[3][3] = 0.0;
T[2][0] = 0.0;T[3][2] = 0.0;T[3][0] = 0.0;
T[1][1] = cos(THETA)*cos(PHI);
T[2][1] = cos(THETA)*sin(PHI);
T[3][1] = -sin(THETA);
T[1][2] = -sin(PHI);
T[2][2] = cos(PHI);

/* Main loop over all pixels of the object*/
index_v=0;
index_p=0;

XIN1 = 0.0;
XIN2 = 0.0;
YIN = 0.0;

for (iz=1;iz<=NZ;iz++){
              
    RADD = (iz-(NZ/2+1-z_offset))*T[3][1];
              
    for (iy=1;iy<=NY;iy++){
                   
         RXADD = (iy-(NY/2+1-y_offset))*T[2][1] + RADD;
         RYADD = (iy-(NY/2+1-y_offset))*T[2][2];
                   
         for (ix=1;ix<=NX;ix++){  
                   
         RX = (ix-(NX/2+1-x_offset))*T[1][1] + RXADD + (NXP/2+1);
         RY = (ix-(NX/2+1-x_offset))*T[1][2] + RYADD + (NYP/2+1);
                   
         IPX = floor(RX);
         IPY = floor(RY);    
         
         XD = RX - IPX;
         YD = RY - IPY;
                   
         if (IPX < 1 || IPY < 1 || IPX > NXP || IPY > NYP) {index_v = index_v + 1;}
                   
         /*  Bilinear interpolation */
         else{
                   
         index_p = (IPY-1)*NXP + IPX-1;
                   
         if (IPX==NXP) {XD=0.0;}
                    
              XIN1 = PROJ[index_p] + (PROJ[index_p+1]-PROJ[index_p])*XD;
              
                   if (IPY==NYP){YD=0.0; XIN2=0.0;}
                   else{XIN2 = PROJ[index_p+NXP] + (PROJ[index_p+NXP+1]-PROJ[index_p+NXP]) *XD;}
                   
              YIN = XIN1 + (XIN2-XIN1) *YD;
                   
              VOL[index_v] = VOL[index_v] + YIN;
              
              index_v = index_v + 1;
                   
              }
         
           }
         }
       }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]){ 
    
    /* Get projection */
    const mxArray *mxPROJ = prhs[0];
    const float *pPROJ = (const float *)mxGetData(mxPROJ);
    
    /* Get dimensions of projection*/
    const mwSize *dimensions_proj = mxGetDimensions(mxPROJ);
    const int NXP = dimensions_proj[0];
    const int NYP = dimensions_proj[1];
    
    /* Get angles and convert to radians*/
    const mxArray *mxANGLES = prhs[1];
    const double *pANGLES = (const double *)mxGetData(mxANGLES);
    const double PHI = (PI/180)*pANGLES[0];
    const double THETA = (PI/180)*pANGLES[1];
    
    /* Get dimvol*/
    const mxArray *mxDIMVOL = prhs[2];
    const double *pDIMVOL = (const double *)mxGetData(mxDIMVOL);
    const int NX = pDIMVOL[0];
    const int NY = pDIMVOL[1];
    const int NZ = pDIMVOL[2];
    
    /* Get offset*/
    const mxArray *mxOFFSET = prhs[3];
    const double *pOFFSET = (const double *)mxGetData(mxOFFSET);
    const double OFFSET_X = pOFFSET[0];
    const double OFFSET_Y = pOFFSET[1];
    const double OFFSET_Z = pOFFSET[2];
    
    /* Create volume */
    float *pVOL = NULL;
    mxArray *mxVOL = NULL;
    mwSize dims[3];
    dims[0] = NX;
    dims[1] = NY;
    dims[2] = NZ;
    mxVOL = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
    pVOL = mxGetData(mxVOL);
    
    /* Calculate volume */
    tom_proj2volc(pVOL,pPROJ,NX,NY,NZ,PHI,THETA,NXP,NYP,OFFSET_X,OFFSET_Y,OFFSET_Z);
    
    /* Return volume */
    plhs[0] = mxVOL;
}
