/*
%TOM_VOL2PROJ calculates a 2D projection of a 3D volume. C implementation.
%
%   proj=tom_vol2proj(vol,angles)
%
%PARAMETERS
%
%  INPUT
%   vol            3D volume (single or double, in C code handled as float), 
%   angeles    a 2-vector [phi theta] (single or double,  in C code handled as double) 
%                    in degrees, which is describing the orientation of the projection in 3d-space.
%                    For example, phi is the projection direction and theta the tiltangle.
%  
%  OUTPUT
%   proj          2D projection (single, in C code handled as float)
%
%EXAMPLE
%   cyl = tom_cylinder(8, 10, [32 32 32]); % Creates volume
%   [cylproj] = tom_vol2proj(cyl,[45 45]); % Calculates projection 
%
%REFERENCES
%   EM package, R. Hegerl
%
%SEE ALSO
%   TOM_BACKPROJ3D
%
%   created by ME 01/02/08
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
*/

#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

#define PI 3.141592653589793238

void tom_vol2projc(int NX, int NY, int NZ, int NXP, int NYP, const float *VOL, double PHI, double THETA, float *PROJ) {
    
    /* Define variables */
    double R[3][3];
    size_t index_v;
    size_t index_p;
    int iz, iy, ix;
    float X[3], XP[2], DIPX, DIPY, W;
    int IPX, IPY;
    
    /*  Calculate projection direction */
    R[0][0] = cos(THETA)*cos(PHI)*cos(PHI) + sin(PHI)*sin(PHI);
    R[1][0] = cos(THETA)*cos(PHI)*sin(PHI) - cos(PHI)*sin(PHI);
    R[2][0] = sin(THETA)*cos(PHI);
    
    R[0][1] = cos(THETA)*cos(PHI)*sin(PHI) - cos(PHI)*sin(PHI);
    R[1][1] = cos(THETA)*sin(PHI)*sin(PHI) + cos(PHI)*cos(PHI);
    R[2][1] = sin(THETA)*sin(PHI);
    
    R[0][2] = -sin(THETA)*cos(PHI);
    R[1][2] = -sin(THETA)*sin(PHI);
    R[2][2] = cos(THETA);
    
    /* Main loop over all pixels of the object*/
    index_v=0;
    index_p=0;
    
    X[2] = -1.0;
    for (iz=1; iz<=NZ; iz++) {
         X[2] = X[2] + 1.0;
         X[1] = -1.0;
         for (iy=1; iy<=NY; iy++) {
              X[1] = X[1] + 1;
              X[0] = -1.0;
              for (ix=1; ix<=NX; ix++) {
                   X[0] = X[0] + 1.0;
                   
                   XP[0] = (X[0] - (NX/2))*R[0][0] + (X[1] - (NY/2))*R[0][1] + (X[2] - (NZ/2))*R[0][2] + (NXP/2);
                   XP[1] = (X[0] - (NX/2))*R[1][0] + (X[1] - (NY/2))*R[1][1] + (X[2] - (NZ/2))*R[1][2] + (NYP/2);
                   
                   if ((XP[0] <= 0) | (XP[0] >= NXP-1) | (XP[1] <= 0) | (XP[1] >= NYP-1)) {index_v = index_v + 1;}
                     
                   /*  Bilinear interpolation */
                   else{
                   
                        IPX = floor(XP[0]);
                        IPY = floor(XP[1]);
                        
                        DIPX = XP[0] - IPX;
                        DIPY = XP[1] - IPY;
                        
                        W = VOL[index_v];
                        
                        index_p = IPX + 1 + IPY*NXP;
                        
                        PROJ[index_p-1] = PROJ[index_p-1] + (1.-DIPX) * (1.-DIPY) * W;
                        PROJ[index_p] = PROJ[index_p] + DIPX * (1.-DIPY) * W;
                        
                        index_p = index_p + NXP;
                        
                        PROJ[index_p-1] = PROJ[index_p-1] + (1.-DIPX) * DIPY * W;
                        PROJ[index_p] = PROJ[index_p] + DIPX * DIPY * W;
                        
                        index_v = index_v + 1;}
                               
            }
        }
    }
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Get volume */
    const mxArray *mxVOL = prhs[0];
    const float *pVOL = (const float *)mxGetData(mxVOL);
    
    /* Get dimensions of volume*/
    const mwSize *dimensions = mxGetDimensions(mxVOL);
    const int NX = dimensions[0];
    const int NY = dimensions[1];
    const int NZ = dimensions[2];
    
    /* Set dimensions of projection */
    const int NXP = NX;
    const int NYP = NY;
    
    /* Get angles and convert to radians*/
    const mxArray *mxANGLES = prhs[1];
    const double *pANGLES = (const double *)mxGetData(mxANGLES);
    const double PHI = (PI/180)*pANGLES[0];
    const double THETA = (PI/180)*pANGLES[1];
    
    /* Create projection */
    float *pPROJ = NULL;
    mxArray *mxPROJ = NULL;
    mwSize dims[2];
    dims[0] = NXP;
    dims[1] = NYP;
    mxPROJ = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    pPROJ = mxGetData(mxPROJ);
    
    /* Calculate projection */
    tom_vol2projc(NX, NY, NZ, NXP, NYP, pVOL, PHI, THETA, pPROJ);
    
    /* Return projection */
    plhs[0] = mxPROJ;

}
