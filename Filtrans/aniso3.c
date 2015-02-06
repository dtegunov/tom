/*  Matlab MEX file.

result = aniso3(input,psiFunction,sigma,iterations,lambda)

*/

/*
 *Copyright 2000, Stanford University

  Permission to use, copy, modify, and distribute this software and its
  documentation for any purpose and without fee is hereby granted,
  provided that the above copyright notice appear in all copies and that
  both that copyright notice and this permission notice appear in
  supporting documentation, and that the name of Stanford University not
  be used in advertising or publicity pertaining to distribution of the
  software without specific, written prior permission.  Stanford
  University makes no representations about the suitability of this
  software for any purpose.  It is provided "as is" without express or
  implied warranty.
 
 */


#include "mex.h"
#include <stdlib.h>
#include <string.h> 

#define TUKEYPSI 0
#define LORENTZPSI 1
#define LINEARPSI 2

#define DEFSIGMA 1.0
#define DEFITER 10
#define DEFLAMBDA 0.25

#define MAX(A,B) (((A)>(B)) ? (A) : (B) )
#define MIN(A,B) (((A)<(B)) ? (A) : (B) )

/* Psi functions */
double tukeyPsi(double x, double sigma)
{
double x2;
if (x<(-sigma))
	return 0.0;

if (x>sigma)
	return 0.0;

x2 = (x/sigma);
x2 *= x2;
x2 = 1-x2;
x2 *= x2;
return (x*x2);
}

double linearPsi(double x)
{
return (2.0*x);
}

double lorentzianPsi(double x, double sigma)
{
return (2.0*x/(2*sigma*sigma + x*x));
}

/* Caniso3 - MAIN C routine.
*/
void Caniso3(register double *vol, int *dims,
             int funtype, double sigma, int Niter, double lambda)
{
    register int x, y, z, iter, index; /* indexes and sizes */
    register int Zind, ZindU, ZindD, Xind, XindL, XindR;
    register double psi;
    register double *prevvol;
    int Nbytes;
    register int Ymax, Xmax, Zmax;
    register int indexYS, indexYN, indexXL, indexXR, indexZU, indexZD;
    
    Nbytes = dims[0]*dims[1]*dims[2]*sizeof(double);
    Ymax = dims[0] - 1;
    Xmax = dims[1] - 1;
    Zmax = dims[2] - 1;
    prevvol = (double *)malloc(Nbytes);
    
    for (iter=0; iter<Niter; iter++) {
     /* copying the current vol to the previous vol */
        memcpy(prevvol, vol, Nbytes);
        for (z=0;z<dims[2];z++) {
            Zind = z*dims[0]*dims[1];
            ZindU = MAX(0, z-1)*dims[0]*dims[1];
            ZindD = MIN(Zmax, z+1)*dims[0]*dims[1];
            for (x=0;x<dims[1];x++) {
                Xind = x*dims[0];
                XindL = MAX(0, x-1)*dims[0];
                XindR = MIN(Xmax, x+1)*dims[0];
                for (y=0;y<dims[0];y++) {
                    index = Zind + Xind + y;
                    indexYS = Zind + Xind + MIN(Ymax, y+1);
                    indexYN = Zind + Xind + MAX(0, y-1);
                    indexXL = Zind + XindL + y;
                    indexXR = Zind + XindR + y;
                    indexZU = ZindU + Xind + y;
                    indexZD = ZindD + Xind + y;
                    switch (funtype) {
                        case TUKEYPSI:{
                            psi = tukeyPsi(prevvol[indexXL]-prevvol[index], sigma) +
                            tukeyPsi(prevvol[indexXR]-prevvol[index], sigma) +
                            tukeyPsi(prevvol[indexZU]-prevvol[index], sigma) +
                            tukeyPsi(prevvol[indexZD]-prevvol[index], sigma) +
                            tukeyPsi(prevvol[indexYN]-prevvol[index], sigma) +
                            tukeyPsi(prevvol[indexYS]-prevvol[index], sigma);
                            break;
                        }
                        case LORENTZPSI:{
                            psi = lorentzianPsi(prevvol[indexXL]-prevvol[index], sigma) +
                            lorentzianPsi(prevvol[indexXR]-prevvol[index], sigma) +
                            lorentzianPsi(prevvol[indexZU]-prevvol[index], sigma) +
                            lorentzianPsi(prevvol[indexZD]-prevvol[index], sigma) +
                            lorentzianPsi(prevvol[indexYN]-prevvol[index], sigma) +
                            lorentzianPsi(prevvol[indexYS]-prevvol[index], sigma);
                            break;
                        }
                        case LINEARPSI:{
                            psi = linearPsi(prevvol[indexXL]-prevvol[index]) +
                            linearPsi(prevvol[indexXR]-prevvol[index]) +
                            linearPsi(prevvol[indexZU]-prevvol[index]) +
                            linearPsi(prevvol[indexZD]-prevvol[index]) +
                            linearPsi(prevvol[indexYN]-prevvol[index]) +
                            linearPsi(prevvol[indexYS]-prevvol[index]);
                            break;
                        }
                        default: {break;}
                    }
                    vol[index] = prevvol[index] + lambda * psi;
                } /* end y loop */
            }  /* end x loop */
        }   /* end z loop */
   }    /* end main iterations loop */

	free(prevvol);
}


void mexFunction(int nlhs,   /* number of arguments on lhs */
		 mxArray	*plhs[],   /* Matrices on lhs      */
		 int nrhs,	   /* no. of mat on rhs    */
		 const mxArray	*prhs[]    /* Matrices on rhs      */
		 )
{
  const mxArray *vol;        /* input volume array */
  mxArray *result;        /* output volume array */
  double *voldata;    /* pointer to input volume data */
  double *resultdata;      /* pointer to the output matrix data */
  int ndimsvol;   /* number of dimensions of the input volume */
  int *dimsvol;   /* dimensions of the vol, array form */
  int funType; /* function type */
  char functionstr[256]; /* function string */
  double sigma, lambda;  /* input parameters */
  int Niter;             /* input parameter */
  int buflen;

  /* Check for proper number of arguments */
  if (nrhs==0) { /* help */
   printf("result = aniso3(vol,<psiFunction>,<sigma>,<iterations>,<lambda>)\n");
   printf("\n  vol: volume data (3D array)\n");
   return;
  }

  /* reading input parameters */
  vol = prhs[0];
  voldata = mxGetPr(vol);
  if (nrhs>1) {
   buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1]) * sizeof(mxChar)) + 1;
   mxGetString(prhs[1], functionstr, buflen);
   if (strcmp(functionstr, "tukeyPsi")==0)
     funType = TUKEYPSI;
   else if (strcmp(functionstr, "lorentzianPsi")==0)
     funType = LORENTZPSI;
   else if (strcmp(functionstr, "linearPsi")==0)
     funType = LINEARPSI;
   else
      mexErrMsgTxt("Invalid psi function");    
  }
  else {
   funType = TUKEYPSI;
  }
  if (nrhs>2)
   sigma = *mxGetPr(prhs[2]);
  else
   sigma = DEFSIGMA;
  if (nrhs>3)
   Niter = (int)*mxGetPr(prhs[3]);
  else
   Niter = DEFITER;
  if (nrhs>4)
   lambda = *mxGetPr(prhs[4]);
  else
   lambda = DEFLAMBDA;
  
  lambda = lambda/6.0; /* to do exactly the same as the m-file */

  /* geting volume size */
  ndimsvol = mxGetNumberOfDimensions(vol);
  if (ndimsvol!=3) {
      mexErrMsgTxt("vol must be a 3-D array.");
  }
  else {
      dimsvol = (int*) mxGetDimensions(vol);
  }

  /* creating output array */
  plhs[0] = mxDuplicateArray(prhs[0]);
  resultdata = mxGetPr(plhs[0]);

/* DEBUG */	
/*mexPrintf("Values:\n");
	mexPrintf("function: %i\n", funType);
	mexPrintf("function: %s\n", functionstr);
	mexPrintf("sigma: %f\n", sigma);
	mexPrintf("lambda: %f\n", lambda);
	mexPrintf("niter: %i\n", Niter);*/


  /* main routine */
  Caniso3(resultdata,  dimsvol, funType, sigma, Niter, lambda);
}
