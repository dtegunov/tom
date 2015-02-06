/*=================================================================
*
* tom_comparec.c    Comparsion of 2 fourier transforms 
* The syntax is:
*
*        comp(VOL1,VOL2)
*
*
* Last changes: Oct. 17, 2003
* M. Riedlberger
*
*=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#define	 PI ((double)3.14159265358979323846264338327950288419716939937510)

/* Input Arguments */
#define    VOL1    prhs[0]
#define    VOL2    prhs[1]
#define    VOLD    prhs[2]
#define    NUM     prhs[3]
#define    START   prhs[4]
#define    STOP    prhs[5]
#define    OUT     prhs[6]

void tom_comparec( 
  float	*vol1, 
  float	*vol2, 
  float	*vold, 
  long	size,
  long  num,
  long  start,
  long  stop,
  float	*out)
{
long  center=(long) (size+1)/2;	/* center  */
long  i, j, k;
long  size2 = size*size;
float *ampl1 = (float *) malloc(sizeof(float)*num);
float *ampl2 = (float *) malloc(sizeof(float)*num);
float *ampl_diff = (float *) malloc(sizeof(float)*num);
float *phares1 = (float *) malloc(sizeof(float)*num);
float *phares2 = (float *) malloc(sizeof(float)*num);
float ccc;
float amplitude1, amplitude2, amplitude_diff;
long  *n_shell = (long *) malloc(sizeof(long)*size);
float length = (float)(size-center) / (float)num;	/* length of shells */
float shell_float;
int   shell;
int   ps=0;
long  index;
float mean,rmsd,phares,phares_v,res,arg,delta,sqrn;

for (i=0; i<num; i++) {
     ampl1[i] = 0;
     ampl2[i] = 0;
     ampl_diff[i] = 0;
     phares1[i] = 0;
     phares2[i] = 0;
     n_shell[i] = 0;
    }
     
for (k=0; k<size; k++)
     for (j=0; j<size; j++)
          for (i=0; i<=(size/2); i++) {
	       /* distance to center, number of shell */
	       shell=(int) (sqrt((i-center)*(i-center)+(j-center)*(j-center)+(k-center)*(k-center))/length);
/*	       shell = (int)shell_float;
	       if ((float) shell == shell_float)
	           shell--;*/
	       if (shell>=num || shell<0 || (i==center&&j==center&&k==center)) continue;
	       n_shell[shell]++;
	       
	       index = i+j*size+k*size2;
	       amplitude1=vol1[index];
	       amplitude2=vol2[index];
	       amplitude_diff=vold[index];
	       phares_v = sqrt(amplitude1)+sqrt(amplitude1);
	       ampl1[shell] += amplitude1;
	       ampl2[shell] += amplitude2;
	       ampl_diff[shell] += amplitude_diff;
	       phares1[shell] += phares_v;

	       arg = 2*sqrt(amplitude1*amplitude2);			/* calculation of cos(theta) for Phase Residuum */
	       if (arg>0) {
	           arg=(amplitude1+amplitude2-amplitude_diff)/arg;
		   if (arg>1) arg=1;
		   if (arg<-1) arg=-1;
		   }
	       delta = acos(arg);   				/* Phaseshift */
	       phares2[shell] += phares_v*delta*delta;		/* Phase Residuum */
	  }
	  printf("Shell  N  Res.  RMS-Dev, Abs.  Rel.  Amplitude1  Amplitude2  Mean F  CCC  2/SQRT(N)  Phares\n");
for (i=0; i<num; i++) {
     ccc = (ampl1[i]+ampl2[i]-ampl_diff[i])/sqrt(ampl1[i]*ampl2[i])/2;
     mean = (ampl1[i]+ampl2[i]-ampl_diff[i]) / (n_shell[i]*2 );
     if(mean<0) mean=-mean;
     mean=sqrt(mean);
     rmsd = sqrt(ampl_diff[i]/(2*(ampl1[i]+ampl2[i])-ampl_diff[i]));
/*     a=sqrt((abs(ampl_diff1[i]+ampl2[i]-ampl_diff[i])/(n_shell[i]*2 )));*/
     ampl1[i] = sqrt(ampl1[i]/n_shell[i]);
     ampl2[i] = sqrt(ampl2[i]/n_shell[i]);
     ampl_diff[i] = sqrt(ampl_diff[i]/n_shell[i]);
     sqrn = 2/sqrt(n_shell[i]);					/* 2 / sqrt(N) */
     phares = sqrt(phares2[i]/phares1[i])*180/PI;		/* Phase Residuum in degrees */
     res = 2*(float)num/((float)i+1);				/* pixel resolution */

     out[i] = i+1;
     out[i+num] = n_shell[i];
     out[i+num*2] = res;
     out[i+num*3] = ampl_diff[i];
     out[i+num*4] = rmsd;
     out[i+num*5] = ampl1[i];
     out[i+num*6] = ampl2[i];
     out[i+num*7] = mean;
     out[i+num*8] = ccc;
     out[i+num*9] = sqrn;
     out[i+num*10] = phares;
     printf("%f %f %5f %5f %5f %5f %5f %f %f %f %f\n",out[i],out[i+num],out[i+num*2],out[i+num*3],out[i+num*4],out[i+num*5],out[i+num*6],out[i+num*7],out[i+num*8],out[i+num*9],out[i+num*10]);
     
     ps+=n_shell[i];
    }
     
printf("%d \n",ps);
free(ampl1);
free(ampl2);
free(ampl_diff);
free(phares1);
free(phares2);
free(n_shell);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   const int *dims;

   /* Check for proper number of arguments */
   if (nrhs != 7) {
       mexErrMsgTxt("7 input arguments required.\n Syntax: proj3d(vol1, vol2, vold, num, start, stop, out)");    }
   else if (nlhs > 1) {printf("%d\n",nrhs);
       mexErrMsgTxt("Too many output arguments.");    }
 
   /* Check the dimensions */
   if (mxGetNumberOfDimensions(VOL1)!=3 || mxGetNumberOfDimensions(VOL1)!=3) {
       mexErrMsgTxt("3 dimensional images as input required.\n");    }

   /* Check data types */
   if (!mxIsSingle(VOL1) || !mxIsSingle(VOL2)) {
       mexErrMsgTxt("Input volumes must be single.\n"); }

   dims=mxGetDimensions(VOL1);

   /* Do the actual computations in a subroutine */
   comp(mxGetData(VOL1),mxGetData(VOL2),mxGetData(VOLD),dims[0],mxGetScalar(NUM),mxGetScalar(START),mxGetScalar(STOP),mxGetData(OUT));
   return;
   }

