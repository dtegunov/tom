/*=================================================================

*
* tom_dist.c  calculates the time-onsuming part of the exact weighting
* ... used to speed up the reconstruction with exact weighting
* The syntax is:
*
*        tom_dist(THETA,ACT_TH,PSI,D,XST,YST,ZST,W);
* 
*        THETA: used tilt anlges as single
*        ACT_TH: actual tilt angle as single
*        PSI: angle with the y-axis as single
*        XST: x-grid of the rotated plane as single
*        YST: y-grid of the rotated plane as single
*        ZST: z-grid of the rotated plane as single
*          W: output Variable with the same size as (XST,YST,ZST9) as single
*             allocate and initialized with 0 (in matlab use single(zeros(size(img))) )  
*
*          
* FB
*
*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h" 
#include "tom_dist.h"
#ifdef MATLAB  
	#include "mex.h" 
#endif

/*dist calculation */
void dist (float *theta,int num_theta,float *act_theta,float *psi,float *d,float *xst,float *yst,float *zst,int num_st_x,int num_st_y,float *w)  
{
	

register int i_theta; /* outer loop var*/
long i; 
/* temporary variables */ 
float tmp; /* temporary var */ 
float sin_psi,sin_theta,cos_psi,cos_theta,sin_dp;
float q_tmp_d;


/* ini */
sin_psi=(float)sin(*psi*PI/180);
cos_psi=(float)cos(*psi*PI/180);
tmp=0.0;
*act_theta=*act_theta-1; 


for (i_theta=0; i_theta < num_theta -1  ; i_theta++)
{
	
	sin_theta=(float)sin(theta[i_theta]*PI/180);
	cos_theta=(float)cos(theta[i_theta]*PI/180);
	
	if (i_theta != (int)(*act_theta))
	{
		for (i=0; i < (num_st_x*num_st_y)-1 ; i++ )
		{
			/*  calculate the distance from the plane */
			tmp=(xst[i]*sin_theta*sin_psi)- (yst[i]*sin_theta*cos_psi) + zst[i]*cos_theta; 
			if (fabs(tmp)>*d)
			{
				tmp=*d;
			}
			/* sinc function ==> spaghetti implementation to get some speed  */
			if ( (tmp/(*d)) != 0 )
			{
				q_tmp_d=(tmp/(*d)*PI );
				sin_dp=(float)sin(q_tmp_d);
				tmp=(float) sin_dp / (q_tmp_d) ; 
			}
			else
			{
				tmp=1;
			}	
			w[i]+=tmp;
		}  /*for (i=0; i < (num_st_x_1*num_st_y_1)-1; i++) */
	} /* if (i_theta != (int)(*act_theta))*/
}/* for (i_theta=0; i_theta < num_theta  ; i_theta++) */

}

#ifdef MATLAB 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
const int *num_theta, *dim_st;
const int *tmp_xst,*tmp_yst,*tmp_zst;

/* parse inputs */

/* Check for proper number of arguments */
if (nrhs != 8) 
{ 
	printf("%d",nrhs);
	mexErrMsgTxt("8 input arguments required.\n Syntax: rot3d(input_image,output_image,[Angle(s)],Interpolation_type,[center])");    
}
else 
{
	if (nlhs > 1) 
	{
		mexErrMsgTxt("Too many output arguments.");    
	}
}

/* Check data types */
if (!mxIsSingle(THETA) || !mxIsSingle(ACT_TH) || !mxIsSingle(PSI) || !mxIsSingle(D) || 
    !mxIsSingle(XST) || !mxIsSingle(YST) || !mxIsSingle(ZST) || !mxIsSingle(W)) 
{
	mexErrMsgTxt("Input Variables must be single.\n"); 
}
num_theta=mxGetDimensions(THETA);

/* Check Dimensions */
tmp_xst=mxGetDimensions(XST);
tmp_yst=mxGetDimensions(YST);
tmp_zst=mxGetDimensions(ZST);

if ( (tmp_xst[0]!= tmp_yst[0]) || (tmp_yst[0]!= tmp_zst[0]) || (tmp_xst[0]!= tmp_zst[0]) || 
	 (tmp_xst[1]!= tmp_yst[1]) || (tmp_yst[1]!= tmp_zst[1]) || (tmp_xst[1]!= tmp_zst[1]) )
{
	mexErrMsgTxt("Dimensions of xst yst and zst must agree !.\n");    
}
else
{
	dim_st=mxGetDimensions(XST);
}

/* Do the actual computations in a subroutine */
dist(mxGetData(THETA),num_theta[0],mxGetData(ACT_TH),mxGetData(PSI),mxGetData(D),mxGetData(XST),mxGetData(YST),mxGetData(ZST),dim_st[0],dim_st[1],mxGetData(W));

}

#endif


