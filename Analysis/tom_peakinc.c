/*=================================================================
*
* rot3d.c    Performs max
* The syntax is:
*
*        
*
*
* Last changes: Oct. 20, 2005
* fb
*
*=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include "tom_peakinc.h"

#ifdef MATLAB  
    #include "matrix.h"	
    #include "mex.h" 
#endif

/* findint the maximum*/
void find_max (float *stream,int *size_stream ,float *max_value,int *max_position)
{
    int k,i,j,sz,sy,sx;
    long index;
    
    sx=size_stream[0];
    sy=size_stream[1];
    sz=size_stream[2];
    
   *max_value=-99999.0;
    index=0;
    
    for (k=0; k < sz; k++)
	{	
		for (j=0; j < sy; j++)
		{
			for (i=0; i < sx; i++) 
			{
               if (stream[index]>*max_value)
               {
                    max_position[0]=i;
                    max_position[1]=j;
                    max_position[2]=k; 
                   *max_value=stream[index];  
               }
               index++;  
            }
        }
    } 
  }





#ifdef MATLAB 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int *dims;
float  *max;
int *max_pos;

dims=mxGetData(DIMS);
max=mxGetData(MAX);
max_pos=mxGetData(POS_MAX);

find_max(mxGetData(STREAM),dims,max,max_pos);
/* Matlab starts counting at 1!!!*/
max_pos[0]=max_pos[0]+1;
max_pos[1]=max_pos[1]+1;
max_pos[2]=max_pos[2]+1;

}
#endif



