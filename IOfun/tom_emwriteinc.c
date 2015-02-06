/*=================================================================
 *
 * tom_emwriteinc.c	writes a file in EM-format
 *
 * The calling syntax is:
 *
 *		tom_emwriteinc('Name',out,nr,nrarea)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 09.12.2002
 * By: Stephan Nickell
 * Revision: 1.00 by 
 *
 *=================================================================*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "mex.h"
#include "matrix.h"

void swap(char *x, char size)
{
  unsigned char c;
  unsigned short s;
  unsigned long l;

  switch (size)
  {
    case 2: /* swap two bytes */
      c = *x;
      *x = *(x+1);
      *(x+1) = c;
      break;
    case 4: /* swap two shorts (2-byte words) */
      s = *(unsigned short *)x;
      *(unsigned short *)x = *((unsigned short *)x + 1);
      *((unsigned short *)x + 1) = s;
      swap ((char *)x, 2);
      swap ((char *)((unsigned short *)x+1), 2);
      break;
    case 8: /* swap two longs (4-bytes words) */
      l = *(unsigned long *)x;
      *(unsigned long *)x = *((unsigned long *)x + 1);
      *((unsigned long *)x + 1) = l;
      swap ((char *)x, 4);
      swap ((char *)((unsigned long *)x+1), 4);
      break;
  }
}
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
    unsigned char magic[1];
    int f_magic;
    char dummya[2];
    unsigned char type[1];
    int dims[3],dims2[2];
    char comment[80];
    int emdata[40];
    char dummyb[256];
    float *floatdata;
    FILE *output = 0;
    long size;
    char *infile;
    int nr[3];
    double *p_nr;
    int area[3];
    int area_d[3];
    double *p_area;
    int n;    
    int ndim;
    unsigned int lauf;
    int laufy, ilaufx, ilaufz;
    long int s1, s2, s3;
    int size_area;
    int count;
    int all=0;
    long int xy_dims,fseek_merker;
    int mysystem = 1;
    
    /* system determination, little or big endian??? */
    
               if(*(char *)&mysystem == 1)
                        mysystem=1;
                else    mysystem=0;

    /*************************************************/
    n = mxGetN(prhs[0])+1;
    infile = mxCalloc(n,sizeof(char));
    mxGetString(prhs[0],infile,n);
     
 if (nrhs == 2) /* create new, empty file */
 { 
  if (mysystem==1)
  {
    if ((output = fopen (infile, "wb")) == 0)
    {
	mexErrMsgTxt("Could not create file in tom_emwritec.\n");		
    }
    magic[0]=6; /* for PC */
    fwrite (magic, 1, 1, output);
    dummya[0]=0; 
    dummya[1]=0; 
    fwrite (dummya, 1, 2, output);
    type[0]=5; /* for float */
    fwrite (type, 1, 1, output);
    p_nr=mxGetData(prhs[1]);
    dims[0]=p_nr[0];
    dims[1]=p_nr[1];
    dims[2]=p_nr[2];
    fwrite (dims, 4, 3, output);
    for (lauf=0;lauf<80;lauf++){comment[lauf]=0.0;}
    fwrite (comment, 1, 80, output);
    for (lauf=0;lauf<40;lauf++){emdata[lauf]=0.0;}
    fwrite (emdata, 4, 40, output);
    for (lauf=0;lauf<256;lauf++){dummyb[lauf]=0.0;}
    fwrite (dummyb, 1, 256, output);
    dims2[0]=dims[0];
    dims2[1]=dims[1];
/*    if ((floatdata = mxCreateNumericArray(2,dims2,mxSINGLE_CLASS,mxREAL))==NULL) */
    if ((floatdata=(float *)(malloc(dims[1]*dims[0]*4))) == 0)
                    mexErrMsgTxt("Memory allocation problem in tom_emwritec.\n");
    for (lauf=0;lauf<dims[1]*dims[0];lauf++){floatdata[lauf]=0.0;};
    for (lauf=0;lauf<dims[2];lauf++){fwrite(&floatdata[0],4,dims[0]*dims[1],output);};
    fflush(output);
    fclose(output);
    free(floatdata);
    return;
                    
  }
 }
    
    if ((output = fopen (infile, "r+b")) == 0)
    {
	mexErrMsgTxt("Could not open file for reading in tom_emwritec.\n");	
    }
    fread (magic, 1, 1, output);
    fread (dummya, 1, 2, output);
    fread (type, 1, 1, output);

 if (nrhs >2)
 { 
  if (mysystem==1)
  {
   
   
    /* for PC on PC*/
    if (magic[0]==1 || magic[0]==2 || magic[0]==6)
    {
        fread (dims, 4, 3, output);
        fread (comment, 1, 80, output);
        fread (emdata, 4, 40, output);
        fread (dummyb, 1, 256, output);fflush(output);
        size=dims[0]*dims[1]*dims[2];
        floatdata=mxGetData(prhs[1]);
        p_nr=mxGetData(prhs[2]);
        nr[0]=p_nr[0];
        nr[1]=p_nr[1];
        nr[2]=p_nr[2];
        p_area=mxGetData(prhs[3]);
        area[0]=p_area[0];
        area[1]=p_area[1];
        area[2]=p_area[2];
        if ((nr[0]+area[0])>dims[0]+1 | nr[1]+area[1]>dims[1]+1 | nr[2]+area[2]>dims[2]+1)
         {mexErrMsgTxt("Subregion dimensions plus offset larger than volume dimensions."); return;}

        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        area_d[0]=area[0]+1;
        area_d[1]=area[1]+1;
        area_d[2]=area[2]+1;
        switch(type[0])
        {
	    case 1:     printf("only floats !\n");break;
	    case 2:     printf("only floats !\n");break;
	    case 5:    
	                fseek(output,4*(nr[0]-1),SEEK_CUR);
                    fseek(output,4*(dims[0]*(nr[1]-1)),SEEK_CUR);
                    fseek(output,4*(dims[0]*dims[1]*(nr[2]-1)),SEEK_CUR);
                    ilaufx=0;
                    ilaufz=0;
                    count=0;
                    fseek_merker=0;
                    s1=0;s2=0;s3=0;
                    xy_dims=area[0]*area[1];
                    for (lauf=nr[2];lauf<nr[2]+area[2];lauf++)
                    {
                    for (laufy=nr[1];laufy<nr[1]+area[1];laufy++)
                        {
                       
                        count=fwrite (&floatdata[(ilaufz*xy_dims)+ilaufx*area[0]], 4,area[0], output);
/*        printf("count: %i \n",count);fflush(stdout);
        printf("counter: %i \n",(ilaufz*xy_dims)+ilaufx*area[0]);fflush(stdout);
        printf("Data: %f \n",floatdata[(ilaufz*xy_dims)+ilaufx*area[0]]);fflush(stdout);
*/
                        fseek(output,(long)(sizeof(float)*(dims[0]-area[0])),SEEK_CUR);
                        fseek_merker=fseek_merker+dims[0];
                        ilaufx=ilaufx+1;

                        }
                        ilaufz=ilaufz+1;
                        ilaufx=0;
                        size_area=0;
                        fseek(output,4*(dims[1]*dims[0]-fseek_merker),SEEK_CUR);
                        fseek_merker=0;

                    }

                    break;
	    }

	    }
	    }
    fclose (output);
    /*free(floatdata);*/
 }
    
}

