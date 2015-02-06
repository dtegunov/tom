/*=================================================================
 *
 * tom_emreadinc.c	reads a file in EM-format
 *
 * The calling syntax is:
 *
 *		[OUT] = tom_emreadinc('Name')
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 20.09.2002
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
    int dims[3];
    char comment[80];
    int emdata[40];
    char dummyb[256];
    unsigned char *bytedata;
    short *intdata;
    float *floatdata;
    FILE *input = 0;
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
    int laufy, ilaufx, ilaufz, ilaufy;
    long int s1, s2, s3;
    int size_area;
    short *pointer;
    long ft;
    int count;
    int lauff,all=0;
    float *testp;
    float waps;
    float testa[11];
    long int xy_dims,fseek_merker;
    unsigned int laufin;
    int mysystem = 1;
    
    /* system determination, little or big endian??? */
    
               if(*(char *)&mysystem == 1)
                        mysystem=1;
                else    mysystem=0;

    /*************************************************/
    
    n = mxGetN(prhs[0])+1;
    infile = mxCalloc(n,sizeof(char));
    mxGetString(prhs[0],infile,n);
    if ((input = fopen (infile, "rb")) == 0)
    {
	mexErrMsgTxt("Could not open file in tom_emreadc.\n");   
    }
    plhs[1] = mxCreateDoubleMatrix(4, 1, mxREAL); 
    fread (magic, 1, 1, input);
    f_magic=magic[0];
        *mxGetPr(plhs[1]) = (double)f_magic; 
    fread (dummya, 1, 2, input);
    f_magic=dummya[0];
        *(mxGetPr(plhs[1])+1) = (double)f_magic; 
    f_magic=dummya[1];
        *(mxGetPr(plhs[1])+2) = (double)f_magic; 
    fread (type, 1, 1, input);
        f_magic=type[0];
        *(mxGetPr(plhs[1])+3) = (double)f_magic;
        
 if (nrhs >1)
 { 
  if (mysystem==1)
  {
   
   
    /* for PC on PC*/
    if ( magic[0]==1 || magic[0]==2 || magic[0]==6)
    {

    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL); 
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);fflush(input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }
        size=dims[0]*dims[1]*dims[2];
        p_nr=mxGetData(prhs[1]);
        nr[0]=p_nr[0];
        nr[1]=p_nr[1];
        nr[2]=p_nr[2];
        p_area=mxGetData(prhs[2]);
        area[0]=p_area[0];
        area[1]=p_area[1];
        area[2]=p_area[2];
        if ((nr[0]+area[0])>dims[0] | nr[1]+area[1]>dims[1] | nr[2]+area[2]>dims[2])
         mexErrMsgTxt("Subregion dimensions larger than volume dimensions.");

        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        area_d[0]=area[0]+1;
        area_d[1]=area[1]+1;
        area_d[2]=area[2]+1;
        
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,area_d,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
        
        switch(type[0])
        {
	    case 1:     printf("only floats !\n");break;
	    case 2:     fseek(input,sizeof(short)*(nr[0]-1),SEEK_CUR);
                    fseek(input,sizeof(short)*(dims[0]*(nr[1]-1)),SEEK_CUR);
                    fseek(input,sizeof(short)*(dims[0]*dims[1]*(nr[2]-1)),SEEK_CUR);
                    ilaufx=0;
                    ilaufz=0;
                    count=0;
                    fseek_merker=0;
                    s1=0;s2=0;s3=0;
                    xy_dims=area_d[0]*area_d[1];
                    for (lauf=nr[2];lauf<=nr[2]+area[2];lauf++)
                    {
                    for (laufy=nr[1];laufy<=nr[1]+area[1];laufy++)
                        {
                            count=fread (&intdata[(ilaufz*xy_dims)+ilaufx*area_d[0]], sizeof(short),area_d[0], input); 
                            fseek(input,(long)(sizeof(short)*(dims[0]-area[0]-1)),SEEK_CUR);
                            fseek_merker=fseek_merker+dims[0];
                            ilaufx=ilaufx+1;
                        }
                        ilaufz=ilaufz+1;
                        ilaufx=0;
                        size_area=0;
                        fseek(input,(long)(sizeof(short)*(dims[1]*dims[0]-fseek_merker)),SEEK_CUR);
                        fseek_merker=0;

                    }

                    break;
	    case 5:    
	                fseek(input,4*(nr[0]-1),SEEK_CUR);
                    fseek(input,4*(dims[0]*(nr[1]-1)),SEEK_CUR);
                    fseek(input,4*(dims[0]*dims[1]*(nr[2]-1)),SEEK_CUR);
                    ilaufx=0;
                    ilaufz=0;
                    count=0;
                    fseek_merker=0;
                    s1=0;s2=0;s3=0;
                    xy_dims=area_d[0]*area_d[1];
                    for (lauf=nr[2];lauf<=nr[2]+area[2];lauf++)
                    {
                    for (laufy=nr[1];laufy<=nr[1]+area[1];laufy++)
                        {
                       
                        count=fread (&floatdata[(ilaufz*xy_dims)+ilaufx*area_d[0]], 4,area_d[0], input); 
                        fseek(input,(long)(sizeof(float)*(dims[0]-area[0]-1)),SEEK_CUR);
                        fseek_merker=fseek_merker+dims[0];
                        ilaufx=ilaufx+1;
                        }
                        ilaufz=ilaufz+1;
                        ilaufx=0;
                        size_area=0;
                        fseek(input,4*(dims[1]*dims[0]-fseek_merker),SEEK_CUR);
                        fseek_merker=0;

                    }

                    break;
	    }
    }
    
    /* for SGI on PC    */
    if (magic[0]==0 || magic[0]==3 || magic[0]==4 || magic[0]==5)
    {
    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        swap(&dims[0],4);
        swap(&dims[1],4);
        swap(&dims[2],4);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL);
    
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        swap(&emdata[lauf],4);
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);fflush(input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }

                size=dims[0]*dims[1]*dims[2];
        p_nr=mxGetData(prhs[1]);
        nr[0]=p_nr[0];
        nr[1]=p_nr[1];
        nr[2]=p_nr[2];
        p_area=mxGetData(prhs[2]);
        area[0]=p_area[0];
        area[1]=p_area[1];
        area[2]=p_area[2];
        if ((nr[0]+area[0])>dims[0] | nr[1]+area[1]>dims[1] | nr[2]+area[2]>dims[2])
         mexErrMsgTxt("Subregion dimensions larger than volume dimensions.");

        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        area_d[0]=area[0]+1;
        area_d[1]=area[1]+1;
        area_d[2]=area[2]+1;
        
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,area_d,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
        
        switch(type[0])
        {
	    case 5:    
	                fseek(input,4*(nr[0]-1),SEEK_CUR);
                    fseek(input,4*(dims[0]*(nr[1]-1)),SEEK_CUR);
                    fseek(input,4*(dims[0]*dims[1]*(nr[2]-1)),SEEK_CUR);
                    ilaufx=0;
                    ilaufz=0;
                    count=0;
                    fseek_merker=0;
                    s1=0;s2=0;s3=0;
                    xy_dims=area_d[0]*area_d[1];
                    for (lauf=nr[2];lauf<=nr[2]+area[2];lauf++)
                    {
                    for (laufy=nr[1];laufy<=nr[1]+area[1];laufy++)
                        {
                        for (laufin=0;laufin<area_d[0];laufin++)
                        {
                        count=fread (&floatdata[(ilaufz*xy_dims)+ilaufx+laufin], 4,1, input); 
                        swap(&floatdata[(ilaufz*xy_dims)+ilaufx+laufin],4);
                        }
                        fseek(input,(long)(sizeof(float)*(dims[0]-area[0]-1)),SEEK_CUR);
                        fseek_merker=fseek_merker+dims[0];
                        ilaufx=ilaufx+area_d[0];
                        }
                        ilaufz=ilaufz+1;
                        ilaufx=0;
                        size_area=0;
                        fseek(input,4*(dims[1]*dims[0]-fseek_merker),SEEK_CUR);
                        fseek_merker=0;

                    }

                    break;
	    }

	     
    }
   }
   else
	   {
	/* for SGI on SGI*/
	if ( magic[0]==3 || magic[0]==4 || magic[0]==5)

    {

    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL); 
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);fflush(input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }
        size=dims[0]*dims[1]*dims[2];
        p_nr=mxGetData(prhs[1]);
        nr[0]=p_nr[0];
        nr[1]=p_nr[1];
        nr[2]=p_nr[2];
        p_area=mxGetData(prhs[2]);
        area[0]=p_area[0];
        area[1]=p_area[1];
        area[2]=p_area[2];
        if ((nr[0]+area[0])>dims[0] | nr[1]+area[1]>dims[1] | nr[2]+area[2]>dims[2])
         mexErrMsgTxt("Subregion dimensions larger than volume dimensions.");

        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        area_d[0]=area[0]+1;
        area_d[1]=area[1]+1;
        area_d[2]=area[2]+1;
        
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,area_d,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetPr(plhs[0]);
                    break;
        }
        
        switch(type[0])
        {
	    case 1:     printf("only floats !\n");break;
	    case 2:     printf("only floats !\n");break;
	    case 5:    
	                fseek(input,4*(nr[0]-1),SEEK_CUR);
                    fseek(input,4*(dims[0]*(nr[1]-1)),SEEK_CUR);
                    fseek(input,4*(dims[0]*dims[1]*(nr[2]-1)),SEEK_CUR);
                    ilaufx=0;
                    ilaufz=0;
                    count=0;
                    fseek_merker=0;
                    s1=0;s2=0;s3=0;
                    xy_dims=area_d[0]*area_d[1];
                    for (lauf=nr[2];lauf<=nr[2]+area[2];lauf++)
                    {
                    for (laufy=nr[1];laufy<=nr[1]+area[1];laufy++)
                        {
                       
                        count=fread (&floatdata[(ilaufz*xy_dims)+ilaufx*area_d[0]], 4,area_d[0], input); 
                        fseek(input,(long)(sizeof(float)*(dims[0]-area[0]-1)),SEEK_CUR);
                        fseek_merker=fseek_merker+dims[0];
                        ilaufx=ilaufx+1;
                        }
                        ilaufz=ilaufz+1;
                        ilaufx=0;
                        size_area=0;
                        fseek(input,4*(dims[1]*dims[0]-fseek_merker),SEEK_CUR);
                        fseek_merker=0;

                    }

                    break;
	    }
    }
    
    /* for PC on SGI    */
    if (magic[0]==1 || magic[0]==2 || magic[0]==6)
    {
    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        swap(&dims[0],4);
        swap(&dims[1],4);
        swap(&dims[2],4);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL);
    
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        swap(&emdata[lauf],4);
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);fflush(input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }

                size=dims[0]*dims[1]*dims[2];
        p_nr=mxGetData(prhs[1]);
        nr[0]=p_nr[0];
        nr[1]=p_nr[1];
        nr[2]=p_nr[2];
        p_area=mxGetData(prhs[2]);
        area[0]=p_area[0];
        area[1]=p_area[1];
        area[2]=p_area[2];
        if ((nr[0]+area[0])>dims[0] | nr[1]+area[1]>dims[1] | nr[2]+area[2]>dims[2])
         mexErrMsgTxt("Subregion dimensions larger than volume dimensions.");

        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        area_d[0]=area[0]+1;
        area_d[1]=area[1]+1;
        area_d[2]=area[2]+1;
        
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,area_d,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,area_d,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
        
        switch(type[0])
        {
	    case 5:    
	                fseek(input,4*(nr[0]-1),SEEK_CUR);
                    fseek(input,4*(dims[0]*(nr[1]-1)),SEEK_CUR);
                    fseek(input,4*(dims[0]*dims[1]*(nr[2]-1)),SEEK_CUR);
                    ilaufx=0;
                    ilaufz=0;
                    count=0;
                    fseek_merker=0;
                    s1=0;s2=0;s3=0;
                    xy_dims=area_d[0]*area_d[1];
                    for (lauf=nr[2];lauf<=nr[2]+area[2];lauf++)
                    {
                    for (laufy=nr[1];laufy<=nr[1]+area[1];laufy++)
                        {
                        for (laufin=0;laufin<area_d[0];laufin++)
                        {
                        count=fread (&floatdata[(ilaufz*xy_dims)+ilaufx+laufin], 4,1, input); 
                        swap(&floatdata[(ilaufz*xy_dims)+ilaufx+laufin],4);
                        }
                        fseek(input,(long)(sizeof(float)*(dims[0]-area[0]-1)),SEEK_CUR);
                        fseek_merker=fseek_merker+dims[0];
                        ilaufx=ilaufx+area_d[0];
                        }
                        ilaufz=ilaufz+1;
                        ilaufx=0;
                        size_area=0;
                        fseek(input,4*(dims[1]*dims[0]-fseek_merker),SEEK_CUR);
                        fseek_merker=0;

                    }

                    break;
	    }

	     
    }
  }
 
 
}
/* normal read of a single layer */
else{

   if (mysystem==1)
  {
   /* for PC on PC */
    if ( magic[0]==1 || magic[0]==2 || magic[0]==6)
    {
    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL); 
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }
        size=dims[0]*dims[1]*dims[2];
        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = (int *) mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
        switch(type[0])
        {
	    case 1:     fread (bytedata, 1, size, input);
	    case 2:     fread (intdata, 2, size, input);
	    case 5:     fread (floatdata, 4, size, input);     
        }
    }
    /* for SGI on PC    */
    if (magic[0]==0 || magic[0]==3 || magic[0]==4 || magic[0]==5)
    {
        plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        swap(&dims[0],4);
        swap(&dims[1],4);
        swap(&dims[2],4);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL);
    
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        swap(&emdata[lauf],4);
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);fflush(input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }
        size=dims[0]*dims[1]*dims[2];
        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
       switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
         
        switch(type[0])
        {
        case 1:
                    for (lauf=0;lauf<dims[0]*dims[1]*dims[2];lauf++)
                        {               
                        fread (&bytedata[lauf], 1, 1, input);
                        swap(&bytedata[lauf],1);
                        }
                    break;
        case 2: 
                      fread (intdata, 2,  size, input);
                      for (lauf=0;lauf<dims[0]*dims[1]*dims[2];lauf++)
                        {
                         swap(&intdata[lauf],2); 
                        }
                      break;
        case 5: 
                    for (lauf=0;lauf<dims[0]*dims[1]*dims[2];lauf++)
                        {
                        fread (&floatdata[lauf], 4,1, input); 
                        swap(&floatdata[lauf],4);
                        }
                    break;

            }
    }
    }

    else
 {
 
 
    /* for SGI on SGI */
    if ( magic[0]==3 || magic[0]==4 || magic[0]==5)

    {
    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL); 
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }
        size=dims[0]*dims[1]*dims[2];
        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = (int *) mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
        switch(type[0])
        {
	    case 1:     fread (bytedata, 1, size, input);
	    case 2:     fread (intdata, 2, size, input);
	    case 5:     fread (floatdata, 4, size, input);     
        }
    }
    /* for PC on SGI    */
     if (magic[0]==0 || magic[0]==1 || magic[0]==2 || magic[0]==6)
   {
        plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        swap(&dims[0],4);
        swap(&dims[1],4);
        swap(&dims[2],4);
        *mxGetPr(plhs[2]) = (double)dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)dims[2]; 
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL);
    
        fread (comment, 1, 80, input);
        for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic; }
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL); 
        fread (emdata, 4, 40, input);
        for (lauf=0;lauf<40;lauf++){
        swap(&emdata[lauf],4);
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic; }
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL); 
        fread (dummyb, 1, 256, input);fflush(input);
        for (lauf=0;lauf<256;lauf++){
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; }
        size=dims[0]*dims[1]*dims[2];
        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
       switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxCHAR_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    bytedata = mxGetData(plhs[0]);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,dims,mxINT16_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    intdata = mxGetData(plhs[0]);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL))==NULL)
                    mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n");    
                    floatdata = mxGetData(plhs[0]);
                    break;
        }
         
        switch(type[0])
        {
        case 1:
                    for (lauf=0;lauf<dims[0]*dims[1]*dims[2];lauf++)
                        {               
                        fread (&bytedata[lauf], 1, 1, input);
                        swap(&bytedata[lauf],1);
                        }
                    break;
        case 2: 
                      fread (intdata, 2,  size, input);
                      for (lauf=0;lauf<dims[0]*dims[1]*dims[2];lauf++)
                        {
                         swap(&intdata[lauf],2); 
                        }
                      break;
        case 5: 
                    for (lauf=0;lauf<dims[0]*dims[1]*dims[2];lauf++)
                        {
                        fread (&floatdata[lauf], 4,1, input); 
                        swap(&floatdata[lauf],4);
                        }
                    break;

            }
    }
 
 
 }
    
 
    
    }
    
    /* for SGI implement swapping asap!!!     if (magic[0]==0 || magic[0]==3 || magic[0]==5)*/
    fclose (input);

/*printf("Single layer"); fflush(stdout);*/


}




