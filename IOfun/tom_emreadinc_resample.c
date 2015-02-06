/*=================================================================
 *
 * tom_emreadinc_resample.c	reads a file in EM-format and resamples the file
 *
 * The calling syntax is:
 *
 *		[OUT] = tom_emreadinc('Name',binning)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 22.12.2005
 * By: Andreas Korinek
 *
 *=================================================================*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "mex.h"
#include "matrix.h"


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
    unsigned char magic[1];
    int f_magic;
    char dummya[2];
    unsigned char type[1];
    int dims[3];
    int binned_dims[3];
    char comment[80];
    int emdata[40];
    char dummyb[256];
    FILE *input = 0;
    char *infile;
    int n;    
    unsigned short ndim = 3;
    unsigned int lauf, i, ii;
    unsigned short varsize;
    double *p_binning;
    unsigned short binning[3];
    unsigned int lines;
    unsigned long outfile_position = 0;
    mxArray *data;
    void *linebuffer;
    unsigned int pointsperline;
    unsigned int pointcounter = 0;
    unsigned int linesperimage;
    unsigned int linescounter = 0;
    
    /*get binning factors for every dimension*/
    p_binning=mxGetData(prhs[1]);
    binning[0]=(short)p_binning[0];
    binning[1]=(short)p_binning[1];
    binning[2]=(short)p_binning[2];
        
    n = mxGetN(prhs[0])+1;
    infile = mxCalloc(n,sizeof(char));
    mxGetString(prhs[0],infile,n);
    
    if ((input = fopen (infile, "rb")) == 0){
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
   
   /* for PC on PC */
   
    if (magic[0]==1 || magic[0]==2 || magic[0]==6) {
    
        plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL); 
        fread (dims, 4, 3, input);
        
        if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
        
        /* calculate new dimensions */
        binned_dims[0] = dims[0] / binning[0];
        binned_dims[1] = dims[1] / binning[1];
        if (ndim == 3) { binned_dims[2] = dims[2] / binning[2]; } else { binned_dims[2] = 1; }
        
        /* Write binned dimensions to header */
        *mxGetPr(plhs[2]) = (double)binned_dims[0]; 
        *(mxGetPr(plhs[2])+1) = (double)binned_dims[1]; 
        *(mxGetPr(plhs[2])+2) = (double)binned_dims[2]; 
        
        /* Write comment to header */
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
        for (lauf=0;lauf<256;lauf++) {
            f_magic=dummyb[lauf];
            *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf]; 
        }
        
        /*create output arrays and set the type of the data that was found in the file*/
        switch(type[0])
        {
    	case 1:     if ((plhs[0] = mxCreateNumericArray(3,binned_dims,mxCHAR_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
                    varsize = sizeof(char);
                    break;
	    case 2:     if ((plhs[0] = mxCreateNumericArray(3,binned_dims,mxINT16_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
                    varsize = sizeof(short);
                    break;
 	    case 5:	    if ((plhs[0] = mxCreateNumericArray(3,binned_dims,mxSINGLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
                    varsize = sizeof(float);
                    break;
        }

    	    /*allocate output array */
            data =  mxGetData (plhs[0]);
            
            /*number of lines to read */
            if (ndim == 2) { lines = dims[1]; } else {lines = dims[1]*(dims[2]/binning[2]); }
            pointsperline = (int)((float)dims[0]/(float)binning[0]);
            linesperimage = (int)((float)dims[1]/(float)binning[1]);
            
            /* x: dim[1]
             * y: dim[0]
             * z: dim[2]
             */
            
            /*mexPrintf("lines to read: %i, lines per image: %i, points per line: %i, dim0: %i, dim1: %i, dim2: %i\n",lines,linesperimage, pointsperline, dims[0], dims[1], dims[2]); */
            
            /*allocate line buffer */
            linebuffer = mxCalloc (dims[0],varsize);                 
                      
            /* loop over lines */
            for (i=0;i<lines;i++) {

                /*skip pages if 3D image */
                if (ndim == 3 && i%dims[1]==0) {
                    fseek (input, (long)((float)(dims[0]*dims[1]*varsize)*(float)(binning[2]-1)), SEEK_CUR);
                    linescounter = 0;

                }    
                
                if (i%binning[1]==0 && linescounter < linesperimage) {
                    /*read one line */
                    fread (linebuffer, varsize, dims[0], input);
                    pointcounter = 0;
                    
                    /*throw out points */
                    for (ii=0; ii < dims[0]; ii++) {
                        if (ii%binning[0]==0 && pointcounter < pointsperline) {

                            switch(type[0]) {
                                case 1:     ((char*)data)[outfile_position] = ((char*)linebuffer)[ii]; 
                                            break;
                                case 2:     ((short*)data)[outfile_position] = ((short*)linebuffer)[ii];
                                            break;
                                case 5:     ((float*)data)[outfile_position] = ((float*)linebuffer)[ii]; 
                                            break;
                            }
                            outfile_position++;
                            pointcounter++;
                        }
                    }
                    linescounter++;
                }
                /* skip lines */
                else { fseek (input, (long)(dims[0]*varsize), SEEK_CUR); }
            }

            
            /*free line buffer */
            mxFree (linebuffer); 
            
   }
   else {
       mexErrMsgTxt ("Only little endian file format on little endian system supported.\n");   
   }
       
    fclose (input);
   
}
