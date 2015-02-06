/*=================================================================
 *
 * tom_emreadinc_resample.c	reads a file in EM-format and resamples the file
 *
 * The calling syntax is:
 *
 *		[OUT] = tom_emreadinc_resample('Name',[resample x y z],[subregion start x y z],[subregion size x y z])
 *
 *
 *  Compile this program with 'mex tom_emreadinc_resample2.c' on linux
 *  and 'mex -DWIN32 tom_emreadinc_resample2.c' on 32 bit Windows
 *
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 29/04/2006
 * By: Andreas Korinek
 *
 *=================================================================*/

#include "io64.h"
#define _FILE_OFFSET_BITS 64
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

void swap (x, size)
unsigned char *x;
int            size;
{
    unsigned char  c;
    unsigned short s;
    unsigned long  l;
    
    switch (size)
    {
        case 2:        /* swap two bytes */
            c = *x;
            *x = *(x+1);
            *(x+1) = c;
            break;
        case 4:        /* swap two shorts (2-byte words) */
            s = *(unsigned short *)x;
            *(unsigned short  *)x = *((unsigned short *)x + 1);
            *((unsigned short *)x + 1) = s;
            swap ((char *)x, 2);
            swap ((char *)((unsigned short *)x+1), 2);
            break;
        case 8:        /* swap two longs (4-bytes words) */
            l = *(unsigned long *)x;
            *(unsigned long  *)x = *((unsigned long *)x + 1);
            *((unsigned long *)x + 1) = l;
            swap ((char *)x, 4);
            swap ((char *)((unsigned long *)x+1), 4);
            break;
    }
    
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    unsigned char magic[1];
    int f_magic;
    char dummya[2];
    unsigned char type[1];
    unsigned int dims[3];
    int resampled_dims[3];
    char comment[80];
    int emdata[40];
    char dummyb[256];
    FILE *input = 0;
    char *infile;
    unsigned int n;
    unsigned int ndim = 3;
    register unsigned int lauf;
    register size_t varsize;
    double *p_nr;
    double *p_area;
    unsigned int nr[3];
    double *p_resampling;
    unsigned int resampling[3];
    #ifdef WIN32
    register long long outfile_position = 0;
    #elif WIN64
    register __int64 outfile_position = 0;
    #else
    register uint64_T outfile_position = 0;
    #endif   
    mxArray *data;
    void *linebuffer;
    register unsigned int pagecounter, linecounter, pointcounter;
    int mysystem = 1;
    bool swapflag;
       
    /*get binning factors for every dimension*/
    if (nrhs > 1) {
        p_resampling=mxGetData(prhs[1]);
        resampling[0]=(unsigned short)p_resampling[0];
        resampling[1]=(unsigned short)p_resampling[1];
        resampling[2]=(unsigned short)p_resampling[2];
    }
    else {
        resampling[0] = 1;
        resampling[1] = 1;
        resampling[2] = 1;
    }
    
    n = mxGetN(prhs[0])+1;
    infile = mxMalloc(n);
    mxGetString(prhs[0],infile,n);
    
    if ((input = fopen (infile, "rb")) == 0){
        mexErrMsgTxt("Could not open file in tom_emreadc.\n");
    }
    
    plhs[1] = mxCreateDoubleMatrix(4, 1, mxREAL);
    fread (magic, 1, 1, input);
    
    /*determine if data in file has to be swapped or not*/
    f_magic=magic[0];
    *mxGetPr(plhs[1]) = (double)f_magic;
    
    /* system determination, little or big endian? */
    /* 1 little endian, 0 big endian*/
    if(*(char *)&mysystem == 1) { mysystem=1;}
    else { mysystem=0; }
    
    /*swap only if system and file endianess differs*/
    if (((magic[0]==1 || magic[0]==2 || magic[0]==6) && mysystem == 1)||((magic[0]==0 || magic[0]==3 || magic[0]==4 || magic[0]==5) && mysystem == 0)) {
        swapflag = 0;
    }
    else {
        swapflag = 1;
    }
    
    fread (dummya, 1, 2, input);
    f_magic=dummya[0];
    *(mxGetPr(plhs[1])+1) = (double)f_magic;
    f_magic=dummya[1];
    *(mxGetPr(plhs[1])+2) = (double)f_magic;
    
    fread (type, 1, 1, input);
    f_magic=type[0];
    *(mxGetPr(plhs[1])+3) = (double)f_magic;
    
    
    
    plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL);
    fread (dims, 4, 3, input);
    if (swapflag == 1) { swap(&dims[0],4); swap(&dims[1],4); swap(&dims[2],4); }
    
    if (dims[2]==1) {if (dims[1]==1) ndim=1; else ndim=2;}
    
    
    /*get subregion data if available*/
    if (nrhs == 4) {
        p_nr=mxGetData(prhs[2]);
        nr[0]=p_nr[0] - 1;
        nr[1]=p_nr[1] - 1;
        nr[2]=p_nr[2] - 1;
        p_area=mxGetData(prhs[3]);
        resampled_dims[0]=p_area[0];
        resampled_dims[1]=p_area[1];
        resampled_dims[2]=p_area[2];
        if (nr[0]+resampled_dims[0]*resampling[0] > dims[0] || nr[1]+resampled_dims[1]*resampling[1] > dims[1] || nr[2]+resampled_dims[2]*resampling[2] > dims[2]) {
            mexErrMsgTxt("Subregion dimensions larger than volume dimensions.\n");
        }
    }
    /*read in full volume*/
    else {
        nr[0] = 0;
        nr[1] = 0;
        nr[2] = 0;
        resampled_dims[0] = (unsigned int)dims[0] / (unsigned int)resampling[0];
        resampled_dims[1] = (unsigned int)dims[1] / (unsigned int)resampling[1];
        if (ndim == 3) { resampled_dims[2] = (unsigned int)dims[2] / (unsigned int)resampling[2]; } else { resampled_dims[2] = 1; }
    }
    
    /* Write binned dimensions to header */
    *mxGetPr(plhs[2]) = (double)resampled_dims[0];
    *(mxGetPr(plhs[2])+1) = (double)resampled_dims[1];
    *(mxGetPr(plhs[2])+2) = (double)resampled_dims[2];
    
    
    /* Write comment to header */
    plhs[3] = mxCreateDoubleMatrix(80, 1, mxREAL);
    fread (comment, 1, 80, input);
    for (lauf=0;lauf<80;lauf++){
        f_magic=comment[lauf];
        *(mxGetPr(plhs[3])+lauf) = (double)f_magic;
    }
    /* write header parameters to output */
    plhs[4] = mxCreateDoubleMatrix(40, 1, mxREAL);
    fread (emdata, 4, 40, input);
    for (lauf=0;lauf<40;lauf++){
        if (swapflag == 1) { swap(&emdata[lauf],4); }
        f_magic=emdata[lauf];
        *(mxGetPr(plhs[4])+lauf) = (double)f_magic;
    }
    /*write fillup to output*/
    plhs[5] = mxCreateDoubleMatrix(256, 1, mxREAL);
    fread (dummyb, 1, 256, input);
    for (lauf=0;lauf<256;lauf++) {
        f_magic=dummyb[lauf];
        *(mxGetPr(plhs[5])+lauf) = (double)dummyb[lauf];
    }
    
    if (resampled_dims[0] < 1 || resampled_dims[1] < 1 || resampled_dims[2] < 1) { mexErrMsgTxt("Empty return matrix in tom_emreadc.\n"); }
    
    /*create output arrays and set the type of the data that was found in the file*/
    switch(type[0]) {
        case 1:     if ((plhs[0] = mxCreateNumericArray(3,resampled_dims,mxCHAR_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
        varsize = 1;
        break;
        case 2:     if ((plhs[0] = mxCreateNumericArray(3,resampled_dims,mxINT16_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
        varsize = 2;
        break;
        case 5:	    if ((plhs[0] = mxCreateNumericArray(3,resampled_dims,mxSINGLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
        varsize = 4;
        break;
        case 9:     if ((plhs[0] = mxCreateNumericArray(3,resampled_dims,mxDOUBLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
        varsize = 8;
        break;
    }
    
    /*allocate output array */
    data = mxGetData (plhs[0]);
    
    /*allocate line buffer */
    linebuffer = mxMalloc(resampled_dims[0]*resampling[0]*varsize);
    
    /*seek to beginning of first page to read*/
    #ifdef WIN32
    fseek(input, (long long)nr[2]*dims[0]*dims[1]*varsize, SEEK_CUR);
    #elif WIN64
    _fseeki64(input, (__int64)nr[2]*dims[0]*dims[1]*varsize, SEEK_CUR);
    #else
    fseeko(input, (off_t)nr[2]*dims[0]*dims[1]*varsize, SEEK_CUR);
    #endif
    /* loop over images */
    for(pagecounter=0;pagecounter<resampled_dims[2];pagecounter++) {
        
        /*seek to beginning of first subregion line in the image*/
        #ifdef WIN32
        fseek(input, (long long)dims[0]*varsize*nr[1], SEEK_CUR);
        #elif WIN64
        _fseeki64(input, (__int64)dims[0]*varsize*nr[1], SEEK_CUR);
        #else
        fseeko(input, (off_t)dims[0]*varsize*nr[1], SEEK_CUR);
        #endif
        /*loop over lines in one image*/
        for(linecounter=0;linecounter<resampled_dims[1]; linecounter++) {
            
            /*read one line, skip first points, skip last points according to subregion definition*/
            #ifdef WIN32
            fseek(input, (long long)nr[0]*varsize,SEEK_CUR);
            #elif WIN64
            _fseeki64(input, (__int64)nr[0]*varsize,SEEK_CUR);
            #else
            fseeko(input, (off_t)nr[0]*varsize,SEEK_CUR);
            #endif
            
            fread(linebuffer, varsize, resampled_dims[0]*resampling[0], input);
            
            #ifdef WIN32
            fseek(input, (long long)(dims[0]-nr[0]-resampled_dims[0]*resampling[0])*varsize,SEEK_CUR);
            #elif WIN64
            _fseeki64(input, (__int64)(dims[0]-nr[0]-resampled_dims[0]*resampling[0])*varsize,SEEK_CUR);
            #else
            fseeko(input, (off_t)(dims[0]-nr[0]-resampled_dims[0]*resampling[0])*varsize,SEEK_CUR);
            #endif
            
            /*throw out points in one line and fill into output array*/
            for(pointcounter=0;pointcounter < resampled_dims[0]*resampling[0];pointcounter++) {
                if (pointcounter % resampling[0] == 0) {
                    
                    switch(type[0]) {
                        case 1:     ((char*)data)[outfile_position] = ((char*)linebuffer)[pointcounter];
                        break;
                        case 2:     ((short*)data)[outfile_position] = ((short*)linebuffer)[pointcounter];
                        break;
                        case 5:     ((float*)data)[outfile_position] = ((float*)linebuffer)[pointcounter];
                        break;
                        case 9:     ((double*)data)[outfile_position] = ((double*)linebuffer)[pointcounter];
                    }
                    
                    outfile_position++;
                }
            }
            
            /*skip the next n lines according to the binning factor*/
            #ifdef WIN32
            fseek(input, (long long)(dims[0]*varsize)*(resampling[1]-1), SEEK_CUR);
            #elif WIN64
            _fseeki64(input, (__int64)(dims[0]*varsize)*(resampling[1]-1), SEEK_CUR);
            #else
            fseeko(input, (off_t)(dims[0]*varsize)*(resampling[1]-1), SEEK_CUR);
            #endif
        }
        /*seek to the end of the page*/
        #ifdef WIN32
        fseek(input, (long long)(dims[1]-nr[1]-resampled_dims[1]*resampling[1])*dims[0]*varsize, SEEK_CUR);
        #elif WIN64
        _fseeki64(input, (__int64)(dims[1]-nr[1]-resampled_dims[1]*resampling[1])*dims[0]*varsize, SEEK_CUR);
        #else
        fseeko(input, (off_t)(dims[1]-nr[1]-resampled_dims[1]*resampling[1])*dims[0]*varsize, SEEK_CUR);
        #endif
        /*skip pages if 3D image*/
        if (ndim == 3) { 
            #ifdef WIN32
            fseek(input, (long long)(dims[0]*dims[1]*varsize)*(resampling[2]-1), SEEK_CUR); 
            #elif WIN64
            _fseeki64(input,(__int64)(dims[0]*dims[1]*varsize)*(resampling[2]-1), SEEK_CUR);
            #else
            fseeko(input, (off_t)(dims[0]*dims[1]*varsize)*(resampling[2]-1), SEEK_CUR); 
            #endif
        }
        
    }
    
    /*free line buffer */
    mxFree(linebuffer);
 
    fclose(input);
    
    /*swap data if necessary*/
    if (swapflag == 1 && varsize > 1) {
        #ifdef WIN32
        register long long swapcounter;
        #elif WIN64
        register __int64 swapcounter;
        #else
        register uint64_T swapcounter;
        #endif
        switch(type[0]) {
            case 2:
                for(swapcounter = 0; swapcounter<outfile_position; swapcounter++) {
                    swap(&((short*)data)[swapcounter],varsize);
                }
                break;
            case 5:
                for(swapcounter = 0; swapcounter<outfile_position; swapcounter++) {
                    swap(&((float*)data)[swapcounter],varsize);
                }
                break;
            case 9:
                for(swapcounter = 0; swapcounter<outfile_position; swapcounter++) {
                    swap(&((double*)data)[swapcounter],varsize);
                }
                break;
        }
    }
    
}
