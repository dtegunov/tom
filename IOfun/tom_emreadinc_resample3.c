/*=================================================================
 *
 * tom_emreadinc_resample.c	reads a file in EM-format and resamples the file
 *
 * The calling syntax is:
 *
 *		[OUT] = tom_emreadinc_resample('Name',[resample x y z],[binning x y z],[subregion start x y z],[subregion size x y z])
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 30/08/2006
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
    int resampled_dims[3], output_dims[3];
    char comment[80];
    int emdata[40];
    char dummyb[256];
    FILE *input = 0;
    char *infile;
    unsigned int n, ndim = 3;
    register unsigned int lauf;
    register size_t varsize;
    double *p_nr, *p_area;
    unsigned int nr[3];
    double *p_resampling, *p_binning;
    unsigned int resampling[3], binning[3];
    #ifdef WIN32
    register long long outfile_position = 0;
    #elif WIN64
    register __int64 outfile_position = 0;
    #else
    register uint64_T outfile_position = 0;
    #endif
    mxArray *data;
    void *linebuffer, *binbuffer, *fieldbuffer;
    register unsigned int pagecounter, linecounter, pointcounter, binbuffer_pagecounter = 0;
    register long bufferlauf, binbuffer_position = 0;
    int mysystem = 1;
    bool swapflag;
    float floatbuffer = 0, bin_vol = 0;
    
    /*get resampling factors for every dimension*/
    if (nrhs > 1) {
        p_resampling=mxGetData(prhs[1]);
        resampling[0]=(unsigned int)p_resampling[0];
        resampling[1]=(unsigned int)p_resampling[1];
        resampling[2]=(unsigned int)p_resampling[2];
    }
    else {
        resampling[0] = 1;
        resampling[1] = 1;
        resampling[2] = 1;
    }
    
    /*get binning factors for every dimension*/
    if (nrhs > 2) {
        p_binning=mxGetData(prhs[2]);
        binning[0]=(unsigned int)p_binning[0];
        binning[1]=(unsigned int)p_binning[1];
        binning[2]=(unsigned int)p_binning[2];
    }
    else {
        binning[0] = 1;
        binning[1] = 1;
        binning[2] = 1;
    }

    n = mxGetN(prhs[0])+1;
    infile = mxMalloc(n);
    mxGetString(prhs[0],infile,n);
    
    if ((input = fopen (infile, "rb")) == 0){
        mexErrMsgTxt("Could not open file in tom_emreadc.\n");
    }
    mxFree(infile);
    
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
    if (nrhs > 3) {
        p_nr=mxGetData(prhs[3]);
        nr[0]=p_nr[0]-1;
        nr[1]=p_nr[1]-1;
        nr[2]=p_nr[2]-1;
        p_area=mxGetData(prhs[4]);
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
    
    /*create output array and set the type of the data that was found in the file*/
    output_dims[0] = (unsigned int) resampled_dims[0] / (unsigned int) binning[0];
    output_dims[1] = (unsigned int) resampled_dims[1] / (unsigned int) binning[1];
    if (ndim == 3) { output_dims[2] = (unsigned int) resampled_dims[2] / (unsigned int) binning[2]; } else { output_dims[2] = 1; binning[2] = 1; }
    
    if (output_dims[0] < 1 || output_dims[1] < 1 || output_dims[2] < 1) { mexErrMsgTxt("Empty return matrix in tom_emreadc.\n"); }
    
    if (binning[0] == 1 && binning[1] == 1 && binning[2] == 1) {
    
        switch(type[0]) {
            case 1:     if ((plhs[0] = mxCreateNumericArray(3,output_dims,mxCHAR_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
            varsize = 1;
            break;
            case 2:     if ((plhs[0] = mxCreateNumericArray(3,output_dims,mxINT16_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
            varsize = 2;
            break;
            case 5:	    if ((plhs[0] = mxCreateNumericArray(3,output_dims,mxSINGLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
            varsize = 4;
            break;
            case 9:     if ((plhs[0] = mxCreateNumericArray(3,output_dims,mxDOUBLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
            varsize = 8;
            break;
        }
    }
    else {
        if ((plhs[0] = mxCreateNumericArray(3,output_dims,mxSINGLE_CLASS,mxREAL))==NULL) { mexErrMsgTxt("Memory allocation problem in tom_emreadc.\n"); }
        switch(type[0]) {
            case 1:   varsize = 1;
            break;
            case 2:   varsize = 2;
            break;
            case 5:	  varsize = 4;
            break;
            case 9:   varsize = 8;
            break;
        }
    }
    bin_vol = (float) binning[0] * binning[1] * binning[2];

    /*allocate output array */
    data = mxGetData (plhs[0]);
    
    /*allocate line buffer */
    linebuffer = mxMalloc(resampled_dims[0]*resampling[0]*varsize);
    
    /*allocate variable to hold binning subregion*/
    binbuffer = mxMalloc(resampled_dims[0]*resampling[0]*resampled_dims[1]*resampling[1]*binning[2]*varsize);
    
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
                        case 1:     ((char*)binbuffer)[binbuffer_position] = ((char*)linebuffer)[pointcounter];
                        break;
                        case 2:     ((short*)binbuffer)[binbuffer_position] = ((short*)linebuffer)[pointcounter];
                        break;
                        case 5:     ((float*)binbuffer)[binbuffer_position] = ((float*)linebuffer)[pointcounter];
                        break;
                        case 9:     ((double*)binbuffer)[binbuffer_position] = ((double*)linebuffer)[pointcounter];
                    }
                    
                    binbuffer_position++;
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
            _fseeki64(input, (__int64)(dims[0]*dims[1]*varsize)*(resampling[2]-1), SEEK_CUR); 
            #else
            fseeko(input, (off_t)(dims[0]*dims[1]*varsize)*(resampling[2]-1), SEEK_CUR); 
            #endif
        }

        binbuffer_pagecounter++;
        /*Check if the binning buffer is full*/
        if (binbuffer_pagecounter == binning[2]) {
            
            /*swap values if endianess differs*/
            if ((binning[0] > 1 || binning[1] > 1 || binning[2] > 1) && swapflag == 1 && varsize > 1) {
                register unsigned int swapcounter2;
                switch(type[0]) {
                    case 2:
                        for(swapcounter2 = 0; swapcounter2<resampled_dims[0]*resampling[0]*resampled_dims[1]*resampling[1]*binning[2]; swapcounter2++) {
                            swap(&((short*)binbuffer)[swapcounter2],varsize);
                        }
                        break;
                    case 5:
                        for(swapcounter2 = 0; swapcounter2<resampled_dims[0]*resampling[0]*resampled_dims[1]*resampling[1]*binning[2]; swapcounter2++) {
                            swap(&((float*)binbuffer)[swapcounter2],varsize);
                        }
                        break;
                    case 9:
                        for(swapcounter2 = 0; swapcounter2<resampled_dims[0]*resampling[0]*resampled_dims[1]*resampling[1]*binning[2]; swapcounter2++) {
                            swap(&((double*)binbuffer)[swapcounter2],varsize);
                        }
                        break;
                }
            }
            
            /*bin and copy the values from the binning buffer to the output matrix*/
            unsigned int idx, idy, bidx, bidy, bidz, ixrun = 0, iyrun = 0;
            for(idy=0;idy<resampled_dims[1]-binning[1]+1;idy = idy + binning[1]) {
                for(idx=0;idx<resampled_dims[0]-binning[0]+1;idx = idx + binning[0]) {
                    /*3d data*/
                    if (ndim == 3) {
                        
                        for(bidz=0;bidz<binning[2];bidz++) {
                            for(bidy=0;bidy<binning[1];bidy++) {
                                for(bidx=0;bidx<binning[0];bidx++) {
                                    switch(type[0]) {
                                    case 1:
                                        floatbuffer = ((char*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0] + bidz*resampled_dims[1]*resampled_dims[0]] + floatbuffer;
                                        break;
                                    case 2:
                                        floatbuffer = ((short*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0] + bidz*resampled_dims[1]*resampled_dims[0]] + floatbuffer;
                                        break;
                                    case 5:    
                                        floatbuffer = ((float*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0] + bidz*resampled_dims[1]*resampled_dims[0]] + floatbuffer;
                                        break;
                                    case 9:
                                        floatbuffer = ((double*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0] + bidz*resampled_dims[1]*resampled_dims[0]] + floatbuffer;
                                        break;

                                    }    
                                }
                            }
                        }
                    
                    } 
                    /*2d data*/
                    else {
                        
                        for(bidy=0;bidy<binning[1];bidy++) {
                            for(bidx=0;bidx<binning[0];bidx++) {
                                switch(type[0]) {
                                case 1:
                                    floatbuffer = ((char*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0]] + floatbuffer;
                                    break;
                                case 2:
                                    floatbuffer = ((short*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0]] + floatbuffer;
                                    break;
                                case 5:
                                    floatbuffer = ((float*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0]] + floatbuffer;
                                    break;
                                case 9:
                                    floatbuffer = ((double*)binbuffer)[idx+bidx + (idy+bidy)*resampled_dims[0]] + floatbuffer;
                                    break;
                                }

                            }  
                        }

                    }
                    
                    ((float*)data)[outfile_position] = (float)floatbuffer / bin_vol;
                    floatbuffer = 0.0;
                    ixrun++;
                    outfile_position++;
                }
                ixrun = 0;
                iyrun++;
            }
            iyrun = 0;
            
            /*reset counters*/
            binbuffer_pagecounter = 0;
            binbuffer_position = 0;
        }
        
    }
    
    /*free line buffer */
    mxFree(linebuffer);
    mxFree(binbuffer);
    fclose(input);
    
    /*swap data if necessary*/
    if (swapflag == 1 && varsize > 1 && binning[0] == 1 && binning[1] == 1 && binning[2] == 1) {
        #ifdef WIN32
        register long long swapcounter;
        #elif WIN64
        register __int64 swapcounter;
        #else
        register uint64_T swapcounter;
        #endif
        switch(type[0]) {
            case 2:
                for(swapcounter = 0; swapcounter<output_dims[0]*output_dims[1]*output_dims[2]; swapcounter++) {
                    swap(&((short*)data)[swapcounter],varsize);
                }
                break;
            case 5:
                for(swapcounter = 0; swapcounter<output_dims[0]*output_dims[1]*output_dims[2]; swapcounter++) {
                    swap(&((float*)data)[swapcounter],varsize);
                }
                break;
            case 9:
                for(swapcounter = 0; swapcounter<output_dims[0]*output_dims[1]*output_dims[2]; swapcounter++) {
                    swap(&((double*)data)[swapcounter],varsize);
                }
                break;
        }
    }
    
}
