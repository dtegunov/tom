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
 * update: 07/08/06 AK
 *
 *=================================================================*/

#include "io64.h"
#define _FILE_OFFSET_BITS 64
#define _LARGE_FILE_API
#include <stddef.h>
#include <stdio.h> 
#include <stdlib.h> 
#include "mex.h"

#if !( defined(WIN32)  || defined(WIN64) || defined(MACOSX))
#include <sys/vfs.h>
#include <sys/statvfs.h>

int check_required_diskspace (char* filename, unsigned int dims[3], unsigned int varsize) {
 
    struct statvfs64 sStats;

    if (statvfs64(filename, &sStats) == -1) { mexErrMsgTxt("statvfs64() failed in tom_emwritec.\n"); }
    
    if ((sStats.f_bavail * sStats.f_bsize) < (dims[0]*dims[1]*dims[2]*varsize+512)/1024) {
        return 0;
    } else {
        return 1;
    }
    
}

#endif



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


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
    unsigned char magic[1];
    char dummya[2];
    unsigned char type[1];
    unsigned int dims[3], dims2[2];
    char comment[80];
    int emdata[40];
    char dummyb[256];
    float *floatdata;
    FILE *output = 0;
    char *infile;
    unsigned int nr[3];
    double *p_nr, *p_area, *p_varsize, *p_magic, *p_emdata;
    unsigned char *p_comment;
    unsigned int area[3], area_d[3];
    unsigned int n, ndim, size_area;    
    long xy_dims;
    register size_t varsize;
    size_t varsize_data;
    #ifdef WIN32
    long long fseek_merker=0, lauf, laufy=0, ilaufx=0, ilaufz=0, s1=0, s2=0, s3=0, size, swapcounter;
    #elif WIN64
    __int64 fseek_merker=0, lauf, laufy=0, ilaufx=0, ilaufz=0, s1=0, s2=0, s3=0, size, swapcounter;
    #else
    uint64_T fseek_merker=0, lauf, laufy=0, ilaufx=0, ilaufz=0, s1=0, s2=0, s3=0, size, swapcounter;
    #endif
    int mysystem = 1;
    bool swapflag = 0;
    
    /* system determination, little or big endian? */
    /* 1 little endian, 0 big endian*/
    if(*(char *)&mysystem == 1) { mysystem=1;}
    else { mysystem=0; }
    
    
    n = mxGetN(prhs[0])+1;
    infile = mxCalloc(n,sizeof(char));
    mxGetString(prhs[0],infile,n);
    
    
    if (nrhs == 3) {/* create new, empty file */
            
            /*get dimensions and data type of new volume*/
            p_nr=mxGetData(prhs[1]);
            dims[0]=p_nr[0];
            dims[1]=p_nr[1];
            dims[2]=p_nr[2];
            p_varsize = mxGetData(prhs[2]);
            varsize = (size_t) p_varsize[0];
            switch (varsize) {
                case 1: type[0] = 1; break;
                case 2: type[0] = 2; break;
                case 4: type[0] = 5; break;
                case 8: type[0] = 9; break;
            }    
        
            if ((output = fopen (infile, "wb")) == 0) {	mexErrMsgTxt("Could not create file in tom_emwritec.\n"); }

            /*check if disk space if sufficient for new file*/
            #if !( defined(WIN32)  || defined(WIN64) || defined(MACOSX))
            if (check_required_diskspace(infile,dims,varsize) == 0) { fclose(output); mexErrMsgTxt("Not enough disk space in tom_emwritec."); }       
            #endif
            /*write endianess information to file*/
            if (mysystem == 1) { magic[0] = 6; } else { magic[0] = 5;}
            fwrite (magic, 1, 1, output);
            dummya[0]=0;
            dummya[1]=0;
            fwrite (dummya, 1, 2, output);
            
            fwrite (type, 1, 1, output);
            fwrite (dims, 4, 3, output);
            for (lauf=0;lauf<80;lauf++){comment[lauf]=32;}
            fwrite (comment, 1, 80, output);
            for (lauf=0;lauf<40;lauf++){emdata[lauf]=0.0;}
            fwrite (emdata, 4, 40, output);
            for (lauf=0;lauf<256;lauf++){dummyb[lauf]=0.0;}
            fwrite (dummyb, 1, 256, output);
            /*write zeros pagewise to file*/
            if ((floatdata=mxCalloc(dims[1]*dims[0],varsize)) == 0) { fclose(output); mexErrMsgTxt("Memory allocation problem in tom_emwritec.\n"); }
            for (lauf=0;lauf<dims[2];lauf++) { if (fwrite(&floatdata[0], varsize,(size_t) dims[0]*dims[1], output) != dims[0]*dims[1]) { fclose(output); mexErrMsgTxt("Could not write to file in tom_emwritec.\n"); } }
            fflush(output);
            fclose(output);
            mxFree(floatdata);
            return;
            
   }  
    
   else if (nrhs == 7) { /*create new file with given data */

        /*get dimensions and data type of new volume*/
        p_nr=mxGetData(prhs[4]);
        dims[0]=p_nr[0];
        dims[1]=p_nr[1];
        dims[2]=p_nr[2];
        
        p_magic=mxGetData(prhs[3]);
        magic[0] = p_magic[0];
        dummya[0]=0;
        dummya[1]=0;
        type[0] = p_magic[3];

        p_varsize = mxGetData(prhs[2]);
        varsize = (size_t) p_varsize[0];
        
        p_comment = (unsigned char *)mxGetData(prhs[5]);
        p_emdata = mxGetData(prhs[6]);
        floatdata = mxGetData(prhs[1]);
        
        /*open the file*/
        if ((output = fopen (infile, "wb")) == 0) {	mexErrMsgTxt("Could not create file in tom_emwritec.\n"); }
        
         /*check if disk space if sufficient for new file*/
        #if !( defined(WIN32)  || defined(WIN64) || defined(MACOSX))
        if (check_required_diskspace(infile,dims,varsize) == 0) { fclose(output); mexErrMsgTxt("Not enough disk space in tom_emwritec."); }
        #endif

        /*write header to file*/
        fwrite (magic, 1, 1, output);
        fwrite (dummya, 1, 2, output);
        fwrite (type, 1, 1, output);
        fwrite (dims, 4, 3, output);
        
        for (lauf=0;lauf<80;lauf++) { comment[lauf] = p_comment[lauf]; }
        fwrite (comment, 1, 80, output);
        for (lauf=0;lauf<40;lauf++) { emdata[lauf] = p_emdata[lauf]; }
        fwrite (emdata, 4, 40, output);
        for (lauf=0;lauf<256;lauf++){ dummyb[lauf]=0.0; }
        fwrite (dummyb, 1, 256, output);
        /*write data to file*/
        if (fwrite(&floatdata[0], varsize, (size_t) dims[0]*dims[1]*dims[2], output) != dims[0]*dims[1]*dims[2]) { fclose(output); mexErrMsgTxt("Could not write to file in tom_emwritec.\n"); }
        fflush(output);
        fclose(output);
        return;
        
        
   }

 /*open existing file and write to it*/
 else { 
     if ((output = fopen(infile, "r+b")) == 0) { mexErrMsgTxt("Could not open file for reading in tom_emwritec.\n");	}
     fread (magic, 1, 1, output);
     
     /*swap only if system and file endianess differs*/
     if (((magic[0]==1 || magic[0]==2 || magic[0]==6) && mysystem == 1)||((magic[0]==0 || magic[0]==3 || magic[0]==4 || magic[0]==5) && mysystem == 0)) {
         swapflag = 0;
     }
     else {
         swapflag = 1;
     }
     
     fread (dummya, 1, 2, output);
     fread (type, 1, 1, output);
     fread (dims, 4, 3, output);
     if (swapflag == 1) { swap(&dims[0],4); swap(&dims[1],4); swap(&dims[2],4); }
     fread (comment, 1, 80, output);
     fread (emdata, 4, 40, output);
     fread (dummyb, 1, 256, output);
     size=dims[0]*dims[1]*dims[2];
     
     /*determine var size of existing file*/
     switch(type[0]) {
         case 1: varsize = 1;
         break;
         case 2: varsize = 2;
         break;
         case 5: varsize = 4;
         break;
         case 9: varsize = 8;
         break;
     }
     
     /*get data to write from matlab*/
     floatdata=mxGetData(prhs[1]);
         
     /*swap data if necessary*/
     if (swapflag == 1 && varsize > 1) {
         
         switch(type[0]) {
             case 2:
                 for(swapcounter = 0; swapcounter<size; swapcounter++) {
                     swap(&((short*)floatdata)[swapcounter],varsize);
                 }
                 break;
             case 5:
                 for(swapcounter = 0; swapcounter<size; swapcounter++) {
                     swap(&((float*)floatdata)[swapcounter],varsize);
                 }
                 break;
             case 9:
                 for(swapcounter = 0; swapcounter<size; swapcounter++) {
                     swap(&((double*)floatdata)[swapcounter],varsize);
                 }
                 break;
         }
     }
     
     
     p_nr=mxGetData(prhs[2]);
     nr[0]=p_nr[0];
     nr[1]=p_nr[1];
     nr[2]=p_nr[2];
     p_area=mxGetData(prhs[3]);
     area[0]=p_area[0];
     area[1]=p_area[1];
     area[2]=p_area[2];
     if ((nr[0]+area[0])>dims[0]+1 | nr[1]+area[1]>dims[1]+1 | nr[2]+area[2]>dims[2]+1) { fclose(output); mexErrMsgTxt("Subregion dimensions plus offset larger than volume dimensions."); }
     
     if (dims[2]==1) { if (dims[1]==1) ndim=1; else ndim=2; }
     area_d[0]=area[0]+1;
     area_d[1]=area[1]+1;
     area_d[2]=area[2]+1;
     xy_dims=area[0]*area[1];
         
     
     /*seek to beginning of data*/
     #ifdef WIN32
     fseek(output, (long long)varsize*((nr[0]-1)+(dims[0]*(nr[1]-1))+(dims[0]*dims[1]*(nr[2]-1))), SEEK_CUR);
     #elif WIN64
     _fseeki64(output, (__int64)varsize*((nr[0]-1)+(dims[0]*(nr[1]-1))+(dims[0]*dims[1]*(nr[2]-1))), SEEK_CUR);
     #else
     fseeko(output, (off_t)varsize*((nr[0]-1)+(dims[0]*(nr[1]-1))+(dims[0]*dims[1]*(nr[2]-1))), SEEK_CUR);
     #endif
     
     for (lauf=nr[2];lauf<nr[2]+area[2];lauf++) {
         for (laufy=nr[1];laufy<nr[1]+area[1];laufy++) {
             if (fwrite(&floatdata[(ilaufz*xy_dims)+ilaufx*area[0]], varsize, (size_t) area[0], output) != area[0]) { fclose(output); mexErrMsgTxt("Could not write to file in tom_emwritec.\n"); }
             #ifdef WIN32
             fseek(output, (long long)(varsize*(dims[0]-area[0])),SEEK_CUR);
             #elif WIN64
             _fseeki64(output, (__int64)(varsize*(dims[0]-area[0])),SEEK_CUR);
             #else
             fseeko(output, (off_t)(varsize*(dims[0]-area[0])),SEEK_CUR);
             #endif
             fseek_merker=fseek_merker+dims[0];
             ilaufx=ilaufx+1;
         }
         
         ilaufz=ilaufz+1;
         ilaufx=0;
         size_area=0;
         #ifdef WIN32
         fseek(output, (long long) varsize*(dims[1]*dims[0]-fseek_merker), SEEK_CUR);
         #elif WIN64
         _fseeki64(output, (__int64) varsize*(dims[1]*dims[0]-fseek_merker), SEEK_CUR);
         #else
         fseeko(output, (off_t) varsize*(dims[1]*dims[0]-fseek_merker), SEEK_CUR);
         #endif
         fseek_merker=0;
     }
     fflush(output);
     fclose(output);
 }
}

