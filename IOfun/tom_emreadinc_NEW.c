
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "tom_emreadinc_NEW.h"

#ifdef MATLAB
    #include "io64.h"
    #include "mex.h"
#endif

#ifndef MATLAB
    #define _GNU_SOURCE 
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

emHeader* readEmHeader(char* filename){

    FILE *input = 0;
    unsigned char magic[4];
    emHeader* header;
    int mysystem = 1;
    int i=0;

    if ((input = fopen (filename, "rb")) == 0){
        printf("EMREAD ERROR\nCould not open file %c in tom_emreadc.\n",filename);
        return NULL;
    }
    
    header = (emHeader*) calloc(1,sizeof(emHeader));
    
    fread (magic, 1, 4, input);
    
    /*determine if data in file has to be swapped or not*/
    (*header).variableType = magic[3];
    
    /*is the value of complex type?*/
    (*header).complexFlag = (*header).variableType == 8 || (*header).variableType == 10;    
        
    /* system determination, little or big endian? */
    /* 1 little endian, 0 big endian*/
    if(*(char *)&mysystem == 1) { mysystem=1;}
    else { mysystem=0; }
    
    /*swap only if system and file endianess differs*/
    if (((magic[0]==1 || magic[0]==2 || magic[0]==6) && mysystem == 1)||((magic[0]==0 || magic[0]==3 || magic[0]==4 || magic[0]==5) && mysystem == 0)) {
        (*header).swapFlag = 0;
    }
    else {
        (*header).swapFlag = 1;
    }
    
    fread((*header).size, 4, 3, input);
    if ((*header).swapFlag == 1) { swap(&(*header).size[0],4); swap(&(*header).size[1],4); swap(&(*header).size[2],4); }
    
    
    /*read comment of header in chars*/
    fread ((*header).headerChars, 1, 80, input);
    
    /*read user defined parameters into 40 ints*/
    fread ((*header).headerInts, 4, 40, input);
    
    for (i=0;i<40 && (*header).swapFlag == 1;i++)
        swap(&(*header).headerInts[i],4); 
    
    /*set the dimension attribute*/
    (*header).dimension = 3;
    if((*header).size[2] <2)
        (*header).dimension = 2;
    
    fclose(input);
    
    return header;
}

bool checkScaling(int* size,int* resample,int* binning,int dimension){
    bool returnValue;
    
    returnValue = resample[0] + binning[0]<=size[0];
    returnValue = returnValue && resample[1] + binning[1]<=size[1];
    returnValue = (returnValue && resample[2] + binning[2]<=size[2]) || (returnValue && dimension == 2);
    
    if(!returnValue)
        printf("EMREAD ERROR\nResampling or binning parameters are too high.\n");
    
    return returnValue;

}

bool checkSubregion(int* size,int* subregion){

    bool returnValue;
    
    returnValue = subregion[0]>=0;
    returnValue = returnValue && subregion[1]>=0;
    returnValue = returnValue && subregion[2]>=0;
    
    if(!returnValue)
        printf("EMREAD ERROR\nError in subregion. Subregion startpoint is out of bounds.\n");
    
    returnValue = returnValue && (subregion[0] + subregion[3]) <= size[0];
    returnValue = returnValue && (subregion[1] + subregion[4]) <= size[1];
    returnValue = returnValue && (subregion[2] + subregion[5]) <= size[2];    
    
    if(!returnValue)
        printf("EMREAD ERROR\nError in subregion. The subregion dimension are larger than volume dimension. Check subregion, resample and binning.\n");
    
    return returnValue;
    
}

bool checkInputParameters(emHeader* header,int* binning,int* resample,int* subregion){
    
    int realSubregion[6];
    bool returnValue;
    
    /*determine the real size of the unresampled and unbinned volume*/
    realSubregion[0] = subregion[0];
    realSubregion[1] = subregion[1];
    realSubregion[2] = subregion[2];
    
    realSubregion[3] = (int)(((*header).size[0] / resample[0]) / binning[0]);
    realSubregion[4] = (int)(((*header).size[1] / resample[1]) / binning[1]);
    if((*header).dimension == 2)
        realSubregion[5] = 1;
    else
        realSubregion[5] = (int)(((*header).size[2] / resample[2]) / binning[2]);
    
    returnValue = checkSubregion((*header).size,realSubregion);
    returnValue = returnValue && checkScaling((*header).size,resample,binning,(*header).dimension);
    
    return returnValue;
    
}

char returnVarsize(emHeader* header){
/*
    byte            1               1
    short           2               2
    long int        4               4
    float           4               5
    float complex   8               8
    double          8               9
    double complex  16              10
*/    
    char type = (*header).variableType;
    switch(type){
        case(1):{
                return 1;}
        case(2):{
                return 2;}
        case(4):{
                return 4;}
        case(5):{
                return 4;}
        case(8):{
                return 8;}         
        case(9):{
                return 8;}    
        case(10):{
                return 16;}
    }
}

long calculateNumberOfElements(emHeader* header,int* binning,int* resample,int* subregion){

    int elements[3];
    
    if(!(subregion == NULL)){
        return subregion[3]*subregion[4]*subregion[5];
    }
    else{
    
         elements[0] =  (*header).size[0] / (resample[0] + binning[0]);
         elements[1] =  (*header).size[1] / (resample[1] + binning[1]);
         elements[2] =  (*header).size[2] / (resample[2] + binning[2]);
         
         return elements[0]*elements[1]*elements[2];
    }
    
}

emVolume* readEmFile(char* filename,int* binning,int* resample,int* subregion,char* data,char filetype){
/*
 
 
 */
    FILE *input = 0;
    emVolume* volume;
    emHeader* header;
    char varSize;    
    int resampled_dims[3];
    void *linebuffer, *binbuffer;

    volume = (emVolume*) calloc(1,sizeof(emVolume));
    register unsigned int pagecounter, linecounter, pointcounter, binbuffer_pagecounter = 0,binbuffer_position = 0;
    float floatbuffer = 0, bin_vol = 0;
    char* newPointer;
    
    #ifdef WIN32
    register long long outfile_position = 0;
    #elif WIN64
    register __int64 outfile_position = 0;
    #else
    #ifdef MATLAB
    register uint64_T outfile_position = 0;
    #else
    register long long outfile_position = 0;
    #endif
    #endif
    
    subregion[0] -=1; 
    subregion[1] -=1;
    subregion[2] -=1;
    
    (*volume).header = readEmHeader(filename);
    
    if((*volume).header == NULL){
        free(volume);
        return NULL;
    }
    header = (*volume).header;
    
    
    /*check if user parameters are valid options for the em file*/
    if(!checkInputParameters(header,binning,resample,subregion)){
        free(volume);
        printf("EMREAD ERROR\n Skipping emRead.\n");
        return NULL;
    }

    (*header).numberElements = calculateNumberOfElements(header,binning,resample,subregion);  
   
    if ((input = fopen (filename, "rb")) == 0){
        printf("Could not open file %c in tom_emreadc.\n",filename);
    }
    
    /*calculate resampled dims*/
    resampled_dims[0] = (*header).size[0] / resample[0];
    resampled_dims[1] = (*header).size[1] / resample[1];
    if ((*header).dimension == 3) { resampled_dims[2] = (*header).size[2] / resample[2]; } else { resampled_dims[2] = 1; }
    
    
    /*allocate memory*/
    varSize = returnVarsize(header);
    if(data == NULL)
        (*volume).data = malloc(varSize * (*header).numberElements);
    else
        (*volume).data = data;
        
    linebuffer = malloc(subregion[3]*resample[0]*varSize);
    binbuffer = malloc(resampled_dims[0]*resample[0]*resampled_dims[1]*resample[1]*binning[2]*varSize);
    bin_vol = (float) binning[0] * binning[1] * binning[2]; 
    
    /*seek to the first page (in z)*/
    #ifdef WIN32
    fseek(input, (long long)512+ (subregion[2]*(*header).size[0]*(*header).size[1]) * varSize, SEEK_SET);
    #elif WIN64
    _fseeki64(input, (__int64)512+ (subregion[2]*(*header).size[0]*(*header).size[1]) * varSize, SEEK_SET);
    #else
    fseeko(input, (off_t)512+ (subregion[2]*(*header).size[0]*(*header).size[1]) * varSize, SEEK_SET);
    #endif 
    
    /* loop over pages */
    for(pagecounter=0;pagecounter<resampled_dims[2];pagecounter++) {
        
        #ifdef WIN32
        fseek(input, (long long) ((*header).size[0]*subregion[1])*varSize, SEEK_CUR);
        #elif WIN64
        _fseeki64(input, (__int64) ((*header).size[0]*subregion[1])*varSize, SEEK_CUR);
        #else
        fseeko(input, (off_t) ((*header).size[0]*subregion[1])*varSize, SEEK_CUR);
        #endif
        
        /*loop over lines in one image*/
        for(linecounter=0;linecounter<resampled_dims[1]; linecounter++) { 
        
            #ifdef WIN32
            fseek(input, (long long) (subregion[0])*varSize,SEEK_CUR);
            #elif WIN64
            _fseeki64(input, (__int64) (subregion[0])*varSize,SEEK_CUR);
            #else
            fseeko(input, (off_t) (subregion[0])*varSize,SEEK_CUR);
            #endif
            
            fread(linebuffer, varSize, resampled_dims[0]*resample[0], input);
            /*hier vielleicht fehler, nochmal binnen?*/            
            /*printf("%f \n",(*(float*)linebuffer));*/

            #ifdef WIN32
            fseek(input, (long long)((*header).size[0]-(subregion[0])-resampled_dims[0]*resample0])*varSize,SEEK_CUR);
            #elif WIN64
            _fseeki64(input, (__int64)((*header).size[0]-(subregion[0])-resampled_dims[0]*resample0])*varSize,SEEK_CUR);
            #else
            fseeko(input, (off_t)((*header).size[0]-(subregion[0])-resampled_dims[0]*resample[0])*varSize,SEEK_CUR);
            /*printf("hier %d \n %d \n",((*header).size[0]-(subregion[0])-resampled_dims[0]*resample[0])*varSize,resample[0]);*/
            #endif
            
            /*throw out points in one line and fill into output array*/
            for(pointcounter=0;pointcounter < resampled_dims[0]*resample[0];pointcounter++) {
                if (pointcounter % resample[0] == 0) {
                    
                    switch((*header).variableType) {
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
            /*skip the next n lines according to the resample factor*/
            #ifdef WIN32
            fseek(input, (long long)((*header).size[0]*varSize)*(resample[1]-1), SEEK_CUR);
            #elif WIN64
            _fseeki64(input, (__int64)((*header).size[0]*varSize)*(resample[1]-1), SEEK_CUR);
            #else
            fseeko(input, (off_t)((*header).size[0]*varSize)*(resample[1]-1), SEEK_CUR);
/*            printf("hier2 %d \n ",((*header).size[0]*varSize)*(resample[1]-1));*/
            #endif
        }
        
        #ifdef WIN32
        fseek(input, (long long)((*header).size[1]-(subregion[1])-resampled_dims[1]*resample[1])*(*header).size[0]*varSize, SEEK_CUR);
        #elif WIN64
        _fseeki64(input, (__int64)((*header).size[1]-(subregion[1])-resampled_dims[1]*resample[1])*(*header).size[0]*varSize, SEEK_CUR);
        #else
        fseeko(input, (off_t)((*header).size[1]-(subregion[1])-resampled_dims[1]*resample[1])*(*header).size[0]*varSize, SEEK_CUR);
        #endif
        
        if ((*header).dimension == 3) { 
            #ifdef WIN32
            fseek(input, (long long)((*header).size[0]*(*header).size[1]*varSize)*(resample[2]), SEEK_CUR); 
            #elif WIN64
            _fseeki64(input, (__int64)((*header).size[0]*(*header).size[1]*varSize)*(resample[2]), SEEK_CUR); 
            #else
            fseeko(input, (off_t)((*header).size[0]*(*header).size[1]*varSize)*(resample[2]), SEEK_CUR); 
            #endif
        }
        
        binbuffer_pagecounter++;
        
        if (binbuffer_pagecounter == binning[2]) {
            
            /*swap values if endianess differs*/
            if ((binning[0] > 1 || binning[1] > 1 || binning[2] > 1) && (*header).swapFlag == 1 && varSize > 1) {
                register unsigned int swapcounter2;
                switch((*header).variableType) {
                    case 2:
                        for(swapcounter2 = 0; swapcounter2<resampled_dims[0]*resample[0]*resampled_dims[1]*resample[1]*binning[2]; swapcounter2++) {
                            swap(&((short*)binbuffer)[swapcounter2],varSize);
                        }
                        break;
                    case 5:
                        for(swapcounter2 = 0; swapcounter2<resampled_dims[0]*resample[0]*resampled_dims[1]*resample[1]*binning[2]; swapcounter2++) {
                            swap(&((float*)binbuffer)[swapcounter2],varSize);
                        }
                        break;
                    case 9:
                        for(swapcounter2 = 0; swapcounter2<resampled_dims[0]*resample[0]*resampled_dims[1]*resample[1]*binning[2]; swapcounter2++) {
                            swap(&((double*)binbuffer)[swapcounter2],varSize);
                        }
                        break;
                }
            }
        
            /*bin and copy the values from the binning buffer to the output matrix*/
            unsigned int idx, idy, bidx, bidy, bidz, ixrun = 0, iyrun = 0;
            for(idy=0;idy<resampled_dims[1]-binning[1]+1;idy = idy + binning[1]) {
                for(idx=0;idx<resampled_dims[0]-binning[0]+1;idx = idx + binning[0]) {
                    /*3d data*/
                    if ((*header).dimension == 3) {
                        
                        for(bidz=0;bidz<binning[2];bidz++) {
                            for(bidy=0;bidy<binning[1];bidy++) {
                                for(bidx=0;bidx<binning[0];bidx++) {
                                    switch((*header).variableType) {
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
                                switch((*header).variableType) {
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
            
                    
                    ((*volume).data)[outfile_position*varSize] = floatbuffer / bin_vol;
                    /*printf("%f %f\n",(float)(*volume).data[outfile_position*varSize],bin_vol);*/
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
    
    free(linebuffer);
    free(binbuffer);
    fclose(input);
    
    /*swap data if necessary*/
    if ((*header).swapFlag == 1 && varSize > 1 && binning[0] == 1 && binning[1] == 1 && binning[2] == 1 ) {

        #ifdef WIN32
        register long long swapcounter;
        #elif WIN64
        register __int64 swapcounter;
        #else
        register long long swapcounter;
        #endif
        switch((*header).variableType) {
            case 2:
                for(swapcounter = 0; swapcounter<(*header).numberElements; swapcounter++) {
                    swap(&((short*)(*volume).data)[swapcounter],varSize);
                }
                break;
            case 5:
                for(swapcounter = 0; swapcounter<(*header).numberElements; swapcounter++) {
                    swap(&((float*)(*volume).data)[swapcounter],varSize);
                }
                break;
            case 9:
                for(swapcounter = 0; swapcounter<(*header).numberElements; swapcounter++) {
                    swap(&((double*)(*volume).data)[swapcounter],varSize);
                }
                break;
        }
    }
    return volume;
}


#ifndef MATLAB
void main(int argc, char** argv)
{
    double *p_resampling, *p_binning;
    unsigned int resampling[3], binning[3];
    unsigned int n, ndim = 3;
    char *infile;
    emVolume* volume; 
    emHeader* header;
    int tmp[3] = {1,1,1};
    int tmp2[6]={1,1,1,48,48,1};
    /*int tmp2[6]={1,1,1,3,3,3};*/
    int i;
    

    volume = readEmFile("./test.em",tmp,tmp,tmp2,0);
   /* for(i=0;i<9;i++)
        printf("%f \n",(float)(*volume).data[i*4]);
    */

    header = (*volume).header;
    

}
#endif

#ifdef MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    double *p_resampling, *p_binning;
    unsigned int resampling[3], binning[3];
    unsigned int n, ndim = 3;
    char *infile;
    emVolume* volume; 
    emHeader* header;
    int tmp[3] = {1,1,1};
    int tmp1[3] = {2,2,2};
    /*int tmp2[6]={1,1,1,48,48,48};*/
    int tmp2[6]={1,1,1,3,3,3};
    int i;
    const mwSize s[3] ={3,3,1};
    
    plhs[0] = mxCreateNumericArray(3, s,mxSINGLE_CLASS,mxREAL);
    /*printf("%d \n",mxGetData(plhs[0])); */
    volume = readEmFile("./small.em",tmp,tmp,tmp2,(char*) mxGetData(plhs[0]),0);
    for(i=0;i<3*3;i++)
        printf("%f \n",(float)(*volume).data[i*4]);
    
    /*printf("%d \n",(*volume).data);*/

  
}
#endif

