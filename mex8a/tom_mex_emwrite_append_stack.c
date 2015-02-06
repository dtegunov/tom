/***********************************************************************//**
 * \file tom_mex_emwrite_append_stack.c
 * \brief MEX-Function for tom_io_em_write_append_stack
 * \author  Thomas Haller
 * \version 0.1
 * \date    25.11.2007
 **************************************************************************/




#include <mex.h>
#include <io64.h>


#include <math.h>
#include <string.h>

#include <tom/io/io.h>
#include <tom/tools/mex_helpfcn.h>





#ifdef _WIN32
    #define snprintf _snprintf
#endif



/***********************************************************************//**
 * \brief mexFunction to append a stack to an existing em-file.
 **************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {


    size_t filename_length = 0;
    #define __MAX_FILENAME_LENGTH__ ((int)2048)
    char filename[__MAX_FILENAME_LENGTH__+1];
    #define __MAX_S__LENGTH__ (__MAX_FILENAME_LENGTH__+1024)
    char s[__MAX_S__LENGTH__+1];

    tom_io_em_header header;
    #define __BUFFERLENGTH_SYNOPSIS__ 1024
    char synopsis[__BUFFERLENGTH_SYNOPSIS__];
    void *data;

    mxArray *mxData;
    size_t sizex, sizey, sizez;
    mxClassID mxType;
    mxComplexity mxIsComplexVolume;

    size_t i;
    int res;
    int header_read;
    int allow_conversion = 1;

    mxArray *plhs_tmp[5] = { NULL, NULL, NULL, NULL, NULL };





    snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "[size, magic, comment, emdata, userdata] = %s(filename, volume, [allow_conversion])", mexFunctionName());

    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis);
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ ".\n");
        return;
    }

    if (nrhs < 2 || nrhs > 3) {
        snprintf(s, __MAX_S__LENGTH__, "%s: call function with up to 3 parameters: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }
    if (nlhs>5) {
        snprintf(s, __MAX_S__LENGTH__, "%s: Too many output parameters: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }


    {
        const mxArray *PRHS_FILENAME = prhs[0];
        const mxArray *PRHS_VOLUME = prhs[1];
        const mwSize *size;


        /* Check input parameters! */
        if (!mxIsChar(PRHS_FILENAME) || mxGetNumberOfDimensions(PRHS_FILENAME)!=2 || mxGetM(PRHS_FILENAME)!=1 || (filename_length=mxGetN(PRHS_FILENAME))<1) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs the file name as first parameter: %s.", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        if (filename_length >= __MAX_FILENAME_LENGTH__) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Maximal length of file name exceeded. (Recompile with larger buffer :)", mexFunctionName());
            mexErrMsgTxt(s);
        }
        mxGetString(PRHS_FILENAME, filename, __MAX_FILENAME_LENGTH__);

        if (!mxIsNumeric(PRHS_VOLUME) || mxGetNumberOfDimensions(PRHS_VOLUME)>3) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs a numerical volume as second parameter: %s", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        data = mxGetData(PRHS_VOLUME);
        size = mxGetDimensions(PRHS_VOLUME);
        mxType = mxGetClassID(PRHS_VOLUME);
        mxIsComplexVolume = mxIsComplex(PRHS_VOLUME);
        sizex = size[0];
        sizey = size[1];
        sizez = mxGetNumberOfDimensions(PRHS_VOLUME)==3 ? size[2] : 1;

        if (mxIsComplexVolume) {
            if (mxType==mxSINGLE_CLASS || mxType==mxDOUBLE_CLASS) {
                mwSize size[3];
                size_t x, y, z;

                size[0] = sizex*2;
                size[1] = sizey;
                size[2] = sizez;
                if (!(mxData = mxCreateNumericArray(sizez==1 ? 2 : 3, size, mxType, mxREAL))) {
                    mexErrMsgTxt("%s: Error allocating temporary buffer for complex data.");
                }
                data = mxGetData(mxData);
                if (mxType == mxSINGLE_CLASS) {
                    float *data_as_real = (float *)data;
                    const float *data_real = (const float *)mxGetData(PRHS_VOLUME);
                    const float *data_complex = (const float *)mxGetImagData(PRHS_VOLUME);
                    for (z=0; z<sizez; z++) {
                        for (y=0; y<sizey; y++) {
                            for (x=0; x<sizex; x++) {
                                *data_as_real++ = *data_real++;
                                *data_as_real++ = *data_complex++;
                            }
                        }
                    }
                } else {
                    double *data_as_real = (double *)data;
                    const double *data_real = (const double *)mxGetData(PRHS_VOLUME);
                    const double *data_complex = (const double *)mxGetImagData(PRHS_VOLUME);
                    for (z=0; z<sizez; z++) {
                        for (y=0; y<sizey; y++) {
                            for (x=0; x<sizex; x++) {
                                *data_as_real++ = *data_real++;
                                *data_as_real++ = *data_complex++;
                            }
                        }
                    }
                }
            } else {
                snprintf(s, __MAX_S__LENGTH__, "%s: Complex data for this type not supported (currently only for single and double).", mexFunctionName());
                mexErrMsgTxt(s);
            }
        }

        if (nrhs >= 3) {
            mwSize numel;
            numel = mxGetM(prhs[2]) * mxGetN(prhs[2]);
            if (!(mxIsNumeric(prhs[2]) || mxIsLogical(prhs[2])) || numel>1) {
                snprintf(s, __MAX_S__LENGTH__, "%s: allow_conversion must be one of true, false or [].", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (numel == 1) {
                allow_conversion = mxGetScalar(prhs[2]) != 0.;
            }
        }
    }


    {
        /* Allocate the memory for the ouput, so that in case of successfully writing, no
           error can happen afterwards in the mexfunction. */
        switch (nlhs) {
            case 5:
                if (!(plhs_tmp[4] = mxCreateNumericMatrix(1, 256, mxINT8_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
            case 4:
                if (!(plhs_tmp[3] = mxCreateNumericMatrix(1, 40, mxINT32_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
            case 3:
                if (!(plhs_tmp[2] = mxCreateNumericMatrix(1, 80, mxINT8_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
            case 2:
                if (!(plhs_tmp[1] = mxCreateNumericMatrix(1, 4, mxINT8_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
            case 1:
                if (!(plhs_tmp[0] = mxCreateNumericMatrix(3, 1, mxUINT32_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
        }
    }


    res = tom_io_em_write_append_stack(filename, data, getIOTypeFromMxClassID(mxType, mxIsComplexVolume), sizex, sizey, sizez, 0, 0, 0, &header, &header_read, allow_conversion);

    if (res != TOM_ERR_OK) {
        if (res ==TOM_ERR_WRITE_FILE) {
            snprintf(s, __MAX_S__LENGTH__, "%s: IO-error occured. The em-file may now be damaged :(", mexFunctionName());
            mexErrMsgTxt(s);
        } else if (res==TOM_ERR_OPEN_FILE) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Error opening the file \"%s\" for writing.", mexFunctionName(), filename);
            mexErrMsgTxt(s);
        } else if (header_read && res==TOM_ERR_WRONG_IOTYPE_CONVERSION) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Wrong typeconversion: can not convert volume of type %d to type %d as in the em-file.", mexFunctionName(), getIOTypeFromMxClassID(mxType, mxIsComplexVolume), tom_io_em_get_iotype_header(&header));
            mexErrMsgTxt(s);
        } else if (res == TOM_ERR_IOTYPE_NOT_SUPPORTED) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Saving data of type %d to em-file is (currently) not supported.", mexFunctionName(), getIOTypeFromMxClassID(mxType, mxIsComplexVolume));
            mexErrMsgTxt(s);
        } else if (header_read && TOM_ERR_WRONG_DATA_SIZE && (sizex!=header.dims[0] || sizey!=header.dims[1])) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Size mismatch: The volume has dimension %dx%dx%d, but the em-file has size %dx%dx%d.", mexFunctionName(), sizex, sizey, sizez, header.dims[0], header.dims[1], header.dims[2]);
            mexErrMsgTxt(s);
        } else {
            snprintf(s, __MAX_S__LENGTH__, "%s: Unexpected error. Is the \"%s\" a valid em-file? (%d)", mexFunctionName(), filename, res);
            mexErrMsgTxt(s);
        }
    }


    {
        /* Construct the output. */
        void *pdata;
        switch (nlhs) {
            case 5:
                pdata = mxGetData(plhs[4] = plhs_tmp[4]);
                for (i=0; i<256; i++) { ((int8_t *)pdata)[i] = header.userdata[i]; }
            case 4:
                pdata = mxGetData(plhs[3] = plhs_tmp[3]);
                for (i=0; i<40; i++) { ((int32_t *)pdata)[i] = header.emdata[i]; }
            case 3:
                pdata = mxGetData(plhs[2] = plhs_tmp[2]);
                for (i=0; i<80; i++) { ((int8_t *)pdata)[i] = header.comment[i]; }
            case 2:
                pdata = mxGetData(plhs[1] = plhs_tmp[1]);
                ((int8_t *)pdata)[0] = header.machine;
                ((int8_t *)pdata)[1] = header.byte2;
                ((int8_t *)pdata)[2] = header.byte3;
                ((int8_t *)pdata)[3] = header.type;
            case 1:
                pdata = mxGetData(plhs[0] = plhs_tmp[0]);
                ((uint32_t *)pdata)[0] = header.dims[0];
                ((uint32_t *)pdata)[1] = header.dims[1];
                ((uint32_t *)pdata)[2] = header.dims[2];
                break;
            case 0:
            default:
                break;
        }
    }
}



