/***********************************************************************//**
 * \file tom_mex_emwrite_new.c
 * \brief MEX-Function for creating a new em-file.
 * \author  Thomas Haller
 * \version 0.2
 * \date    21.01.2008
 **************************************************************************/

#include <mex.h>
#include <io64.h>

#include <assert.h>
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

    #define __MAX_S__LENGTH__ (4096)
    char s[__MAX_S__LENGTH__+1];

    #define __BUFFERLENGTH_SYNOPSIS__ 1024
    char synopsis[__BUFFERLENGTH_SYNOPSIS__];

    int iotype, element_size;

    size_t k;
    FILE *f;

    tom_io_em_header header;
    uint32_t dims[3];

    const size_t buffer_size = 16384;
    char *buffer;

    memset(&header, 0, sizeof(header));


    snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "%s(filename, size, magic, [comment, emdata, userdata])", mexFunctionName());
    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis);
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ ".\n");
        return;
    }

    if (nlhs>0 || nrhs<3 || nrhs>6) {
        snprintf(s, __MAX_S__LENGTH__, "%s: Wrong number of in/output arguments: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }

    {
        const mxArray *arg;
        mxArray *tmpArray = NULL;
        double *pdouble;
        const int8_t *pint8;
        const int32_t *pint32;


        /* filename */
        arg = prhs[0];
        if (!mxIsChar(arg) || mxGetNumberOfDimensions(arg)!=2 || mxGetM(arg)!=1 || (filename_length=mxGetN(arg))<1 || filename_length>=__MAX_FILENAME_LENGTH__) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs file name as first parameter: %s.", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        mxGetString(arg, filename, __MAX_FILENAME_LENGTH__);

        /* size */
        arg = prhs[1];
        if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2 || mxGetNumberOfElements(arg)!=3) {
            snprintf(s, __MAX_S__LENGTH__, "%s: The size must be a 1x3 vector: %s", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        if (!(tmpArray=getDoubleArray(arg))) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the subregion. Maybe out of memory.", mexFunctionName());
            mexErrMsgTxt(s);
        }
        pdouble = mxGetPr(tmpArray);
        for (k=0; k<3; k++) {
            dims[k] = header.dims[k] = (uint32_t)pdouble[k];
            if (pdouble[k]!=dims[k] || dims[k]<1) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The size must contain positive integer values: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
        }
        mxDestroyArray(tmpArray);



        /* magic */
        arg = prhs[2];
        if (mxIsChar(arg)) {
            mxGetString(arg, s, __MAX_S__LENGTH__);
            header.machine = 6;
            header.byte2 = 0;
            header.byte3 = 0;
            if (strcmp(s, "int8")) {
                header.type = 1;
            } else if (strcmp(s, "int16")) {
                header.type = 2;
            } else if (strcmp(s, "int32")) {
                header.type = 4;
            } else if (strcmp(s, "single") || strcmp(s, "float")) {
                header.type = 5;
            } else if (strcmp(s, "complex32") || strcmp(s, "single_complex") || strcmp(s, "float_complex")) {
                header.type = 8;
            } else if (strcmp(s, "double")) {
                header.type = 9;
            } else if (strcmp(s, "complex64") || strcmp(s, "double_complex")) {
                header.type = TOM_IO_TYPE_COMPLEX64;
            } else {
                snprintf(s, __MAX_S__LENGTH__, "%s: magic must be a string, or a 4-vector of int8 (e.g. possible values are single, double, int32, ...).", mexFunctionName());
                mexErrMsgTxt(s);
            }
        } else {
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2 || mxGetNumberOfElements(arg)!=4 || (mxGetN(arg)!=1&&mxGetM(arg)!=1)) {
                snprintf(s, __MAX_S__LENGTH__, "%s: magic must be a string, or a 4-vector of int8.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            pint8 = (const int8_t *)mxGetData(arg);
            header.machine = pint8[0];
            header.byte2 = pint8[1];
            header.byte3 = pint8[2];
            header.type = pint8[3];
        }



        /* comment */
        if (nrhs >= 4) {
            arg = prhs[3];
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2) {
                snprintf(s, __MAX_S__LENGTH__, "%s: comment must be a 80-vector of int8 or [].", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (mxGetNumberOfElements(arg) > 0) {
                if (mxGetClassID(arg)!=mxINT8_CLASS || mxGetNumberOfElements(arg)!=80 || (mxGetN(arg)!=1&&mxGetM(arg)!=1)) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: comment must be a 80-vector of int8 or [].", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pint8 = (int8_t *)mxGetData(arg);
                for (k=0; k<80; k++) { header.comment[k] = pint8[k]; }
            }
        }

        /* emdata */
        if (nrhs >= 5) {
            arg = prhs[4];
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2) {
                snprintf(s, __MAX_S__LENGTH__, "%s: emdata must be a 40-vector of int32 or [].", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (mxGetNumberOfElements(arg) > 0) {
                if (mxGetClassID(arg)!=mxINT32_CLASS || mxGetNumberOfElements(arg)!=40 || (mxGetN(arg)!=1&&mxGetM(arg)!=1)) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: comment must be a 40-vector of int32 or [].", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pint32 = (int32_t *)mxGetData(arg);
                for (k=0; k<40; k++) { header.emdata[k] = pint32[k]; }
            }
        }


        /* userdata */
        if (nrhs >= 6) {
            arg = prhs[5];
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2) {
                snprintf(s, __MAX_S__LENGTH__, "%s: userdata must be a 256-vector of int8 or [].", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (mxGetNumberOfElements(arg) > 0) {
                if (mxGetClassID(arg)!=mxINT8_CLASS || mxGetNumberOfElements(arg)!=256 || (mxGetN(arg)!=1&&mxGetM(arg)!=1)) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: userdata must be a 256-vector of int8 or [].", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pint8 = (int8_t *)mxGetData(arg);
                for (k=0; k<256; k++) { header.userdata[k] = pint8[k]; }
            }
        }
    }

    assert(sizeof(tom_io_em_header) == 512);

    if (!tom_io_em_is_valid_header(&header, sizeof(header))) {
        snprintf(s, __MAX_S__LENGTH__, "%s: The given header is not a valid EM-header.", mexFunctionName());
        mexErrMsgTxt(s);
    }

    iotype = tom_io_em_get_iotype_header(&header);
    element_size = tom_io_iotype_datasize(iotype);

    if (!(buffer = (char *)calloc(buffer_size, sizeof(char)))) {
        snprintf(s, __MAX_S__LENGTH__, "%s: Error allocating memory for temporary buffer.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }

    if (!(f = fopen(filename, "wb"))) {
        free(buffer);
        snprintf(s, __MAX_S__LENGTH__, "%s: Could not open file \"%s\".", mexFunctionName(), filename);
        mexErrMsgTxt(s);
    }
    if (tom_io_em_is_swaped(&header)) {
        tom_io_em_swap_header(&header);
    }
    if (fwrite(&header, sizeof(header), 1, f) != 1) {
        free(buffer);
        fclose(f);
        snprintf(s, __MAX_S__LENGTH__, "%s: Error writing the header to \"%s\"", mexFunctionName(), filename);
        mexErrMsgTxt(s);
    }

    {
        uint64_t remaining = (uint64_t)dims[0] * (uint64_t)dims[1] * (uint64_t)dims[2] * (uint64_t)element_size;
        k = buffer_size;
        while (remaining > 0) {
            if (remaining < k) {
                k = remaining;
            }
            remaining -= k;
            if (fwrite(buffer, sizeof(buffer[0]), k, f) != k) {
                free(buffer);
                fclose(f);
                snprintf(s, __MAX_S__LENGTH__, "%s: Error writing to file \"%s\"", mexFunctionName(), filename);
                mexErrMsgTxt(s);
            }
        }
    }

    free(buffer);
    fclose(f);
}



