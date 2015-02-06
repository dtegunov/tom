/***********************************************************************//**
 * \file tom_mex_emwrite.c
 * \brief MEX-Function for writing an em-file.
 * \author  Thomas Haller
 * \version 0.2
 * \date    27.11.2008
 *
 *
 **************************************************************************/


#include <mex.h>
#include <io64.h>


#include <math.h>
#include <errno.h>
#include <string.h>

#include <tom/io/io.h>
#include <tom/tools/mex_helpfcn.h>



#if defined _WIN32 || defined _WIN64
    #include <system/pstdint.h>
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

    size_t i, j, k;

    const void *w_data;
    int w_iotype;

    size_t dims[3];
    uint32_t *subregion = NULL;
    uint32_t subregion_field[6];
    tom_io_em_header header;

    void *complexCopy = 0;
    char autoval_magic = 1;
    size_t stridey=0, stridez=0;

    memset(&header, 0, sizeof(header));

    snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "%s(filename, data, [subregion_in, magic, comment, emdata, userdata])", mexFunctionName());
    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis);
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ ".\n");
        return;
    }

    if (nlhs>0 || nrhs<2 || nrhs>7) {
        snprintf(s, __MAX_S__LENGTH__, "%s: Wrong number of in/output arguments: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }

    {
        #define __MAX_UINT32_T__ 0xFFFFFFFF
        const mxArray *arg;
        mxArray *tmpArray = NULL;
        const mwSize *size;
        double *pdouble;
        int8_t *pint8;
        int32_t *pint32;
        mxClassID type;

        /* filename */
        arg = prhs[0];
        if (!mxIsChar(arg) || mxGetNumberOfDimensions(arg)!=2 || mxGetM(arg)!=1 || (filename_length=mxGetN(arg))<1 || filename_length>=__MAX_FILENAME_LENGTH__) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs file name as first parameter: %s.", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        mxGetString(arg, filename, __MAX_FILENAME_LENGTH__);


        /* subregion */
        if (nrhs >= 3) {
            arg = prhs[2];
            i = mxGetM(arg); j = mxGetN(arg);
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2 || !((i==0&&j==0) || (i==1&&j==6) || (i==6&&j==1))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The subregion must be either [] or a 1x6 vector: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (i!=0) {
                if (!(tmpArray=getDoubleArray(arg))) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the subregion. Maybe out of memory.", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pdouble = mxGetPr(tmpArray);
                for (k=0; k<6; k++) {
                    if (pdouble[k] != floor(pdouble[k]) || pdouble[k]<0 || (k>=3&&pdouble[k]<=0) || pdouble[k]>__MAX_UINT32_T__) {
                        snprintf(s, __MAX_S__LENGTH__, "%s: The subregion must contain positive integer values: %s", mexFunctionName(), synopsis);
                        mexErrMsgTxt(s);
                    }
                    subregion_field[k] = (uint32_t)pdouble[k];
                }
                subregion = subregion_field;
                mxDestroyArray(tmpArray);
            }
        }



        /* magic */
        if (nrhs >= 4) {
            arg = prhs[3];
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2) {
                snprintf(s, __MAX_S__LENGTH__, "%s: magic must be a 4-vector of int8 or [].", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (mxGetNumberOfElements(arg) > 0) {
                if (mxGetClassID(arg)!=mxINT8_CLASS || mxGetNumberOfElements(arg)!=4 || (mxGetN(arg)!=1&&mxGetM(arg)!=1)) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: magic must be a 4-vector of int8 or [].", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pint8 = (int8_t *)mxGetData(arg);
                header.machine = pint8[0];
                header.byte2 = pint8[1];
                header.byte3 = pint8[2];
                header.type = pint8[3];
                autoval_magic = 0;
            }
        }

        /* comment */
        if (nrhs >= 5) {
            arg = prhs[4];
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
                for (i=0; i<80; i++) { header.comment[i] = pint8[i]; }
            }
        }

        /* emdata */
        if (nrhs >= 6) {
            arg = prhs[5];
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
                for (i=0; i<40; i++) { header.emdata[i] = pint32[i]; }
            }
        }


        /* userdata */
        if (nrhs >= 7) {
            arg = prhs[6];
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
                for (i=0; i<256; i++) { header.userdata[i] = pint8[i]; }
            }
        }

        arg = prhs[1];
        if (!mxIsNumeric(arg) || mxGetNumberOfDimensions(arg)>3 || mxGetNumberOfElements(arg)<1) {
            snprintf(s, __MAX_S__LENGTH__, "%s: The volume must be a non-empty numeric array.", mexFunctionName());
            mexErrMsgTxt(s);
        }
        size = mxGetDimensions(arg);
        dims[0] = size[0];
        dims[1] = size[1];
        dims[2] = mxGetNumberOfDimensions(arg)==3 ? size[2] : 1;

        type = mxGetClassID(arg);
        if (mxIsComplex(arg)) {
            const size_t numel = dims[0]*dims[1]*dims[2];
            if (!autoval_magic && tom_io_em_get_iotype_header(&header)!=TOM_IO_TYPE_COMPLEX32 && tom_io_em_get_iotype_header(&header)!=TOM_IO_TYPE_COMPLEX64) {
                snprintf(s, __MAX_S__LENGTH__, "%s: Saving complex data is only possible as complex (single) floating point.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (!(complexCopy = malloc(numel*2*(type==mxDOUBLE_CLASS?sizeof(double):sizeof(float))))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: Error allocating memory for temporary copy of complex data.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (type == mxSINGLE_CLASS) {
                float *pdst = (float *)complexCopy;
                const float *psrc_re, *psrc_im;
                w_iotype = TOM_IO_TYPE_COMPLEX32;
                psrc_re = (float *)mxGetData(arg);
                psrc_im = (float *)mxGetImagData(arg);
                for (i=0; i<numel; i++) {
                    *pdst++ = psrc_re[i];
                    *pdst++ = psrc_im[i];
                }
            } else if (type == mxDOUBLE_CLASS) {
                double *pdst = (double *)complexCopy;
                const double *psrc_re, *psrc_im;
                w_iotype = TOM_IO_TYPE_COMPLEX64;
                psrc_re = mxGetPr(arg);
                psrc_im = mxGetPi(arg);
                for (i=0; i<numel; i++) {
                    *pdst++ = psrc_re[i];
                    *pdst++ = psrc_im[i];
                }
            } else {
                free(complexCopy);
                snprintf(s, __MAX_S__LENGTH__, "%s: EM supports only floating point complex data (single or double).", mexFunctionName(), type);
                mexErrMsgTxt(s);
            }
            w_data = complexCopy;
        } else {
            switch (type) {
                case mxINT8_CLASS:
                    w_iotype = TOM_IO_TYPE_INT8;
                    break;
                case mxINT16_CLASS:
                    w_iotype = TOM_IO_TYPE_INT16;
                    break;
                case mxINT32_CLASS:
                    w_iotype = TOM_IO_TYPE_INT32;
                    break;
                case mxSINGLE_CLASS:
                    w_iotype = TOM_IO_TYPE_FLOAT;
                    break;
                case mxDOUBLE_CLASS:
                    w_iotype = TOM_IO_TYPE_DOUBLE;
                    break;
                default:
                    snprintf(s, __MAX_S__LENGTH__, "%s: The volume has type %s: This can not be saved to EM without conversion.", mexFunctionName(), type);
                    mexErrMsgTxt(s);
            }
            w_data = mxGetData(arg);
        }
        if (subregion) {
            if (subregion[0]+subregion[3]>dims[0] || subregion[1]+subregion[4]>dims[1] || subregion[2]+subregion[5]>dims[2]) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The subregion is out of the volume.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            header.dims[0] = subregion[3];
            header.dims[1] = subregion[4];
            header.dims[2] = subregion[5];
            w_data = (const char *)w_data + ((((subregion[2]*dims[1] + subregion[1]))*dims[0] + subregion[0])*tom_io_iotype_datasize(w_iotype));
            stridey = dims[0] * tom_io_iotype_datasize(w_iotype);
            stridez = dims[1] * stridey;
        } else {
            if (dims[0]>=__MAX_UINT32_T__ || dims[1]>=__MAX_UINT32_T__ || dims[2]>=__MAX_UINT32_T__) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The volume is too large to save it to EM.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            header.dims[0] = dims[0];
            header.dims[1] = dims[1];
            header.dims[2] = dims[2];
        }
        if (autoval_magic) {
            header.machine = 6;
            tom_io_em_set_iotype_header(&header, w_iotype);
        }
        if (!tom_io_em_is_valid_header(&header, sizeof(header))) {
            if (complexCopy) { free(complexCopy); complexCopy=0; }
            snprintf(s, __MAX_S__LENGTH__, "%s: The given EM-Header is not valid. Check the content of magic.", mexFunctionName());
            mexErrMsgTxt(s);
        }
        #undef __MAX_UINT32_T__
    }

    {
        /*printf("WRITE: swapped: %d,    %d->%d %20.15f\n", tom_io_em_is_swaped(&header), w_iotype, tom_io_em_get_iotype_header(&header), *((double *)w_data));*/
        int i;
        if ((i=tom_io_em_write(filename, &header, w_data, w_iotype, 0, stridey, stridez)) != TOM_ERR_OK) {
            const int errno_ = errno;
            char s2[__MAX_S__LENGTH__];
            if (complexCopy) { free(complexCopy); complexCopy=0; }

            tom_io_strerror(i, errno_, s2, __MAX_S__LENGTH__);
            s2[__MAX_S__LENGTH__-1] = 0;
            snprintf(s, __MAX_S__LENGTH__, "%s: %s (writing '%s').", mexFunctionName(), s2, filename);
            mexErrMsgTxt(s);
        }
    }

    if (complexCopy) { free(complexCopy); complexCopy=0; }
}



