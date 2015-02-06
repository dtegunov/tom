/***********************************************************************//**
 * \file tom_mex_emwrite_paste.c
 * \brief MEX-Function for writing an em-file by pasting it.
 * \author  Thomas Haller
 * \version 0.2
 * \date    22.01.2008
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

    #define __MAX_S__LENGTH__ (4096)
    char s[__MAX_S__LENGTH__+1];

    #define __BUFFERLENGTH_SYNOPSIS__ 1024
    char synopsis[__BUFFERLENGTH_SYNOPSIS__];

    size_t i, j, k;

    const void *w_data;
    int w_iotype;

    uint32_t subregion_start[3] = { 0,0,0 };
    uint32_t reverse_sampling[3] = { 0,0,0 };
    tom_io_em_header header;
    int allow_conversion = 1;
    mxArray *plhs_tmp[5] = { NULL, NULL, NULL, NULL, NULL };
    void *complexCopy = NULL;

    uint32_t w_sizex, w_sizey, w_sizez;
    size_t w_stridex=0, w_stridey=0, w_stridez=0;




    snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "[size, magic, comment, emdata, userdata] = %s(filename, volume, [subregion_start, reverse_sampling, allow_conversion])", mexFunctionName());
    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis);
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ ".\n");
        return;
    }

    if (nlhs>5 || nrhs<2 || nrhs>5) {
        snprintf(s, __MAX_S__LENGTH__, "%s: Wrong number of in/output arguments: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }

    {
        #define __MAX_UINT32_T__ 0xFFFFFFFF
        const mxArray *arg;
        mxArray *tmpArray = NULL;
        const mwSize *size;
        double *pdouble;
        mxClassID type;

        /* filename */
        arg = prhs[0];
        if (!mxIsChar(arg) || mxGetNumberOfDimensions(arg)!=2 || mxGetM(arg)!=1 || (filename_length=mxGetN(arg))<1 || filename_length>=__MAX_FILENAME_LENGTH__) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs file name as first parameter: %s.", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        mxGetString(arg, filename, __MAX_FILENAME_LENGTH__);


        /* subregion_start */
        if (nrhs >= 3) {
            arg = prhs[2];
            i = mxGetM(arg); j = mxGetN(arg);
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2 || !((i==0&&j==0) || (i==1&&j==3) || (i==3&&j==1))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The subregion_start must be either [] or a 1x3 vector: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (i!=0) {
                if (!(tmpArray=getDoubleArray(arg))) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the subregion_start. Maybe out of memory.", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pdouble = mxGetPr(tmpArray);
                for (k=0; k<3; k++) {
                    subregion_start[k] = (uint32_t)pdouble[k];
                    if (pdouble[k] != subregion_start[k]) {
                        snprintf(s, __MAX_S__LENGTH__, "%s: subregion_start must contain integer values: %s", mexFunctionName(), synopsis);
                        mexErrMsgTxt(s);
                    }
                }
                mxDestroyArray(tmpArray);
            }
        }

        /* reverse_sampling */
        if (nrhs >= 4) {
            arg = prhs[3];
            i = mxGetM(arg); j = mxGetN(arg);
            if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)!=2 || !((i==0&&j==0) || (i==1&&j==3) || (i==3&&j==1))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The reverse_sampling must be either [] or a 1x3 vector: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (i!=0) {
                if (!(tmpArray=getDoubleArray(arg))) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the reverse_sampling. Maybe out of memory.", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pdouble = mxGetPr(tmpArray);
                for (k=0; k<3; k++) {
                    reverse_sampling[k] = (uint32_t)pdouble[k];
                    if (pdouble[k] != reverse_sampling[k]) {
                        snprintf(s, __MAX_S__LENGTH__, "%s: reverse_sampling must contain integer values: %s", mexFunctionName(), synopsis);
                        mexErrMsgTxt(s);
                    }
                }
                mxDestroyArray(tmpArray);
            }
        }
        if (reverse_sampling[0] < 1) { reverse_sampling[0] = 1; }
        if (reverse_sampling[1] < 1) { reverse_sampling[1] = 1; }
        if (reverse_sampling[2] < 1) { reverse_sampling[2] = 1; }


        /* allow_conversion */
        if (nrhs >= 5) {
            const mwSize numel = mxGetNumberOfElements(arg = prhs[4]);
            if (!(mxIsNumeric(arg) || mxIsLogical(arg)) || numel>1) {
                snprintf(s, __MAX_S__LENGTH__, "%s: allow_conversion must be one of true, false or [].", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (numel == 1) {
                allow_conversion = mxGetScalar(arg) != 0;
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


        /* data */
        arg = prhs[1];
        if (!mxIsNumeric(arg) || mxGetNumberOfDimensions(arg)>3) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs a numerical volume as second parameter: %s", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        type = mxGetClassID(arg);
        size = mxGetDimensions(arg);
        w_sizex = size[0];
        w_sizey = size[1];
        w_sizez = mxGetNumberOfDimensions(arg)==3 ? size[2] : 1;
        if (w_sizex!=size[0] || w_sizey!=size[1] || w_sizez!=(mxGetNumberOfDimensions(arg)==3?size[2]:1)) {
            snprintf(s, __MAX_S__LENGTH__, "%s: The volume is too large to save it as EM.", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }

        if (mxIsComplex(arg)) {
            const size_t numel = w_sizex * w_sizey * w_sizez;
            if (!allow_conversion && type!=mxSINGLE_CLASS && type!=mxDOUBLE_CLASS) {
                snprintf(s, __MAX_S__LENGTH__, "%s: Can not convert complex volume (not allow_conversion).", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (type!=mxSINGLE_CLASS && type!=mxDOUBLE_CLASS) {
                snprintf(s, __MAX_S__LENGTH__, "%s: complex volume must be floating point (single or double).", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (!(complexCopy = malloc(numel*2*(type==mxSINGLE_CLASS?sizeof(float):sizeof(double))))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: Error allocating memory for temporary copy of complex data.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            if (type == mxSINGLE_CLASS) {
                float *pdst = complexCopy;
                const float *psrc_re = (float *)mxGetData(arg);
                const float *psrc_im = (float *)mxGetData(arg);
                for (i=0; i<numel; i++) {
                    *pdst++ = psrc_re[i];
                    *pdst++ = psrc_im[i];
                }
                w_iotype = TOM_IO_TYPE_COMPLEX32;
            } else if (type == mxDOUBLE_CLASS) {
                double *pdst = complexCopy;
                const double *psrc_re = mxGetPr(arg);
                const double *psrc_im = mxGetPi(arg);
                for (i=0; i<numel; i++) {
                    *pdst++ = psrc_re[i];
                    *pdst++ = psrc_im[i];
                }
                w_iotype = TOM_IO_TYPE_COMPLEX64;
            } else {
                free(complexCopy); complexCopy = NULL;
                snprintf(s, __MAX_S__LENGTH__, "%s: EM supports only floating point complex data (single or double).", mexFunctionName(), type);
                mexErrMsgTxt(s);
            }
            w_data = complexCopy;
        } else {
            w_data = mxGetData(arg);
            w_iotype = getIOTypeFromMxClassID(type, false);
        }
    }

    {
        int header_read;
        int fcn_res = tom_io_em_write_paste(filename, w_data, w_iotype, w_sizex, w_sizey, w_sizez, w_stridex, w_stridey, w_stridez, &header, &header_read, allow_conversion, subregion_start, reverse_sampling);
        if (fcn_res != TOM_ERR_OK) {
            if (complexCopy) { free(complexCopy); complexCopy = NULL; }
            switch (fcn_res) {
                case TOM_ERR_VOLUME_TOO_LARGE:
                    snprintf(s, __MAX_S__LENGTH__, "%s: The %lux%lux%lu-volume is to large to paste it into the %lux%lux%lu-EM-file (start at %lux%lux%lu, sample %lux%lux%lu).", mexFunctionName(), (unsigned long)w_sizex, (unsigned long)w_sizey, (unsigned long)w_sizez, (unsigned long)header.dims[0], (unsigned long)header.dims[1], (unsigned long)header.dims[2], (unsigned long)subregion_start[0], (unsigned long)subregion_start[1], (unsigned long)subregion_start[2], (unsigned long)reverse_sampling[0], (unsigned long)reverse_sampling[1], (unsigned long)reverse_sampling[2]);
                    break;
                case TOM_ERR_WRONG_IOTYPE_CONVERSION:
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error converting the volume from type %d to type %d (allow conversion = %s).", mexFunctionName(), (int)w_iotype, (int)tom_io_em_get_iotype_header(&header), allow_conversion ? "yes" : "no");
                    break;
                default:
                    snprintf(s, __MAX_S__LENGTH__, "%s: Unspecified error pasting the volume to file (%d).", mexFunctionName(), fcn_res);
                    break;
            }
            mexErrMsgTxt(s);
        }
        if (complexCopy) { free(complexCopy); complexCopy = NULL;}
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



