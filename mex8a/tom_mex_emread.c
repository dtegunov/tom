/***********************************************************************//**
 * \file tom_mex_emread.c
 * \brief MEX-Function for em_read, replacing tom_emread.
 * \author  Thomas Haller
 * \version 0.2
 * \date    25.10.2007
 **************************************************************************/





#include <mex.h>
#include <io64.h>

#include <math.h>
#include <string.h>
#include <errno.h>

#ifdef _WIN32
    #define __PATH_SEP '\\'
    #define snprintf _snprintf
    #include <direct.h>
#else
    #define __PATH_SEP '/'
    #include <unistd.h>
#endif


#include <tom/io/io.h>
#include <tom/tools/mex_helpfcn.h>







/***********************************************************************//**
 * \brief Splits the filename in directory name and filename.
 *
 * \param[in,out] buffer
 * \param[in] size Size of the buffer (it must be at least the length of
 *   the path+1
 **************************************************************************/
static int getCwd(char *buffer, size_t size) {
    if (!buffer || size<=1) {
        return 0;
    }
#if (defined _UNISTD_H_) || (defined _UNISTD_H)
   return getcwd(buffer, size-1) != NULL;
#elif (defined _DIRECT_H_) || (defined _INC_DIRECT)
   return _getcwd(buffer, (int)size-1) != NULL;
#else
    #error HALLO
#endif
}


static char *getAbsoluteName(const char *filename) {
    size_t lfilename;
    char *res;
    if (!filename) {
        return NULL;
    }
    lfilename = strlen(filename);
#ifdef _WIN32
    if (lfilename>2 && filename[1] == ':') {
#else
    if (lfilename>1 && filename[0]== __PATH_SEP) {
#endif
        res = (char *)malloc(lfilename+1);
        strcpy(res, filename);
    } else {
        #define __getAbs_maxlength 4096
        char buffer[__getAbs_maxlength];
        if (!getCwd(buffer, __getAbs_maxlength)) {
            res = (char *)malloc(lfilename+1);
            strcpy(res, filename);
        } else {
            size_t bufferlength = strlen(buffer);
            int moved;
            char *sub;
            res = (char *)malloc(bufferlength+lfilename+2);
            strcpy(res, buffer);
            res[bufferlength] = __PATH_SEP;
            strcpy(res+bufferlength+1, filename);

            do {
                moved = 0;

                sprintf(buffer, "%c%c", __PATH_SEP, __PATH_SEP);
                sub = strstr(res, buffer);
                if (sub) {
                    memmove(sub, sub+strlen(buffer)-1, strlen(sub+strlen(buffer)-1)+1);
                    moved = 1;
                }

                sprintf(buffer, "%c%c%c", __PATH_SEP, '.', __PATH_SEP);
                sub = strstr(res, buffer);
                if (sub) {
                    memmove(sub, sub+strlen(buffer)-1, strlen(sub+strlen(buffer)-1)+1);
                    moved = 1;
                }
            } while (moved);
        }
    }
    return res;
}

/***********************************************************************//**
 * \brief Splits the filename in directory name and filename.
 *
 * \param[in] filename
 * \param[out] dir_name
 * \param[out] file_name
 **************************************************************************/
static void split_filename(const char *infilename, char *dir_name, size_t dir_name_size, char *file_name, size_t file_name_size) {
    char *filename;
    int clear_path;
    char *p_char;
    if (!infilename) {
        if (dir_name && dir_name_size) { dir_name[0]  = 0; }
        if (file_name && file_name_size) { file_name[0] = 0; }
        return;
    }
    filename = getAbsoluteName(infilename);
    if (!filename) {
        return;
    }
    if (__PATH_SEP == '\\') {
        char *pfilename = filename;
        while ((pfilename = strchr(pfilename, '/'))) {
            *pfilename++ = __PATH_SEP;
        }
    }

    clear_path = 0;
    {
        p_char = strrchr(filename, __PATH_SEP);
        if (!p_char) {
            clear_path = 1;
            if (strlen(filename)+1 <= file_name_size) {
                strcpy(file_name, filename);
            } else {
                if (file_name_size) {
                    file_name[0] = 0;
                }
            }
        } else {
            if ((size_t)(p_char-filename)+1 < dir_name_size) {
                strncpy(dir_name, filename, (size_t)(p_char-filename)+1);
                dir_name[(size_t)(p_char-filename)+1] = 0;
            } else {
                clear_path = 1;
            }
            if (strlen(p_char+1) && strlen(p_char+1)<file_name_size)  {
                strcpy(file_name, p_char+1);
            } else {
                if (file_name_size) {
                    file_name[0] = 0;
                }
            }
        }
    }
    if (dir_name && clear_path) {
        switch (dir_name_size) {
            case 0:
                break;
            case 1:
                dir_name[0] = 0;
                break;
            case 2:
                dir_name[0] = '.';
                dir_name[1] = 0;
                break;
            default:
                dir_name[0] = '.';
                dir_name[1] = __PATH_SEP;
                dir_name[2] = 0;
        }
    }
    free(filename);
}




/***********************************************************************//**
 * \brief mexFunction to read an em-file.
 *
 * Similar to tom_emread.m this reads an emfile from disk.
 * This mex-interface expects 4 input arguments passed from matlab.
 * The first parameter is the filename, the second one is the subregion,
 * the 3rd is the sampling factor and the last parameter is the binning
 * factor.
 **************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {


    size_t filename_length = 0;
    #define __MAX_FILENAME_LENGTH__ ((int)2048)
    char filename[__MAX_FILENAME_LENGTH__+1];
    #define __MAX_S__LENGTH__ (__MAX_FILENAME_LENGTH__+1024)
    char s[__MAX_S__LENGTH__+1];

    uint32_t subregion_field[6];
    uint32_t binning_field[3];
    uint32_t sampling_field[3];
    const uint32_t *subregion = NULL;
    const uint32_t *binning = NULL;
    const uint32_t *sampling = NULL;
    tom_io_em_header header;
    FILE *f;
    uint32_t dims_l[3];
    int iotype, restype;
    #define __BUFFERLENGTH_SYNOPSIS__ 1024
    char synopsis[__BUFFERLENGTH_SYNOPSIS__+1];


    mxArray *mxData;

    size_t i, j, k;
    int res;


    snprintf(synopsis, __BUFFERLENGTH_SYNOPSIS__, "[vol, size, magic, comment, emdata, userdata] = %s(filename, [subregion, sampling, binning])", mexFunctionName());

    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis);
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ ".\n");
        return;
    }




    if (nrhs < 1 || nrhs > 4) {
        snprintf(s, __MAX_S__LENGTH__, "%s: call function with up to 4 parameters: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }
    if (nlhs>6) {
        snprintf(s, __MAX_S__LENGTH__, "%s: Too many output parameters: %s.", mexFunctionName(), synopsis);
        mexErrMsgTxt(s);
    }


    {
        #define __MAX_UINT32_T__ 0xFFFFFFFF
        const mxArray *PRHS_FILENAME = prhs[0];
        const mxArray *PRHS_SUBREGION = prhs[1];
        const mxArray *PRHS_SAMPLING = prhs[2];
        const mxArray *PRHS_BINNING = prhs[3];

        /* Check input parameters! */

        mxArray *tmpArray = NULL;
        const double *pdouble;


        if (!mxIsChar(PRHS_FILENAME) || mxGetNumberOfDimensions(PRHS_FILENAME)!=2 || mxGetM(PRHS_FILENAME)!=1 || (filename_length=mxGetN(PRHS_FILENAME))<1) {
            snprintf(s, __MAX_S__LENGTH__, "%s: needs the file name as first parameter: %s.", mexFunctionName(), synopsis);
            mexErrMsgTxt(s);
        }
        if (filename_length >= __MAX_FILENAME_LENGTH__) {
            snprintf(s, __MAX_S__LENGTH__, "%s: Maximal length of file name exceeded. (Recompile with larger buffer :)", mexFunctionName());
            mexErrMsgTxt(s);
        }
        mxGetString(PRHS_FILENAME, filename, __MAX_FILENAME_LENGTH__);

        if (nrhs >= 2) {
            i = mxGetM(PRHS_SUBREGION); j = mxGetN(PRHS_SUBREGION);
            if (!mxIsNumeric(PRHS_SUBREGION) || mxIsComplex(PRHS_SUBREGION) || mxGetNumberOfDimensions(PRHS_SUBREGION)!=2 ||
                !((i==0&&j==0) || (i==1&&j==6) || (i==6&&j==1))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The subregion must be either [] or a 1x6 vector: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (i!=0) {
                if (!(tmpArray=getDoubleArray(PRHS_SUBREGION))) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the subregion. Maybe out of memory.", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pdouble = mxGetPr(tmpArray);
                for (k=0; k<6; k++) {
                    if (pdouble[k] != floor(pdouble[k]) || pdouble[k] < 0 || pdouble[k]>__MAX_UINT32_T__) {
                        snprintf(s, __MAX_S__LENGTH__, "%s: The subregion must contain positive integer values: %s", mexFunctionName(), synopsis);
                        mexErrMsgTxt(s);
                    }
                    subregion_field[k] = (uint32_t)pdouble[k];
                }
                subregion = subregion_field;
                mxDestroyArray(tmpArray);
            }
        }


        if (nrhs >= 3) {
            i = mxGetM(PRHS_SAMPLING); j = mxGetN(PRHS_SAMPLING);
            if (!mxIsNumeric(PRHS_SAMPLING) || mxIsComplex(PRHS_SAMPLING) || mxGetNumberOfDimensions(PRHS_SAMPLING)!=2 ||
                !((i==0&&j==0) || (i==1&&j==3) || (i==3&&j==1))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The sampling factor must be either [] or a 1x3 vector: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (i!=0) {
                if (!(tmpArray=getDoubleArray(PRHS_SAMPLING))) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the sampling factor. Maybe out of memory.", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pdouble = mxGetPr(tmpArray);
                for (k=0; k<3; k++) {
                    if (pdouble[k] != floor(pdouble[k]) || pdouble[k] < 0 || pdouble[k]>__MAX_UINT32_T__) {
                        snprintf(s, __MAX_S__LENGTH__, "%s: The sampling must contain positive integer values: %s", mexFunctionName(), synopsis);
                        mexErrMsgTxt(s);
                    }
                    sampling_field[k] = (uint32_t)pdouble[k];
                }
                sampling = sampling_field;
                mxDestroyArray(tmpArray);
            }
        }

        if (nrhs >= 4) {
            i = mxGetM(PRHS_BINNING); j = mxGetN(PRHS_BINNING);
            if (!mxIsNumeric(PRHS_BINNING) || mxIsComplex(PRHS_BINNING) || mxGetNumberOfDimensions(PRHS_BINNING)!=2 ||
                !((i==0&&j==0) || (i==1&&j==3) || (i==3&&j==1))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: The binning factor must be either [] or a 1x3 vector: %s", mexFunctionName(), synopsis);
                mexErrMsgTxt(s);
            }
            if (i!=0) {
                if (!(tmpArray=getDoubleArray(PRHS_BINNING))) {
                    snprintf(s, __MAX_S__LENGTH__, "%s: Error creating a copy of the binning factor. Maybe out of memory.", mexFunctionName());
                    mexErrMsgTxt(s);
                }
                pdouble = mxGetPr(tmpArray);
                for (k=0; k<3; k++) {
                    if (pdouble[k] != floor(pdouble[k]) || pdouble[k] < 0 || pdouble[k]>__MAX_UINT32_T__) {
                        snprintf(s, __MAX_S__LENGTH__, "%s: The sampling must contain positive integer values: %s", mexFunctionName(), synopsis);
                        mexErrMsgTxt(s);
                    }
                    binning_field[k] = (uint32_t)pdouble[k];
                }
                binning = binning_field;
                mxDestroyArray(tmpArray);
            }
        }
        #undef __MAX_UINT32_T__
    }




    if ((res = tom_io_em_read_header(filename, "rb", &header, &f)) != TOM_ERR_OK) {
        char s2[__MAX_S__LENGTH__];
        tom_io_strerror(res, errno, s2, __MAX_S__LENGTH__);
        snprintf(s, __MAX_S__LENGTH__, "%s: Reading file \"%s\": %s", mexFunctionName(), filename, s2);
        mexErrMsgTxt(s);
    }

    if ((res=tom_io_calculate_sizes(header.dims, subregion, subregion_field, sampling, sampling_field, binning, binning_field, NULL, dims_l)) != TOM_ERR_OK) {
        fclose(f);
        if (res == TOM_ERR_SUBREGION_OUT) {
            snprintf(s, __MAX_S__LENGTH__, "%s: The subregion is out of the %ux%ux%u-volume.", mexFunctionName(), header.dims[0], header.dims[1], header.dims[2]);
        } else if (res == TOM_ERR_BINNING_TOO_HIGH) {
            snprintf(s, __MAX_S__LENGTH__, "%s: The binning factor is to high for the %ux%ux%u-volume.", mexFunctionName(), header.dims[0], header.dims[1], header.dims[2]);
        } else {
            char s2[__MAX_S__LENGTH__];
            tom_io_strerror(res, errno, s2, __MAX_S__LENGTH__);
            snprintf(s, __MAX_S__LENGTH__, "%s: Reading file %s: An unexpected error is returned when calculating the resulting size of the subregion: %s", mexFunctionName(), filename, s2);
        }
        mexErrMsgTxt(s);
    }


    iotype = tom_io_em_get_iotype_header(&header);
    restype = iotype;

    switch (iotype) {
        case TOM_IO_TYPE_INT8:
        case TOM_IO_TYPE_INT16:
        case TOM_IO_TYPE_INT32:
            restype = (binning_field[0]>1 || binning_field[1]>1 || binning_field[2]>1) ? TOM_IO_TYPE_DOUBLE : iotype;
            break;
        case TOM_IO_TYPE_FLOAT:
            restype = (binning_field[0]>1 || binning_field[1]>1 || binning_field[2]>1) ? TOM_IO_TYPE_DOUBLE : iotype;
            break;
        case TOM_IO_TYPE_COMPLEX32:
        case TOM_IO_TYPE_DOUBLE:
        case TOM_IO_TYPE_COMPLEX64:
            /* binning is not supported for complex type (up to now...) */
            restype = iotype;
            break;
        default:
            fclose(f);
            snprintf(s, __MAX_S__LENGTH__, "%s: data-type of the header (%d) is not recognized. This is probably an implementation error somewhere.", mexFunctionName(), iotype);
            mexErrMsgTxt(s);
            break;
    }

    {
        const mxClassID classID = getMxClassIDFromEM(restype);
        mwSize mx_dims[3];
        mx_dims[0] = dims_l[0];
        mx_dims[1] = dims_l[1];
        mx_dims[2] = dims_l[2];
        if (restype==TOM_IO_TYPE_COMPLEX32 || restype==TOM_IO_TYPE_COMPLEX64) {
            mx_dims[2] *= 2;
            mxData = mxCreateNumericArray(3, mx_dims, classID, mxREAL);
        } else {
            mxData = mxCreateNumericArray(mx_dims[2]==1 ? 2 : 3, mx_dims, classID, mxREAL);
        }
        if (!mxData) {
            fclose(f);
            snprintf(s, __MAX_S__LENGTH__, "%s: Error allocating memory for the result.", mexFunctionName());
            mexErrMsgTxt(s);
        }

        if ((res=tom_io_read_vol(f, header.dims, iotype, tom_io_em_is_swaped(&header), subregion_field, sampling_field, binning_field, mxGetData(mxData), restype)) != TOM_ERR_OK) {
            const int errno_ = errno;
            char error_text[__MAX_S__LENGTH__];
            fclose(f);
            tom_io_strerror(res, errno_, error_text, __MAX_S__LENGTH__);
            error_text[__MAX_S__LENGTH__-1] = 0;
            snprintf(s, __MAX_S__LENGTH__, "%s: %s (reading %s).", mexFunctionName(), error_text, filename);
            mexErrMsgTxt(s);
        }

        fclose(f);

        if (restype==TOM_IO_TYPE_COMPLEX32 || restype==TOM_IO_TYPE_COMPLEX64) {
            /* The em-format saves the complex data interleaved.
               Split it up, into two new arrays as needed in matlab. */
            mxArray *mxData_as_real = mxData;
            size_t numel;

            mx_dims[2] /= 2;
            if (!(mxData = mxCreateNumericArray(mx_dims[2]==1 ? 2 : 3, mx_dims, classID, mxCOMPLEX))) {
                snprintf(s, __MAX_S__LENGTH__, "%s: Error allocating memory for the complex result.", mexFunctionName());
                mexErrMsgTxt(s);
            }
            numel = mx_dims[0] * mx_dims[1] * mx_dims[2];
            if (restype==TOM_IO_TYPE_COMPLEX32) {
                const float *data_as_real;
                float *data_r, *data_i;
                data_as_real = (const float *)mxGetData(mxData_as_real);
                data_r = (float *)mxGetData(mxData);
                data_i = (float *)mxGetImagData(mxData);
                for (i=0; i<numel; i++) {
                    *data_r++ = *data_as_real++;
                    *data_i++ = *data_as_real++;
                }
            } else {
                const double *data_as_real;
                double *data_r, *data_i;
                data_as_real = (const double *)mxGetData(mxData_as_real);
                data_r = (double *)mxGetData(mxData);
                data_i = (double *)mxGetImagData(mxData);
                for (i=0; i<numel; i++) {
                    *data_r++ = *data_as_real++;
                    *data_i++ = *data_as_real++;
                }
            }
            mxDestroyArray(mxData_as_real);
        }
    }

    {
        mxArray *p;
        void *pdata;
        switch (nlhs) {
            case 6:
                if (!(p = mxCreateNumericMatrix(1, 256, mxUINT8_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
                pdata = mxGetData(plhs[5] = p);
                for (i=0; i<256; i++) { ((uint8_t *)pdata)[i] = header.userdata[i]; }
            case 5:
                if (!(p = mxCreateNumericMatrix(1, 40, mxINT32_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
                pdata = mxGetData(plhs[4] = p);
                for (i=0; i<40; i++) { ((int32_t *)pdata)[i] = header.emdata[i]; }
            case 4:
                if (!(p = mxCreateNumericMatrix(1, 80, mxINT8_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
                pdata = mxGetData(plhs[3] = p);
                for (i=0; i<80; i++) { ((int8_t *)pdata)[i] = header.comment[i]; }
            case 3:
                if (!(p = mxCreateNumericMatrix(1, 4, mxUINT8_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
                pdata = mxGetData(plhs[2] = p);
                ((uint8_t *)pdata)[0] = header.machine;
                ((uint8_t *)pdata)[1] = header.byte2;
                ((uint8_t *)pdata)[2] = header.byte3;
                ((uint8_t *)pdata)[3] = header.type;
            case 2:
                if (!(p = mxCreateNumericMatrix(3, 1, mxUINT32_CLASS, mxREAL))) { mexErrMsgTxt("Error allocating memory"); }
                pdata = mxGetData(plhs[1] = p);
                ((uint32_t *)pdata)[0] = header.dims[0];
                ((uint32_t *)pdata)[1] = header.dims[1];
                ((uint32_t *)pdata)[2] = header.dims[2];
            case 1:
            case 0:
                plhs[0] = mxData;
                break;
            default:
                break; /* not expected. */
        }
    }
}



