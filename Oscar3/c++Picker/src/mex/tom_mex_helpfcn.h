/****************************************************************************//**
 * \file tom_mex_helpfcn.h
 * \author Thomas Haller
 * \date Oct. 26 2007
 * \brief Contains some helper functions used in mex-files.
 *
 * Only to be included in mex-files, thus it does not hurt that there are
 * function definitions.
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__MEX_HELPFCN_H__
#define ___INCLUDE_CORE__MEX_HELPFCN_H__




#ifdef __GNUC__
    #include <stdint.h>
#else
    #include "system/pstdint.h"
#endif


#include "tom/core/defines.h"





/****************************************************************************//**
 * \brief Copies the content of an numerical mxArray into a Double-mxArray
 *
 * \param[in] src A pointer to the soure array.
 * \return A copied version of the original array of type mxDOUBLE_CLASS.
 *   It must be freed later calling mxDestroyArray. In case of an error,
 *   NULL is returned.
 *******************************************************************************/
mxArray *getDoubleArray(const mxArray *src) {

    mxClassID srcClassID;
    mxArray *dst;

    if (!src || mxIsSparse(src) || !mxIsNumeric(src)) {
        return NULL;
    }

    srcClassID = mxGetClassID(src);

    if (srcClassID == mxDOUBLE_CLASS) {
        dst = mxDuplicateArray(src);
    } else {
        double *pdst;
        const void *psrc;
        int cpart;
        size_t i;

        const size_t numel = mxGetM(src) * mxGetN(src);
        dst = mxCreateNumericArray(mxGetNumberOfDimensions(src), mxGetDimensions(src), mxDOUBLE_CLASS, mxIsComplex(src) ? mxCOMPLEX : mxREAL);


        for (cpart=0; cpart<1; cpart++) {
            if (cpart == 0) {
                psrc = mxGetData(src);
                pdst = mxGetPr(dst);
            } else {
                if (!mxIsComplex(src)) {
                    continue;
                }
                psrc = mxGetImagData(src);
                pdst = mxGetPi(dst);
            }
            switch (srcClassID) {
                case mxSINGLE_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const float    *)psrc)[i]; }
                    break;
                case mxINT8_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const int8_t   *)psrc)[i]; }
                    break;
                case mxUINT8_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const uint8_t  *)psrc)[i]; }
                    break;
                case mxINT16_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const int16_t  *)psrc)[i]; }
                    break;
                case mxUINT16_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const uint16_t *)psrc)[i]; }
                    break;
                case mxINT32_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const int32_t  *)psrc)[i]; }
                    break;
                case mxUINT32_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = ((const uint32_t *)psrc)[i]; }
                    break;
                case mxINT64_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = (double)((const int64_t  *)psrc)[i]; }
                    break;
                case mxUINT64_CLASS:
                    for (i=0; i<numel; i++) { pdst[i] = (double)((const uint64_t *)psrc)[i]; }
                    break;
                default:
                    mxDestroyArray(dst);
                    return NULL;
            }
        }
    }

    return dst;
}





/***********************************************************************//**
 * \brief Convert numerical constant to mxClassID
 *
 * Gets the numerical constants representing the a datatype as used in
 * tom_io.h. It is one of the defines from tom_io.h and is named something
 * like __TOM_IOTYPE_XXX
 * \param[in] iotype The numerical representation of the data type.
 * \return The corresponding class id of the MEX-interface.
 * mxUNKNOWN_CLASS if the type is not recognized.
 **************************************************************************/
mxClassID getMxClassIDFromEM(int iotype) {
    switch (iotype) {
        case TOM_IO_TYPE_INT8:
            return mxINT8_CLASS;
        case TOM_IO_TYPE_INT16:
            return mxINT16_CLASS;
        case TOM_IO_TYPE_INT32:
            return mxINT32_CLASS;
        case TOM_IO_TYPE_COMPLEX4:
        case TOM_IO_TYPE_FLOAT:
            return mxSINGLE_CLASS;
        case TOM_IO_TYPE_DOUBLE:
            return mxDOUBLE_CLASS;
    }
    return mxUNKNOWN_CLASS;
}


mxComplexity getMxComplexityFromEM(int iotype) {
    switch (iotype) {
        case TOM_IO_TYPE_COMPLEX4:
            return mxCOMPLEX;
        default:
            return mxREAL;
    }
}

int getIOTypeFromMxClassID(mxClassID id, int isComplex) {
    if (isComplex) {
        switch (id) {
            case mxSINGLE_CLASS:
                return TOM_IO_TYPE_COMPLEX4;
        }
    } else {
        switch (id) {
            case mxINT8_CLASS:
                return TOM_IO_TYPE_INT8;
            case mxINT16_CLASS:
                return TOM_IO_TYPE_INT16;
            case mxINT32_CLASS:
                return TOM_IO_TYPE_INT32;
            case mxSINGLE_CLASS:
                return TOM_IO_TYPE_FLOAT;
            case mxDOUBLE_CLASS:
                return TOM_IO_TYPE_DOUBLE;
        }
    }
    return TOM_IO_TYPE_UNKNOWN;
}


#ifdef __cplusplus
// Declare c++ functions.



#include <sstream>
#include <boost/cast.hpp>
#include <stdexcept>

template<typename T>
inline T getScalar(const mxArray *a) {

    if (mxGetNumberOfElements(a) != 1) {
        throw std::invalid_argument("The mxArray does not contain a single scalar.");
    }
    if (mxIsComplex(a)) {
        throw std::invalid_argument("getScalar does not allow complex data.");
    }

    T res;

    switch (mxGetClassID(a)) {
        case mxINT8_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<int8_t *>(mxGetData(a)));
            break;
        case mxINT16_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<int16_t *>(mxGetData(a)));
            break;
        case mxINT32_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<int32_t *>(mxGetData(a)));
            break;
        case mxINT64_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<int64_t *>(mxGetData(a)));
            break;
        case mxUINT8_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<uint8_t *>(mxGetData(a)));
            break;
        case mxUINT16_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<uint16_t *>(mxGetData(a)));
            break;
        case mxUINT32_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<uint32_t *>(mxGetData(a)));
            break;
        case mxUINT64_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<uint64_t *>(mxGetData(a)));
            break;
        case mxSINGLE_CLASS:
            {
                const float &d = *reinterpret_cast<float *>(mxGetData(a));
                if (std::numeric_limits<T>::is_integer) {
                    if (ceil(d) != d) {
                        throw std::bad_cast();
                    }
                }
                res = boost::numeric_cast<T>(d);
            }
            break;
        case mxDOUBLE_CLASS:
            {
                const double &d = *reinterpret_cast<double *>(mxGetData(a));
                if (std::numeric_limits<T>::is_integer) {
                    if (ceil(d) != d) {
                        throw std::bad_cast();
                    }
                }
                res = boost::numeric_cast<T>(d);
            }
            break;
        case mxLOGICAL_CLASS:
            res = boost::numeric_cast<T>(*reinterpret_cast<mxLogical *>(mxGetData(a)));
            break;
        case mxCHAR_CLASS:
        case mxFUNCTION_CLASS:
        case mxSTRUCT_CLASS:
        case mxCELL_CLASS:
        case mxUNKNOWN_CLASS:
        default:
            throw std::bad_cast();
    }
    return res;
}


#endif

#endif



