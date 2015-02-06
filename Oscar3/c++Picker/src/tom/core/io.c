/***********************************************************************//**
 * \file io.c
 * \brief IO-Functions. Contains functions to save and load volumes from file.
 * \author  Thomas Haller
 * \version 0.1
 * \date    23.10.2007
 **************************************************************************/
#include "tom/core/io.h"



#include <string.h>
#include <math.h>
#include <errno.h>
#include <assert.h>


#if defined _WIN32 && defined __LCC__
    #define int64_t int64_T
    #define uint64_t uint64_T
#endif



#if defined _WIN32 || defined _WIN64
    #define ___fseek_cur(f, setpos, success)                                    \
        if (1) {                                                                \
            success = 0==_fseeki64((f), (setpos), SEEK_CUR);                    \
        } else (void)0
    #define ___fseek_end(f, setpos, success)                                    \
        if (1) {                                                                \
            success = 0==_fseeki64((f), (setpos), SEEK_END);                    \
        } else (void)0
    #define ___fseek_set(f, setpos, success)                                    \
        if (1) {                                                                \
            success = 0==_fseeki64((f), (setpos), SEEK_SET);                    \
        } else (void)0
    #define ___ftell(f, getpos)                                                 \
        if (1) {                                                                \
            (getpos) = _ftelli64((f));                                          \
        } else (void)0
#elif !defined __STRICT_ANSI__ || defined __cplusplus
    #define ___fseek_cur(f, setpos, success)                                    \
        if (1) {                                                                \
            success = 0==fseeko((f), (setpos), SEEK_CUR);                      \
        } else (void)0
    #define ___fseek_end(f, setpos, success)                                    \
        if (1) {                                                                \
            success = 0==fseeko((f), (setpos), SEEK_END);                      \
        } else (void)0
    #define ___fseek_set(f, setpos, success)                                    \
        if (1) {                                                                \
            success = 0==fseeko((f), (setpos), SEEK_SET);                      \
        } else (void)0
    #define ___ftell(f, getpos)                                                 \
        if (1) {                                                                \
            (getpos) = ftello((f));                                             \
        } else (void)0
#else
    fpos_t ___fpos;
    #define ___fseek_cur(f, setpos, success)                                    \
        if (1) {                                                                \
            fgetpos(f, &___fpos);                                               \
            ___fpos.__pos += (setpos);                                          \
            success = 0==fsetpos((f), &___fpos);                                \
        } else (void)0
    #define ___fseek_end(f, setpos, success)                                    \
        if (1) {                                                                \
            fseek(f, 0, SEEK_END);                                              \
            fgetpos(f, &___fpos);                                               \
            ___fpos.__pos += (setpos);                                          \
            success = 0==fsetpos((f), &___fpos);                                \
        } else (void)0
    #define ___fseek_set(f, setpos, success)                                    \
        if (1) {                                                                \
            ___fpos.__pos = (setpos);                                           \
            success = 0==fsetpos((f), &___fpos);                                \
        } else (void)0
    #define ___ftell(f, getpos)                                                 \
        if (1) {                                                                \
            fgetpos(f, &___fpos);                                               \
            (getpos) = ___fpos.__pos;                                           \
        } else (void)0
#endif


#if defined MATLAB_MEX_FILE
    #define printf mexPrintf
    #include "mex.h"
#endif



/***********************************************************************//**
 * \brief Returns the byte-order of the running computer.
 *
 * Checks the byte order of the machine where the program runs. It is
 * either little or big endian.
 * \return One of the predefined numerical constants
 *   MACHINE_BYTEORDER_LITTLE_ENDIAN and MACHINE_BYTEORDER_BIG_ENDIAN.
 **************************************************************************/
int tom_io_machine_byteorder() {
    const uint16_t word = 0x0001;
    const char *byte = (const char *)&word;
    return(byte[0] ? TOM_IO_MACHINE_BYTEORDER_LITTLE_ENDIAN : TOM_IO_MACHINE_BYTEORDER_BIG_ENDIAN);
}


/***********************************************************************//**
 * \brief Returns the size of a file in bytes.
 *
 * Returns the filesize in bytes of an open file. It does this by seeking to the end of the file, telling its position and seeking back to the original position.
 * \param[in] file is a pointer to the open file handle.
 * \return Gets the size of the files in bytes.
 **************************************************************************/
int64_t tom_io_getFileSize(FILE *file) {

     int64_t currpos, size;
     int success;

     ___ftell(file, currpos);
     ___fseek_end(file, 0, success);
     ___ftell(file, size);
     ___fseek_set(file, currpos, success);
     return size;
 }



/***********************************************************************//**
 * \brief Swap 2 bytes in memory.
 *
 * Gets a pointer and exchanges the two bytes to which is pointed to.
 * \param[in,out]    x   pointer to the first of the two bytes.
 **************************************************************************/
void tom_io_swap02(void *x) {
    uint16_t *v = (uint16_t *)x;
    *v = (((*v >> 8)) | (*v << 8));
}
/***********************************************************************//**
 * \brief Swap 4 bytes in memory.
 *
 * Gets a pointer and exchanges the 4 bytes to which is pointed to.
 * \param[in,out]    x   pointer to the first of the two bytes.
 **************************************************************************/
void tom_io_swap04(void *x) {
    uint32_t *v = (uint32_t *)x;
    *v = (((*v & 0x000000FF) << 24) +
          ((*v & 0x0000FF00) <<  8) +
          ((*v & 0x00FF0000) >>  8) +
          ((*v & 0xFF000000) >> 24));
}
/***********************************************************************//**
 * \brief Swap 8 bytes in memory.
 *
 * Gets a pointer and exchanges the 8 bytes to which is pointed to.
 * \param[in,out]    x   pointer to the first of the two bytes.
 **************************************************************************/
void tom_io_swap08(void *x) {
    register int i = 0;
    register int j = 7;
    uint8_t tmp;
    uint8_t *v = (uint8_t *)x;
    while (i < j) {
        tmp  = v[i];
        v[i] = v[j];
        v[j] = tmp;
        i++;
        j--;
    }
}
/***********************************************************************//**
 * \brief Swap position of bytes in memory.
 *
 * Gets a pointer and exchanges the position of the next n byte.
 * \param[in,out] x   pointer to the first byte.
 * \param[in]     n   number of bytes to be swaped.
 **************************************************************************/
void tom_io_swap(void *x, size_t size) {
    switch (size) {
        case 0:
        case 1:
            break;
        case 2:
            tom_io_swap02(x);
            break;
        case 4:
            tom_io_swap04(x);
            break;
        default:
            {
                register size_t i = 0;
                register size_t j = size-1;
                uint8_t tmp;
                uint8_t *v = (uint8_t *)x;
                while (i < j) {
                    tmp  = v[i];
                    v[i] = v[j];
                    v[j] = tmp;
                    i++;
                    j--;
                }
            }
    }
}
/***********************************************************************//**
 * \brief Swaps the bytes of 2-bytes elements of an array
 *
 * Gets a pointer to the first element end exchanges 2 bytes.
 * \param[in,out] x   Pointer to the first element.
 * \param[in]     n   The number of elements to swap.
 **************************************************************************/
void tom_io_nswap02(void *x, size_t n) {
    size_t i;
    uint16_t *xt = (uint16_t *)x;
    for (i=0; i<n; i++) {
        tom_io_swap02(&xt[i]);
    }
}
/***********************************************************************//**
 * \brief Swaps the bytes of 4-bytes elements of an array
 *
 * Gets a pointer to the first element end exchanges 4 bytes.
 * \param[in,out] x   Pointer to the first element.
 * \param[in]     n   The number of elements to swap.
 **************************************************************************/
void tom_io_nswap04(void *x, size_t n) {
    size_t i;
    uint32_t *xt = (uint32_t *)x;
    for (i=0; i<n; i++) {
        tom_io_swap04(&xt[i]);
    }
}
/***********************************************************************//**
 * \brief Swaps the bytes of 8-bytes elements of an array
 *
 * Gets a pointer to the first element end exchanges 8 bytes.
 * \param[in,out] x   Pointer to the first element.
 * \param[in]     n   The number of elements to swap.
 **************************************************************************/
void tom_io_nswap08(void *x, size_t n) {
    size_t i;
    uint64_t *xt = (uint64_t *)x;
    for (i=0; i<n; i++) {
        tom_io_swap08(&xt[i]);
    }
}
/***********************************************************************//**
 * \brief Swap array in memory of tom_io_em_complex type.
 *
 * Gets a pointer to the first element end exchanges the bytes.
 * \param[in,out] x   Pointer to the first element.
 * \param[in]     n   The number of elements to swap.
 **************************************************************************/
void tom_io_nswap_em_complex(void *x, size_t n) {
    size_t i;
    tom_io_em_complex *xt = (tom_io_em_complex *)x;
    for (i=0; i<n; i++) {
        tom_io_swap04(&xt[i].re);
        tom_io_swap04(&xt[i].im);
    }
}
/***********************************************************************//**
 * \brief Swap the bytes of an array with elements of a given size.
 *
 * Gets a pointer to the first element end exchanges the bytes.
 * \param[in,out]    x   Pointer to the first element.
 * \param[in]     size   Size of one element in bytes.
 * \param[in]        n   The number of elements to swap.
 **************************************************************************/
void tom_io_nswap(void *x, size_t size, size_t n) {
    switch (size) {
        case 0:
        case 1:
            break;
        case 2:
            tom_io_nswap02(x, n);
            break;
        case 4:
            tom_io_nswap04(x, n);
            break;
        case 8:
            tom_io_nswap08(x, n);
            break;
        default:
            {
                size_t i;
                for (i=0; i<n; i++) {
                    tom_io_swap(x, size);
                    x = (char *)x + size;
                }
            }
    }
}









/** \cond __HIDE_FUNCTIONS */
#define copy_vector(SRCTYPE, DSTTYPE, NAME)                                                         \
static void NAME(size_t        copy_vector_var_n,                                                   \
                 size_t        copy_vector_var_stridex,                                             \
                 const char  *copy_vector_var_src,                                                  \
                 char  *copy_vector_var_dst) {                                                      \
    size_t i, j;                                                                                    \
    for (i=0, j=0; i<copy_vector_var_n; i++, j+=copy_vector_var_stridex) {                          \
        ((DSTTYPE *)copy_vector_var_dst)[i] =                                                       \
                      (*((const SRCTYPE *)(copy_vector_var_src+j)));                                \
    }                                                                                               \
}
copy_vector(int8_t,    int8_t,      copy_vector_int8_t__int8_t)
copy_vector(int8_t,    int16_t,     copy_vector_int8_t__int16_t)
copy_vector(int8_t,    int32_t,     copy_vector_int8_t__int32_t)
copy_vector(int8_t,    float,       copy_vector_int8_t__float)
copy_vector(int8_t,    double,      copy_vector_int8_t__double)
copy_vector(int16_t,   int16_t,     copy_vector_int16_t__int16_t)
copy_vector(int16_t,   int32_t,     copy_vector_int16_t__int32_t)
copy_vector(int16_t,   float,       copy_vector_int16_t__float)
copy_vector(int16_t,   double,      copy_vector_int16_t__double)
copy_vector(int32_t,   int32_t,     copy_vector_int32_t__int32_t)
copy_vector(int32_t,   float,       copy_vector_int32_t__float)
copy_vector(int32_t,   double,      copy_vector_int32_t__double)
copy_vector(float,     float,       copy_vector_float__float)
copy_vector(float,     double,      copy_vector_float__double)
copy_vector(double,    float,       copy_vector_double__float)
copy_vector(double,    double,      copy_vector_double__double)
#define copy_vector_r2c(SRCTYPE, NAME)                                                                      \
static void NAME(size_t        copy_vector_var_n,                                                           \
                 size_t        copy_vector_var_stridex,                                                     \
                 const char  *copy_vector_var_src,                                                          \
                 char  *copy_vector_var_dst) {                                                              \
    size_t i, j;                                                                                            \
    for (i=0, j=0; i<copy_vector_var_n; i++, j+=copy_vector_var_stridex) {                                  \
        ((tom_io_em_complex *)copy_vector_var_dst)[i].re = (*((const SRCTYPE *)(copy_vector_var_src+j)));   \
        ((tom_io_em_complex *)copy_vector_var_dst)[i].im = 0;                                               \
    }                                                                                                       \
}
copy_vector(tom_io_em_complex, tom_io_em_complex, copy_vector_complex__complex)
copy_vector_r2c(int8_t,      copy_vector_int8_t__complex)
copy_vector_r2c(int16_t,     copy_vector_int16_t__complex)
copy_vector_r2c(int32_t,     copy_vector_int32_t__complex)
copy_vector_r2c(float,       copy_vector_float__complex)
copy_vector_r2c(double,      copy_vector_double__complex)
static void copy_vector_memcpy(size_t        copy_vector_var_n,
                               size_t        copy_vector_var_stridex,
                               const char  *copy_vector_var_src,
                               char  *copy_vector_var_dst) {
    /* Warning. in this case the variable copy_vector_var_n contains the number
       of bytes, not the number of elements. */
    memcpy(copy_vector_var_dst, copy_vector_var_src, copy_vector_var_n);
}

#define copy_vector_add(SRCTYPE, DSTTYPE, NAME)                                 \
static void NAME(size_t        copy_vector_var_n,                               \
                 size_t        copy_vector_var_stridex,                         \
                 size_t        copy_vector_var_binning,                         \
                 const char  *copy_vector_var_src,                              \
                 char  *copy_vector_var_dst) {                                  \
    size_t i, j;                                                                \
    for (i=0, j=0; i<copy_vector_var_n; i++, j+=copy_vector_var_stridex) {      \
        ((DSTTYPE *)copy_vector_var_dst)[i/copy_vector_var_binning] +=          \
            *((const SRCTYPE *)(copy_vector_var_src+j));                        \
    }                                                                           \
}
copy_vector_add(int8_t,   int64_t, copy_vector_add_int8_t__int64_t)
copy_vector_add(int8_t,   double,   copy_vector_add_int8_t__double)
copy_vector_add(int16_t,  int64_t, copy_vector_add_int16_t__int64_t)
copy_vector_add(int16_t,  double,   copy_vector_add_int16_t__double)
copy_vector_add(int32_t,  int64_t, copy_vector_add_int32_t__int64_t)
copy_vector_add(int32_t,  double,   copy_vector_add_int32_t__double)
copy_vector_add(float,    double,    copy_vector_add_float__double)
copy_vector_add(double,   double,    copy_vector_add_double__double)
/** \endcond*/





/***********************************************************************//**
 * \brief Uses fread to read data from file and swapes bytes.
 *
 * \param[out] ptr Pointer to memory (size*nmemb bytes).
 * \param[in]  size size in bytes of each element.
 * \param[in] nmemb Number of elements to read.
 * \param[in,out] stream File where to read from.
 **************************************************************************/
size_t tom_io_fread_swap(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    const size_t res = fread(ptr, size, nmemb, stream);
    if (res) { tom_io_nswap(ptr, size, res); }
    return res;
}
/** \copydoc tom_io_fread_swap */
size_t tom_io_fread_swap02(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    const size_t res = fread(ptr, 2, nmemb, stream);
    if (res) { tom_io_nswap02(ptr, res); }
    return res;
}
/** \copydoc tom_io_fread_swap */
size_t tom_io_fread_swap04(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    const size_t res = fread(ptr, 4, nmemb, stream);
    if (res) { tom_io_nswap04(ptr, res); }
    return res;
}
/** \copydoc tom_io_fread_swap */
size_t tom_io_fread_swap08(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    const size_t res = fread(ptr, 8, nmemb, stream);
    if (res) { tom_io_nswap08(ptr, res); }
    return res;
}
/** \copydoc tom_io_fread_swap */
size_t tom_io_fread_swap_em_complex(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    const size_t res = fread(ptr, sizeof(tom_io_em_complex), nmemb, stream);
    if (res) { tom_io_nswap_em_complex(ptr, res); }
    return res;
}




/***********************************************************************//**
 * \brief Returns the name of the microscope according to the header.
 *
 * \param[in] header The em-header.
 * \param[out] name The name of the microscope.
 * \param[in] bytes_reserved The length of the buffer name. The resulting
 *   name can be at most bytes_reserved-1 characters long.
 **************************************************************************/
void tom_io_em_get_microscope_name(const tom_io_em_header *header, char *name, size_t bytes_reserved) {
    switch (header->emdata[7]) {
        case 1:
            strncpy(name, "EM420          ", bytes_reserved-1);
            break;
        case 2:
            strncpy(name, "CM12           ", bytes_reserved-1);
            break;
        case 3:
            strncpy(name, "CM200          ", bytes_reserved-1);
            break;
        case 4:
            strncpy(name, "CM120/Biofilter", bytes_reserved-1);
            break;
        case 5:
            strncpy(name, "CM300          ", bytes_reserved-1);
            break;
        case 6:
            strncpy(name, "Polara         ", bytes_reserved-1);
            break;
        case 7:
            strncpy(name, "Titan          ", bytes_reserved-1);
            break;
        case 8:
            strncpy(name, "Tecnai F20     ", bytes_reserved-1);
            break;
        default:
            strncpy(name, "extern         ", bytes_reserved-1);
    }
    name[bytes_reserved-1] = 0;
}


/***********************************************************************//**
 * \brief Check whether the memory contains a valid em-header.
 *
 * Gets a pointer to data and checks whether it contains data as
 * defined by the header format. For that some fields (bytes) in the
 * memory are viewed. For that the size of the memory must be given,
 * so that the valid range is not exceeded.
 * \param[in] header is a pointer to the memory containing the header to be checked.
 * \param[in] size of the memory.
 * \return 0 if the format is not recognized and otherwise 1.
 **************************************************************************/
int tom_io_em_is_valid_header(const void *header, size_t size) {
    const tom_io_em_header *header_ = (const tom_io_em_header *)header;
    if (size < 512) {
        return 0;
    }
    switch (header_->machine) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
            break;
        default:
            return 0;
    }
    switch (header_->type) {
        case 1:
        case 2:
        case 4:
        case 5:
        case 8:
        case 9:
            break;
        default:
            return 0;
    }
    return 1;
}







/***********************************************************************//**
 * \brief Swaps the bytes of the em-header if it is necessary.
 *
 * Checks whether the header comes in an other encoding
 * (little vs. big endian) and if it does, swaps the fields. It does
 * not check the content of the header for consistency, except the part
 * needed to decide what to do.
 * \param [in,out] header is a pointer to the header.
 * \returns If the header contains data which does not allow to detect
 * the byte order, TOM_ERR_WRONG_HEADER is returned. Otherwise TOM_ERR_OK
 **************************************************************************/
int tom_io_em_swap_header(tom_io_em_header *header) {

    const int swap_bytes = tom_io_em_is_swaped(header);
    if (swap_bytes < 0) {
        return TOM_ERR_WRONG_HEADER;
    }
    if (swap_bytes) {
        /* Swap the bytes in the header. */
        int i;
        tom_io_swap(&header->dims[0], sizeof(header->dims[0]));
        tom_io_swap(&header->dims[1], sizeof(header->dims[1]));
        tom_io_swap(&header->dims[2], sizeof(header->dims[2]));
        for (i=0; i<40; i++) {
            tom_io_swap(&header->emdata[i], sizeof(header->emdata[i]));
        }
    }
    return TOM_ERR_OK;
}




/***********************************************************************//**
 * \brief Check whether the binary em-data is swaped.
 *
 * Looks into the header of the em-data and returns whether the data
 * is swaped compared to the current machine where the program runs.
 * This has to do with the byte order of the data: little vs. big endian.
 * \param[in] header A pointer to the em-header.
 * \return -1 if the header containes unexpected data, 1 if the data is
 *   in a different order than the current computer works with,
 *   and 0 it the data is in the same order.
 **************************************************************************/
int tom_io_em_is_swaped(const tom_io_em_header *header) {
    switch (header->machine) {
        case 1:
        case 2:
        case 6:
            return (tom_io_machine_byteorder() == TOM_IO_MACHINE_BYTEORDER_BIG_ENDIAN) ? 1 : 0;
        case 0:
        case 3:
        case 4:
        case 5:
            return (tom_io_machine_byteorder() == TOM_IO_MACHINE_BYTEORDER_LITTLE_ENDIAN) ? 1 : 0;
    }
    return -1;
}


/***********************************************************************//**
 * \brief Returns a identifier number for the data type to an em-header.
 *
 * Looks into an em-header and returns an integer constant which defines
 * the datatype. These numbers are declared in tom_defines.h and are named
 * \c TOM_IO_TYPE_XXX.
 * \param[in] header A pointer to the em-header.
 * \returns the numerical constant of the data-type or
 * \c TOM_IO_TYPE_VOID if the the type could not be detected
 * (i.e. the header is non valid).
 **************************************************************************/
int tom_io_em_get_iotype(const tom_io_em_header *header) {
    switch (header->type) {
        case 1:
            return TOM_IO_TYPE_INT8;
        case 2:
            return TOM_IO_TYPE_INT16;
        case 4:
            return TOM_IO_TYPE_INT32;
        case 5:
            return TOM_IO_TYPE_FLOAT;
        case 8:
            return TOM_IO_TYPE_COMPLEX4;
        case 9:
            return TOM_IO_TYPE_DOUBLE;
    }
    return TOM_IO_TYPE_UNKNOWN;
}


/***********************************************************************//**
 * \brief Sets the data-type of an em-header.
 *
 * Sets the datatype of an em-header to the given io-type. It is one
 * of the numerical constants defined in tom_io.h.
 * \param[out] header A pointer to the em-header.
 * \param[in]  iotype The io-type.
 * \return 1 if the iotype is supported by the em-format. Otherwise
 *   0 is returned.
 **************************************************************************/
int tom_io_em_set_iotype(tom_io_em_header *header, int iotype) {
    if (!header) {
        return 0;
    }
    switch (iotype) {
        case  TOM_IO_TYPE_INT8:
            header->type = 1;
            break;
        case TOM_IO_TYPE_INT16:
            header->type = 2;
            break;
        case TOM_IO_TYPE_INT32:
            header->type = 4;
            break;
        case TOM_IO_TYPE_FLOAT:
            header->type = 5;
            break;
        case TOM_IO_TYPE_COMPLEX4:
            header->type = 8;
            break;
        case TOM_IO_TYPE_DOUBLE:
            header->type = 9;
            break;
        default:
            return 0;
    }
    return 1;
}


/***********************************************************************//**
 * \brief Gets the size of the io-datatype in bytes.
 *
 * Gets the number of bytes needed to represent one element of a
 * given io-type. These types are numerical constants defined in
 * tom_io.h and named TOM_IO_TYPE_XXX.
 * \param[in] type the numerical representation of the type,
 *   i.e. one of the predefined constants.
 * \return the size of one data element in bytes. 0 if the type
 *   is not recognized.
 **************************************************************************/
int tom_io_iotype_datasize(int type) {
    switch (type) {
        case TOM_IO_TYPE_INT8:
            return 1;
        case TOM_IO_TYPE_INT16:
            return 2;
        case TOM_IO_TYPE_INT32:
            return 4;
        case TOM_IO_TYPE_FLOAT:
            return 4;
        case TOM_IO_TYPE_COMPLEX4:
            return 8;
        case TOM_IO_TYPE_DOUBLE:
            return 8;
        case TOM_IO_TYPE_INT64:
            return 8;
    }
    return 0;
}







/***********************************************************************//**
 * \brief Reads the header of an em-file.
 *
 * Reads the header of an em-file and returns the open filehandle.
 * \param[in]  filename is the name of the em-file to read.
 * \param[in]  mode The mode how to open the file. The file is opened
 *   calling fopen from ISO C89. This is its second parameter.
 *   NULL defaults to "rb", for reading in binary mode.
 * \param[out] header   memory to store the header of the file.
 *   This memory must already be preallocated. If the byte order
 *   differs, the function swapes the bytes of the header before
 *   returning it.
 * \param[out] fhandle  is a pointer to a (pointer to a) filehandle.
 *   If <tt>fhandle == NULL</tt> the file is closed after reading. Otherwise
 *   the filehandle is returned and must be closed later using fclose.
 *   The open file is seeked to the first position immediately after
 *   the em-header.

 * \return \c TOM_ERR_OK if the file could be successfully opened and a
 *   valid header could be read. If an error occures an error status is
 *   returned. In that case the content of \c tom_io_em_header and \c fhandle is
 *   undefined. However the file is always closed.
 **************************************************************************/
int tom_io_em_read_header(const char *filename, const char *mode, tom_io_em_header *header, FILE **fhandle) {

    FILE *f = NULL;
    int64_t fsize;
    int element_size;
    int fseek_success;

    if (!header || !filename) {
        return TOM_ERR_WRONG_INPUT;
    }
    if (!mode) {
        mode = "rb";
    }


    /* Open the file... */
    if (!(f = fopen(filename, mode))) {
        /*const int errorno = errno;*/
        return TOM_ERR_OPEN_FILE;
    }

    ___fseek_set(f, 0, fseek_success);

    assert(sizeof(tom_io_em_header) == 512);
    /* Does the compiler padd the struct wrong?
     * http://c-faq.com/struct/endpad.html
     * http://c-faq.com/struct/align.esr.html
     * http://c-faq.com/struct/padding.html
     * This should not happen since the largest type of tom_io_em_header is uint32 with
     * 4 bytes, and all the fields of tom_io_em_header are nultiples of 4.
     * Otherwise the reading/writing of the header should be changed to
     * http://c-faq.com/stdio/extconform.html
     */

    /* Read header and check it. */
    if (fread(header, sizeof(tom_io_em_header), 1, f) != 1) {
        fclose(f);
        return TOM_ERR_READ_FILE;
    }

    if (tom_io_em_swap_header(header)!=TOM_ERR_OK || !tom_io_em_is_valid_header(header, sizeof(tom_io_em_header))) {
        fclose(f);
        return TOM_ERR_WRONG_HEADER;
    }

    /* get file size */
    fsize = tom_io_getFileSize(f);
    element_size = tom_io_iotype_datasize(tom_io_em_get_iotype(header));

    /* Check the file size */
    if (fsize-(int64_t)sizeof(tom_io_em_header) != (int64_t)element_size* (int64_t)header->dims[0] * header->dims[1]*header->dims[2]) {
        fclose(f);
        return TOM_ERR_FILESIZE_MISMATCH;
    }

    if (fhandle) {
        *fhandle = f;
    } else {
        fclose(f);
    }
    return TOM_ERR_OK;

}






/***********************************************************************//**
 * \brief Returns the resulting size of a volume after binning and sampling.
 *
 * Takes the size of a volume and parameters specifying a subregion, sampling
 * factor and binning. It returns the size of the volume after sampling and
 * after both binning and sampling.
 * \param[in] voldims a 3 vector with the size of the volume along the
 *   directions x, y, and z respectively. This parameter can not be
 *   ommitted. All other parameters can be set to \c NULL.
 * \param[in] subregion_ The region of interest of the volume.
 *   The positions <tt>subregion_[0,1,2]</tt> is the number of the "first" corner in
 *   the volume. The numbering is zero based. The parameters
 *   <tt>subregion_[3,4.5]</tt> contain the size of the subregion. Passing \c subregion_ as
 *   \c NULL is the same as selecting the whole volume. The subregion must be
 *   entirely inside the dimensions of the volume.
 * \param[out] subregion_out Returns the selected subregion. If it is \c NULL,
 *   nothing happens. Otherwise the content of the valid \c subregion_ is
 *   copied into this already preallocated 6 vector. If \c subregion_ is set
 *   to \c NULL, \c subregion_out contains the whole volume dimension.
 * \param[in] binning_ can either be \c NULL or a 3 vector with the binning
 *   factors along each direction. \c NULL means no binning. Otherwise this
 *   specifies how many data elements along each direction are combined.
 *   A value of 0 or 1 means no binning along that direction.
 * \param[out] binning_out Must be a 3 vector or \c NULL. If given, the binning
 *   factors are copied there. i.e. either the values are copied directly from \c binning_
 *   or set to 1 (= no binning).
 * \param[in] sampling_ similar to \c binning_ this is the sampling factor.
 * \param[out] sampling_out similar to \c binning_out the sampling factor is copied here.
 * \param[out] dims_sampled_out If non set to \c NULL, the resulting dimension
 *   after sampling of the subregion is copied here.
 * \param[out] dims_out Contains the size of the resulting volume after
 *   selecting the region of interest (ROI), sampling and binning.
 * \return If the subregion is inside the volume and the resulting volume
 *   (\c dims_out) after binning is not empty, \c TOM_ERR_OK is returned. If a
 *   different error-value is returned none of the output parameters are changed.
 *
 * Binning and sampling can be combined. In that case the volume is sampled
 * first. That means every <tt>sampling_[i]</tt>th value from the subregion is taken,
 * beginning with the first corner. If the dimension of the volume is not a
 * multiple of the sampling factor, the last reminding elements are also taken into
 * account. Thus the sampled volume contains always the first corner.
 * This is different for binning: there only the whole multiples of the
 * binning factor are taken. By binning (and optional
 * sampling) the volume can became empty.
 **************************************************************************/
int tom_io_calculate_sizes(const uint32_t *voldims,
                    const uint32_t *subregion_, uint32_t *subregion_out,
                    const uint32_t *sampling_, uint32_t *sampling_out,
                    const uint32_t *binning_, uint32_t *binning_out,
                    uint32_t *dims_sampled_out,
                    uint32_t *dims_out) {

    uint32_t binning[3] =   { 1, 1, 1};
    uint32_t sampling[3] =  { 1, 1, 1};
    uint32_t subregion[6] = { 0, 0, 0, 0, 0, 0 };
    uint32_t dims_sampled[3];
    uint32_t dims[3];
    int i;

    if (!voldims || !voldims[0] || !voldims[1] || !voldims[2]) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (subregion_) {
        for (i=0; i<6; i++) {
            subregion[i] = subregion_[i];
        }
    } else {
        subregion[3] = voldims[0];
        subregion[4] = voldims[1];
        subregion[5] = voldims[2];
    }

    if (binning_) {
        binning[0] = binning_[0] ? binning_[0] : 1;
        binning[1] = binning_[1] ? binning_[1] : 1;
        binning[2] = binning_[2] ? binning_[2] : 1;
    }
    if (sampling_) {
        sampling[0] = sampling_[0] ? sampling_[0] : 1;
        sampling[1] = sampling_[1] ? sampling_[1] : 1;
        sampling[2] = sampling_[2] ? sampling_[2] : 1;
    }

    if (subregion[3]<=0 || (((size_t)subregion[0]+(size_t)subregion[3])>voldims[0]) ||
        subregion[4]<=0 || (((size_t)subregion[1]+(size_t)subregion[4])>voldims[1]) ||
        subregion[5]<=0 || (((size_t)subregion[2]+(size_t)subregion[5])>voldims[2])) {
        return TOM_ERR_SUBREGION_OUT;
    }


    dims[0] = dims_sampled[0] = (subregion[3]+sampling[0]-1) / sampling[0];
    dims[1] = dims_sampled[1] = (subregion[4]+sampling[1]-1) / sampling[1];
    dims[2] = dims_sampled[2] = (subregion[5]+sampling[2]-1) / sampling[2];


    if (binning[0]>1 || binning[1]>1 || binning[2]>1) {
        dims_sampled[0] = dims[0] - dims[0]%binning[0];
        dims_sampled[1] = dims[1] - dims[1]%binning[1];
        dims_sampled[2] = dims[2] - dims[2]%binning[2];
        if (!dims_sampled[0] || !dims_sampled[1] || !dims_sampled[2]) {
            return TOM_ERR_BINNING_TOO_HIGH;
        }
        dims[0] = dims_sampled[0] / binning[0];
        dims[1] = dims_sampled[1] / binning[1];
        dims[2] = dims_sampled[2] / binning[2];
    }


    /* Copy the local copies to the output arrays. */
    if (subregion_out)    { for (i=0; i<6; i++) { subregion_out[i]    = subregion[i];    } }
    if (binning_out)      { for (i=0; i<3; i++) { binning_out[i]      = binning[i];      } }
    if (sampling_out)     { for (i=0; i<3; i++) { sampling_out[i]     = sampling[i];     } }
    if (dims_sampled_out) { for (i=0; i<3; i++) { dims_sampled_out[i] = dims_sampled[i]; } }
    if (dims_out)         { for (i=0; i<3; i++) { dims_out[i]         = dims[i];         } }


    return TOM_ERR_OK;
}





/****************************************************************************//**
 * \brief Check the size of stride parameters and set in case of default value.
 *
 * If the stride parameters are NULL they are set to their minimum (default)
 * If 0 is returned, *stride stayes unchanged.
 *******************************************************************************/
int tom_io_check_stride(size_t size_of_type, size_t sizex, size_t sizey, size_t sizez,
                        size_t *stridex, size_t *stridey, size_t *stridez) {

    if (size_of_type<1 || sizex<1 || sizey<1 || sizez<1 || !stridex || !stridey || !stridez) {
        return 0;
    }

    {
        size_t sx = *stridex;
        size_t sy = *stridey;
        size_t sz = *stridez;
        int res;

        /* Set default values of the stride parameter. */
        if (!sx) { sx = size_of_type;   }
        if (!sy) { sy = sizex * sx;     }
        if (!sz) { sz = sizey * sy;     }

        res = sx>=size_of_type && sy>=sizex*sx && sz>=sizey*sy;

        if (res) {
            *stridex = sx;
            *stridey = sy;
            *stridez = sz;
        }
        return res;
    }
}




/***********************************************************************//**
 * \brief reads rawdata from an open file.
 *
 * Reads a subregion (or the entire volume) out of an open file.
 * Maybe you want to use the function tom_io_read_vol_sac which is a wrapper for
 * this one, and allocates the memory by itself. Use this version if you
 * want to allocate the memory on your own. You can use tom_io_calculate_sizes to
 * know the resulting amount of memory to allocate according to the same
 * parameters.
 * \param[in]  f   a handle to the file where to read from. It must be
 *   open and seeked to the position of the first data element of the volume
 *   (not of the subregion!). The file is never closed in this function, but
 *   the position of the filepointer is undefined afterwards.
 * \param[in]  voldims     The size of the volume contained in the file.
 *   The elements along dimension 0 are stored one after another, the dimension
 *   with size voldims[2] is the slowest one.
 * \param[in]  iotype  An identifier of the data type contained in the file.
 *   It is one of the numerical constants defined in tom_io.h and is named
 *   something like TOM_IO_TYPE_XXX.
 * \param[in]  swaped  If true, the data elements have the be swaped after
 *   reading (little vs. big endian).
 * \param[in]  subregion   is the region of interest of the volume. It must
 *   be a 6 vector with the first corner and the dimension of the subregion.
 *   See also the parameter subregion_out of tom_io_calculate_sizes.
 * \param[in]  sampling    is the sampling rate. See the parameter subregion_out
 *   of tom_io_calculate_sizes.
 * \param[in]  binning     analog to sampling it is the binning factor.
 * \param[out] pdata   A pointer to the already allocated memory. The resulting
 *   volume is saved there. Of course it is never freed by this function. Its
 *   size can be determined before using the function tom_io_calculate_sizes and
 *   corresponds the output parameter dims_out of these function (multiplied by
 *   sizeof(restype)).
 * \param[in]  restype     The output data can be converted from the datatype
 *   iotype to this one. pdata must be saved as large to hold data of this type
 *   (again it is one of the numerical constants TOM_IO_TYPE_XXX). Specially
 *   for binning integer data, here must be an floating point type. Conversions
 *   to a smaller data range (for example TOM_IO_TYPE_INT16 to TOM_IO_TYPE_INT8)
 *   is not implmented because it is not clear how to handle loss of precision.
 * \return TOM_ERR_OK if no error occured or an error status otherwise.
 **************************************************************************/
int tom_io_read_vol(FILE *f, const uint32_t voldims[], int iotype, int swaped,
                 const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning,
                 void *pdata, int restype) {

    uint32_t x, y, z;
    const uint32_t element_size = tom_io_iotype_datasize(iotype);
    const uint32_t reselement_size = tom_io_iotype_datasize(restype);
    uint32_t element_size_tmp;
    void *pdata_tmp;
    unsigned char *pdatac;
    size_t numel;
    uint32_t dims[3];
    uint32_t dims_sampled[3];
    uint32_t buffer_length;
    char *buffer = NULL;
    uint32_t subregion_[6];
    uint32_t sampling_[3];
    uint32_t binning_[3];
    int fseek_success;

    size_t (*fcn_fread)(void *ptr, size_t size, size_t nmemb, FILE *stream) = fread;

    /* Check input parameters whether they are consistent. */
    if (!f || !voldims || !subregion || !sampling || !binning || !pdata) {
        return TOM_ERR_WRONG_INPUT;
    }
    switch (iotype) {
        case TOM_IO_TYPE_INT8:
        case TOM_IO_TYPE_INT16:
        case TOM_IO_TYPE_INT32:
        case TOM_IO_TYPE_FLOAT:
        case TOM_IO_TYPE_COMPLEX4:
        case TOM_IO_TYPE_DOUBLE:
            break;
        default:
            return TOM_ERR_IOTYPE_NOT_SUPPORTED;
    }
    switch (restype) {
        case TOM_IO_TYPE_INT8:
        case TOM_IO_TYPE_INT16:
        case TOM_IO_TYPE_INT32:
        case TOM_IO_TYPE_FLOAT:
        case TOM_IO_TYPE_COMPLEX4:
        case TOM_IO_TYPE_DOUBLE:
            break;
        default:
            return TOM_ERR_IOTYPE_NOT_SUPPORTED;
    }


    /* Get the subregion and check the dimension. */
    if ((x=tom_io_calculate_sizes(voldims, subregion, subregion_, sampling, sampling_, binning, binning_, dims_sampled, dims)) != TOM_ERR_OK) {
        return x;
    }


    /* Switch input parameter to clean parameters. */
    subregion = subregion_; sampling = sampling_; binning = binning_;

    numel = (size_t)dims[0] * (size_t)dims[1] * (size_t)dims[2];
    buffer_length = (dims_sampled[0]-1)*sampling[0] + 1;

    if (binning[0]>1 || binning[1]>1 || binning[2]>1) {
        void (*copy_vector_fcn)(size_t, size_t, size_t, const char*, char*) = NULL;

        #if 0 && defined __LCC__ && defined _WIN32
        const double binning3 = (double)binning[0]*(double)binning[1]*(double)binning[2];
        #else
        const uint64_t binning3 = (uint64_t)binning[0]*(uint64_t)binning[1]*(uint64_t)binning[2];
        #endif
        int iotype_tmp = TOM_IO_TYPE_UNKNOWN;

        /*_printf(("READ: binning %s\n", (swaped?"with swap":"no swap")));*/
        if (iotype == TOM_IO_TYPE_COMPLEX4) {
            /* Complex not yet implemented. Maybe it is not meaningfull to do binning of complex type?!?! */
            return TOM_ERR_NO_COMPLEX_BINNING;
        }

        /* Set the line sample function... */
        switch (iotype) {
            case TOM_IO_TYPE_INT8:
                if (binning3 < 72057594037927936.) { /* <= 0x100000000000000 = 72057594037927936 = (2^63)/(2^7)*/
                    iotype_tmp = TOM_IO_TYPE_INT64;
                    copy_vector_fcn = copy_vector_add_int8_t__int64_t;
                } else {
                    iotype_tmp = TOM_IO_TYPE_DOUBLE;
                    copy_vector_fcn = copy_vector_add_int8_t__double;
                }
                break;
            case TOM_IO_TYPE_INT16:
                if (binning3 <= 281474976710656.) { /* <= 0x1000000000000 = 281474976710656 = (2^63)/(2^15)*/
                    iotype_tmp = TOM_IO_TYPE_INT64;
                    copy_vector_fcn = copy_vector_add_int16_t__int64_t;
                } else {
                    iotype_tmp = TOM_IO_TYPE_DOUBLE;
                    copy_vector_fcn = copy_vector_add_int16_t__double;
                }
                if (swaped) { fcn_fread = tom_io_fread_swap02; }
                break;
            case TOM_IO_TYPE_INT32:
                if (binning3 <= 4294967296.) { /* <= 0x100000000 = 4294967296 = (2^63)/(2^31)*/
                    iotype_tmp = TOM_IO_TYPE_INT64;
                    copy_vector_fcn = copy_vector_add_int32_t__int64_t;
                } else {
                    iotype_tmp = TOM_IO_TYPE_DOUBLE;
                    copy_vector_fcn = copy_vector_add_int32_t__double;
                }
                if (swaped) { fcn_fread = tom_io_fread_swap04; }
                break;
            case TOM_IO_TYPE_FLOAT:
                iotype_tmp = TOM_IO_TYPE_DOUBLE;
                copy_vector_fcn = copy_vector_add_float__double;
                if (swaped) { fcn_fread = tom_io_fread_swap04; }
                break;
            case TOM_IO_TYPE_DOUBLE:
                iotype_tmp = TOM_IO_TYPE_DOUBLE;
                copy_vector_fcn = copy_vector_add_double__double;
                if (swaped) { fcn_fread = tom_io_fread_swap08; }
                break;
            case TOM_IO_TYPE_COMPLEX4:                                                     break; /* Not expected. */
            default:                                                                        break; /* Not expected. */
        }

        /* In case of binnings of integers, a temporary buffer of a larger integer type
         * is used to avoid overflows. In case of floating point numbers, the destination
         * buffer is used. */
        if (iotype_tmp == restype) {
            pdata_tmp = pdata;
            element_size_tmp = tom_io_iotype_datasize(iotype_tmp);
            /* Set the sum-buffer to 0. */
            memset(pdata_tmp, 0, numel*element_size_tmp);
        } else {
            /* Sum-buffer is initialized with 0 by calling calloc. */
            element_size_tmp = tom_io_iotype_datasize(iotype_tmp);
            if (!(pdata_tmp = calloc(numel, element_size_tmp))) {
                return TOM_ERR_MALLOC;
            }
        }

        /* A buffer for reading line per line (i.e. along the x-direction). */
        buffer = (char *)malloc(buffer_length * element_size);

        /* Read the elements line per line (along x direction) and sum them up with the copy_vector_fcn. */

        {
            const int64_t fseek_x_width = (int64_t)((int64_t)voldims[0]*(int64_t)sampling[1] - buffer_length) * element_size;
            const int64_t fseek_y_width = (int64_t)((((int64_t)voldims[1]*(int64_t)sampling[2] - (int64_t)dims_sampled[1]*sampling[1]) * voldims[0])) * (int64_t)element_size;
            const uint32_t dims0_bytes = dims[0]*element_size_tmp;
            const uint32_t dims1_bytes = dims0_bytes*dims[1];
            uint32_t binning_1_mod = binning[1];
            uint32_t binning_2_mod = binning[2];
            int fseek_success;

            /* Set the variables for the sampling function. */
            const size_t copy_vector_var_n = dims_sampled[0];
            const size_t copy_vector_var_stridex = sampling[0] * element_size;
            const size_t copy_vector_var_binning = binning[0];
            const char *const copy_vector_var_src = buffer;
            char *copy_vector_var_dst = (char *)pdata_tmp;
            char *copy_vector_var_dst_start = copy_vector_var_dst;

            /* seek to the first element to read... */
            ___fseek_cur(f, (((int64_t)subregion[2]*(int64_t)voldims[1] + (int64_t)subregion[1])*voldims[0] + (int64_t)subregion[0]) * (int64_t)element_size, fseek_success);
            if (!fseek_success) {
                if (pdata != pdata_tmp) {
                    free(pdata_tmp);
                }
                free(buffer);
                return TOM_ERR_READ_FILE;
            }

            /*{long long o; ___ftell(f, o); printf("%d: seek to %llu\n", __LINE__, (long long)o); }*/

            for (z=0; z<dims_sampled[2]; z++) {
                for (y=0; y<dims_sampled[1]; y++) {
                    /* Read the current line into to buffer. */
                    if (fcn_fread(buffer, element_size, buffer_length, f) != buffer_length) {
                        if (pdata != pdata_tmp) {
                            free(pdata_tmp);
                        }
                        free(buffer);
                        return TOM_ERR_READ_FILE;
                    }
                    /* Sample the current line and sum it up into the destination buffer copy_vector_var_dst. */

                    copy_vector_fcn(copy_vector_var_n, copy_vector_var_stridex, copy_vector_var_binning, copy_vector_var_src, copy_vector_var_dst);

                    /* Choose the next line in the destination buffer if enough bytes have been sampled... */
                    if (--binning_1_mod < 1) {
                        binning_1_mod = binning[1];
                        copy_vector_var_dst += dims0_bytes;
                    }

                    ___fseek_cur(f, fseek_x_width, fseek_success); /* Seek to next data position. */
                }
                /* Reset destination buffer or step to next "heigth" in the volume. */
                if (--binning_2_mod < 1) {
                    binning_2_mod = binning[2];
                    copy_vector_var_dst_start += dims1_bytes;
                }
                copy_vector_var_dst = copy_vector_var_dst_start;

                ___fseek_cur(f, fseek_y_width, fseek_success); /* Seek to next data position. */
            }
        }
        free(buffer);

        /* Divide every element of the sum volume to get the mean (binned value). */
        #define __ROUND(x) (floor(((double)x)+0.5))
        {
            #if defined __LCC__ && defined _WIN32
            size_t i;
            #else
            uint64_t i;
            #endif
            const double binning3d = binning3;
            if (pdata==pdata_tmp) {
                switch (restype) {
                    case TOM_IO_TYPE_FLOAT:
                        for (i=0; i<numel; i++) {
                            ((float *)pdata)[i] = (float) ((((float *)pdata)[i]) / binning3d);
                        }
                        break;
                    case TOM_IO_TYPE_DOUBLE:
                        for (i=0; i<numel; i++) {
                            ((double *)pdata)[i] = ((double *)pdata)[i] / binning3d;
                        }
                        break;
                    default:                            break; /* Not expected. */
                }
            } else {
                switch (restype) {
                    case TOM_IO_TYPE_INT8:
                        switch (iotype_tmp) {
                            case TOM_IO_TYPE_INT64:
                                for (i=0; i<numel; i++) {
                                    ((int8_t *)pdata)[i] = (int8_t) __ROUND( ((double)(((int64_t *)pdata_tmp)[i])) / binning3d );
                                }
                                break;
                            case TOM_IO_TYPE_DOUBLE:
                                for (i=0; i<numel; i++) {
                                    ((int8_t *)pdata)[i] = (int8_t) __ROUND( (((double *)pdata_tmp)[i]) / binning3d );
                                }
                                break;
                            default:                    break; /* Not expected. */
                        }
                        break;
                    case TOM_IO_TYPE_INT16:
                        switch (iotype_tmp) {
                            case TOM_IO_TYPE_INT64:
                                for (i=0; i<numel; i++) {
                                    ((int16_t *)pdata)[i] = (int16_t) __ROUND( ((double)(((int64_t *)pdata_tmp)[i])) / binning3d );
                                }
                                break;
                            case TOM_IO_TYPE_DOUBLE:
                                for (i=0; i<numel; i++) {
                                    ((int16_t *)pdata)[i] = (int16_t) __ROUND( (((double *)pdata_tmp)[i]) / binning3d );
                                }
                                break;
                            default:                    break; /* Not expected. */
                        }
                        break;
                    case TOM_IO_TYPE_INT32:
                        switch (iotype_tmp) {
                            case TOM_IO_TYPE_INT64:
                                for (i=0; i<numel; i++) {
                                    ((int32_t *)pdata)[i] = (int32_t) __ROUND( ((double)(((int64_t *)pdata_tmp)[i])) / binning3d );
                                }
                                break;
                            case TOM_IO_TYPE_DOUBLE:
                                for (i=0; i<numel; i++) {
                                    ((int32_t *)pdata)[i] = (int32_t) __ROUND( (((double *)pdata_tmp)[i]) / binning3d );
                                }
                                break;
                            default:                    break; /* Not expected. */
                        }
                        break;
                    case TOM_IO_TYPE_FLOAT:
                        switch (iotype_tmp) {
                            case TOM_IO_TYPE_INT64:
                                for (i=0; i<numel; i++) {
                                    ((float *)pdata)[i] = (float)(((double)(((int64_t *)pdata_tmp)[i])) / binning3d);
                                }
                                break;
                            case TOM_IO_TYPE_DOUBLE:
                                for (i=0; i<numel; i++) {
                                    ((float *)pdata)[i] = (float)((((double *)pdata_tmp)[i]) / binning3d);
                                }
                                break;
                            default:                    break; /* Not expected. */
                        }
                        break;
                    case TOM_IO_TYPE_DOUBLE:
                        switch (iotype_tmp) {
                            case TOM_IO_TYPE_INT64:
                                for (i=0; i<numel; i++) {
                                    ((double *)pdata)[i] = ((double)(((int64_t *)pdata_tmp)[i])) / binning3d;
                                }
                                break;
                            case TOM_IO_TYPE_DOUBLE:
                                for (i=0; i<numel; i++) {
                                    ((double *)pdata)[i] = (((double *)pdata_tmp)[i]) / binning3d;
                                }
                                break;
                            default:                    break; /* Not expected. */
                        }
                        break;
                    case TOM_IO_TYPE_COMPLEX4:         break; /* Not expected. */
                    default:                            break; /* Not expected. */
                }
            }
            if (pdata != pdata_tmp) {
                free(pdata_tmp);
            }
        }
        #undef __ROUND

    } else {
        size_t bytes_read;

        if (iotype==restype &&
            sampling[0]==1 && sampling[1]==1 && sampling[2]==1 &&
            subregion[3]==voldims[0] && subregion[4]==voldims[1] && subregion[5]==voldims[2]) {
            /* Read the whole volume without type conversion and without swaping data. */

            #if defined _WIN64
                #define __READ_ALL_AT_ONCE_AMOUNT 0xFFFFFF
            #else
                #define __READ_ALL_AT_ONCE_AMOUNT 0xFFFFFFFF
            #endif
            size_t remainder = __READ_ALL_AT_ONCE_AMOUNT / element_size;
            /*_printf(("READ: all at once\n"));*/

            if (swaped && element_size>1) {
                if (iotype == TOM_IO_TYPE_COMPLEX4) {
                    fcn_fread = tom_io_fread_swap_em_complex;
                } else {
                    fcn_fread = tom_io_fread_swap;
                }
            }

            pdatac = (unsigned char *)pdata;
            while (numel > 0) {
                if (numel > remainder) {
                    numel -= remainder;
                } else {
                    remainder = numel;
                    numel = 0;
                }

                bytes_read = fcn_fread(pdatac, element_size, remainder, f);
                if (bytes_read != remainder) {
                    if (feof(f)) {
                        return TOM_ERR_FILESIZE_MISMATCH;
                    } else {
                        return TOM_ERR_READ_FILE;
                    }
                }
                pdatac += element_size*remainder;
            }
            #undef __READ_ALL_AT_ONCE_AMOUNT

        } else {
            void (*copy_vector_fcn)(size_t, size_t, const char*, char*) = NULL;

            const uint32_t dims0_bytes = dims[0]*reselement_size;
            const int64_t fseek_x_width = (int64_t)((int64_t)voldims[0]*sampling[1] - (int64_t)buffer_length) * (int64_t)element_size;
            const int64_t fseek_y_width = (int64_t)((((int64_t)voldims[1]*(int64_t)sampling[2] - (int64_t)dims[1]*sampling[1]) * voldims[0])) * (int64_t)element_size;

            /* seek to the first element to read... */
            ___fseek_cur(f, (((int64_t)subregion[2]*(int64_t)voldims[1] + (int64_t)subregion[1])*voldims[0] + (int64_t)subregion[0]) * (int64_t)element_size, fseek_success);
            /*{long long o; ___ftell(f, o); printf("%d: seek to %llu\n", __LINE__, (long long)o); }*/

            if (sampling[0] != 1 || restype!=iotype) {

                size_t copy_vector_var_n;
                size_t copy_vector_var_stridex;
                const char *copy_vector_var_src;
                char *copy_vector_var_dst;

                /* Sample also along the x-deriction... */
                /*_printf(("READ: sample %s\n", swaped?"swaped":"noswaped"));*/

                /* Function to sample the values from the read buffer (specific to each datatype). */
                switch (iotype) {
                    case TOM_IO_TYPE_INT8:
                        switch (restype) {
                            case TOM_IO_TYPE_INT8:     copy_vector_fcn = copy_vector_int8_t__int8_t;            break;
                            case TOM_IO_TYPE_INT16:    copy_vector_fcn = copy_vector_int8_t__int16_t;           break;
                            case TOM_IO_TYPE_INT32:    copy_vector_fcn = copy_vector_int8_t__int32_t;           break;
                            case TOM_IO_TYPE_FLOAT:    copy_vector_fcn = copy_vector_int8_t__float;             break;
                            case TOM_IO_TYPE_DOUBLE:   copy_vector_fcn = copy_vector_int8_t__double;            break;
                            default:
                                return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                        }
                        break;
                    case TOM_IO_TYPE_INT16:
                        switch (restype) {
                            case TOM_IO_TYPE_INT16:    copy_vector_fcn = copy_vector_int16_t__int16_t;          break;
                            case TOM_IO_TYPE_INT32:    copy_vector_fcn = copy_vector_int16_t__int32_t;          break;
                            case TOM_IO_TYPE_FLOAT:    copy_vector_fcn = copy_vector_int16_t__float;            break;
                            case TOM_IO_TYPE_DOUBLE:   copy_vector_fcn = copy_vector_int16_t__double;           break;
                            default:
                                return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                        }
                        if (swaped) { fcn_fread = tom_io_fread_swap02; }
                        break;
                    case TOM_IO_TYPE_INT32:
                        switch (restype) {
                            case TOM_IO_TYPE_INT32:    copy_vector_fcn = copy_vector_int32_t__int32_t;          break;
                            case TOM_IO_TYPE_FLOAT:    copy_vector_fcn = copy_vector_int32_t__float;            break;
                            case TOM_IO_TYPE_DOUBLE:   copy_vector_fcn = copy_vector_int32_t__double;           break;
                            default:
                                return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                        }
                        if (swaped) { fcn_fread = tom_io_fread_swap04; }
                        break;
                    case TOM_IO_TYPE_FLOAT:
                        switch (restype) {
                            case TOM_IO_TYPE_FLOAT:    copy_vector_fcn = copy_vector_float__float;              break;
                            case TOM_IO_TYPE_DOUBLE:   copy_vector_fcn = copy_vector_float__double;             break;
                            default:
                                return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                        }
                        if (swaped) { fcn_fread = tom_io_fread_swap04; }
                        break;
                    case TOM_IO_TYPE_COMPLEX4:
                        switch (restype) {
                            case TOM_IO_TYPE_COMPLEX4: copy_vector_fcn = copy_vector_complex__complex;          break;
                            default:
                                return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                        }
                        if (swaped) { fcn_fread = tom_io_fread_swap_em_complex; }
                        break;
                    case TOM_IO_TYPE_DOUBLE:
                        switch (restype) {
                            case TOM_IO_TYPE_FLOAT:    copy_vector_fcn = copy_vector_double__float;             break;
                            case TOM_IO_TYPE_DOUBLE:   copy_vector_fcn = copy_vector_double__double;            break;
                            default:
                                return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                        }
                        if (swaped) { fcn_fread = tom_io_fread_swap08; }
                        break;
                    default:
                        return TOM_ERR_IOTYPE_NOT_SUPPORTED;
                }
                buffer = (char *)malloc(buffer_length * element_size);
                copy_vector_var_n = dims[0];
                copy_vector_var_stridex = sampling[0]*element_size;
                copy_vector_var_src = buffer;
                copy_vector_var_dst = (char *)pdata;
                for (z=0; z<dims[2]; z++) {
                    for (y=0; y<dims[1]; y++) {
                        /* Read the current line into to buffer. */
                        bytes_read = fcn_fread(buffer, element_size, buffer_length, f);
                        if (bytes_read != buffer_length) {
                            free(buffer);
                            if (feof(f)) {
                                return TOM_ERR_FILESIZE_MISMATCH;
                            } else {
                                return TOM_ERR_READ_FILE;
                            }
                        }
                        copy_vector_fcn(copy_vector_var_n, copy_vector_var_stridex, copy_vector_var_src, copy_vector_var_dst);

                        copy_vector_var_dst += dims0_bytes;
                        ___fseek_cur(f, fseek_x_width, fseek_success);
                    }
                    ___fseek_cur(f, fseek_y_width, fseek_success);
                }
                free(buffer);
            } else {
                switch (iotype) {
                    case TOM_IO_TYPE_INT8:
                        break;
                    case TOM_IO_TYPE_INT16:
                        if (swaped) { fcn_fread = tom_io_fread_swap02; }
                        break;
                    case TOM_IO_TYPE_INT32:
                        if (swaped) { fcn_fread = tom_io_fread_swap04; }
                        break;
                    case TOM_IO_TYPE_FLOAT:
                        if (swaped) { fcn_fread = tom_io_fread_swap04; }
                        break;
                    case TOM_IO_TYPE_COMPLEX4:
                        if (swaped) { fcn_fread = tom_io_fread_swap_em_complex; }
                        break;
                    case TOM_IO_TYPE_DOUBLE:
                        if (swaped) { fcn_fread = tom_io_fread_swap08; }
                        break;
                    default:
                        return TOM_ERR_IOTYPE_NOT_SUPPORTED;
                }


                pdatac = (unsigned char *)pdata;
                /* Read along the x-direction without gaps and without swaping and without type conversion. */

                /*_printf(("READ: sample whole line (%sswaped)\n", swaped ? "" : "not "));*/

                for (z=0; z<dims[2]; z++) {
                    for (y=0; y<dims[1]; y++) {
                        /* No sampling along x. Read directly into memory. */
                        bytes_read = fcn_fread(pdatac, element_size, buffer_length, f);
                        if (bytes_read != buffer_length) {
                            if (feof(f)) {
                                return TOM_ERR_FILESIZE_MISMATCH;
                            } else {
                                return TOM_ERR_READ_FILE;
                            }
                        }
                        pdatac += dims0_bytes;
                        ___fseek_cur(f, fseek_x_width, fseek_success);
                    }
                    ___fseek_cur(f, fseek_y_width, fseek_success);
                }
            }
        }
    }
    return TOM_ERR_OK;
}




/***********************************************************************//**
 * \brief Reads a volume from an open file.
 *
 * This is a wrapper for tom_io_read_vol (stand alone complex :) It allocates the memory on its own using malloc.
 * \param[in]  f see the corresponding parameter in tom_io_read_vol.
 * \param[in]  voldims see the corresponding parameter in tom_io_read_vol.
 * \param[in]  iotype see the corresponding parameter in tom_io_read_vol.
 * \param[in]  swaped see the corresponding parameter in tom_io_read_vol.
 * \param[in]  subregion_ see the corresponding parameter in tom_io_read_vol.
 * \param[in]  sampling_ see the corresponding parameter in tom_io_read_vol.
 * \param[in]  binning_ see the corresponding parameter in tom_io_read_vol.
 * \param[out] data in contrast to tom_io_read_vol this is a pointer to a pointer
 *   to the volume. The memory is allocated by this function and must be
 *   freed later calling free.
 * \param[out] dims_ is a 3 vector where the dimension of the result is
 *   returned. This parameter is not present in tom_io_read_vol because there the
 *   caller must know the resulting size in before.
 * \param[out] restype_ Returns the datatype. In most cases it is the same
 *   as the input parameter iotype. Only in case of binning of integer data,
 *   a floating point type is used.
 * \param[in]  fcn_malloc The function with which the memory is allocated.
 *   Passing NULL, defaults to malloc from ISO C89.
 * \param[in]  fcn_free The free function to release the memory which was allocated
 *   using fcn_malloc. It is only called in case of an exception to clean up.
 *   Passing NULL, defaults to free from ISO C89. Leaving fcn_malloc undefined,
 *   but defining fcn_free results in an error and vice versa.
 * \return TOM_ERR_OK if no error occures. The open file f is never
 *   closed however. In case of error the content of the output parameters
 *   is undefined.
 *
 * The momory is allocated using the function fcn_malloc (and,
 * in case of an error during reading, freed with tom_fcn_vol_free).
 * These function-pointers default to malloc and free from the standard
 * library.
 **************************************************************************/
int tom_io_read_vol_sac(FILE *f, const uint32_t voldims[], int iotype, int swaped,
                 const uint32_t *subregion_, const uint32_t *sampling_, const uint32_t *binning_,
                 void **data, uint32_t *dims_, int *restype_,
                 void *(*fcn_malloc)(size_t size),
                 void  (*fcn_free)(void *ptr)) {

    void *pdata;
    int i;

    uint32_t  binning[3] = {1, 1, 1};
    uint32_t sampling[3] = {1, 1, 1};
    uint32_t subregion[6];
    uint32_t dims[3];
    uint32_t dims_sampled[3];
    int restype;

    if (!f || !voldims || !data || !restype_ || (!fcn_malloc && fcn_free) || (fcn_malloc && !fcn_free)) {
        return TOM_ERR_WRONG_INPUT;
    }

    /* Get the subregion and check the dimension. */
    if ((i=tom_io_calculate_sizes(voldims, subregion_, subregion, sampling_, sampling, binning_, binning, dims_sampled, dims)) != TOM_ERR_OK) {
        return i;
    }

    if (binning[0]>1 || binning[1]>1 || binning[2]>1) {
        switch (iotype) {
            case TOM_IO_TYPE_INT8:
            case TOM_IO_TYPE_INT16:
            case TOM_IO_TYPE_INT32:
                restype = TOM_IO_TYPE_DOUBLE;
                break;
            case TOM_IO_TYPE_FLOAT:
            case TOM_IO_TYPE_DOUBLE:
            case TOM_IO_TYPE_COMPLEX4:
                restype = iotype;
                break;
            default:
                return TOM_ERR_IOTYPE_NOT_SUPPORTED;
                break;
        }
    } else {
        restype = iotype;
    }
    /*restype = TOM_IO_TYPE_DOUBLE;*/


    if (!fcn_malloc || !fcn_free) {
        fcn_malloc = &malloc;
        fcn_free = &free;
    }

    pdata = fcn_malloc((size_t)dims[0] * (size_t)dims[1] * (size_t)dims[2] * tom_io_iotype_datasize(restype));

    if ((i=tom_io_read_vol(f, voldims, iotype, swaped, subregion, sampling, binning, pdata, restype)) != TOM_ERR_OK) {
        fcn_free(pdata);
    } else {
        *data = pdata;
        *restype_ = restype;
        dims_[0] = dims[0];
        dims_[1] = dims[1];
        dims_[2] = dims[2];
    }

    return i;
}



/***********************************************************************//**
 * \brief Read em-data from a file.
 *
 * Reads a volume (or 2D image) from an em-file.
 * \param[in]  filename   The name of the em-file to be read.
 * \param[in]  subregion_ A pointer to a 6-vector specifying the region
 *   of interest. If this parameter is set to \c NULL, the entire volume is
 *   read. Otherwise at positions [0,1,2] the coordinate of the first corner
 *   of the volume is given, at positions [3,4.5] the dimension of the region
 *   of interrest. The ROI is always relative to the original volume dimension
 *   no matter what binning or sampling is selected. An 2D image has heigth 1
 *   and thus subregion_[2] == 0 and subregion_[5] == 1. The subregion must
 *   be entirely inside the dimension of the volume.
 * \param[in]  sampling   Sample the volume down for each direction respectively.
 *   Passing NULL is the same as no sub-sampling. Otherwise it is an array of
 *   3 integers. A sampling factor of 5 for example, means that beginning
 *   with the first corner specified in subregion, each 5th element is taken
 *   into the output. If the volume dimension is not a multiple of the
 *   sampling factor, the last trailing element is taken into the resulting
 *   volume. For example sampling with distance 5 of a dimension of length
 *   22, results in a subsampled volume with width 5 (= ceil(22.0/5.0)). It
 *   follows that at least one element is always in the resulting volume,.
 *   A sampling factor of 0 is the same as 1 (i.e. no subsampling).
 * \param[in]  binning    specifies the binning factor along each direction
 *   x, y, z. Binning combines along each dimension \c i <tt>binning[i]</tt> elemtents
 *   by taking its mean. In contrast to sampling, the trailing elements are
 *   not taken into account. That means, binning with factor 5 along a dimension
 *   containing 22 elements, results in 4 values (= floor(22.0/5.0)). Through
 *   binning the resulting volume can become empty. In that case the function
 *   returns with an error. The parameters subregion_, sampling and binning
 *   can be combined. The subregion is sampled first and then binned. Again
 *   the first sampling step keeps the trailing elements of each dimension.
 *   Binning the sampled values down, throws the remainder of the number of
 *   samples away.
 * \param[out] header     Memory for taking the em-header. It can be set to
 *   \c NULL if it is not needed.
 * \param[out] data       A pointer where to save the allocated memory of the
 *   result. It must be freed falling \c free. If the function returns with an
 *   error, the memory is not allocated.
 * \param[out] dims       A pointer to 3 integers for saving the resulting
 *   dimension of the volume. It can not be ommitted.
 * \param[out] restype    A pointer to retrieve the data type of the volume.
 *   Usually it is the same as saved in the em-file. On case of binning it is
 *   converted to \c TOM_IO_TYPE_DOUBLE to avoid loosing information.
 * \return The function returns error-codes as defined in the tom_defines.h
 *   If it is different from \c TOM_ERR_OK, the output values are undefined.
 *
 * This wrapper combines the functions tom_io_em_read_header and
 * tom_io_read_vol_sac.
 **************************************************************************/
int tom_io_em_read(const char *filename, const uint32_t *subregion_, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void **data, uint32_t *dims, int *restype,
                 void *(*fcn_malloc)(size_t size),
                 void  (*fcn_free)(void *ptr)) {

    tom_io_em_header header_local;
    int i;
    FILE *f;

    if (!filename || !data || !dims || !restype) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (!header) {
        header = &header_local;
    }


    if ((i = tom_io_em_read_header(filename, "rb", header, &f)) != TOM_ERR_OK) {
        return i;
    }


    i = tom_io_read_vol_sac(f, header->dims, tom_io_em_get_iotype(header), tom_io_em_is_swaped(header),
                 subregion_, sampling, binning, data, dims, restype, fcn_malloc, fcn_free);

    fclose(f);

    return i;
}




/***********************************************************************//**
 * \brief Writes a volume to em-file.
 *
 * \param[in] f A file handle to an open file (write modus) seeked to the
 *   current position where to write.
 * \param[in] iotype_write Datatype to write to em-file. "Some" conversions
 *   are implemented.
 * \param[in] data Pointer to the volume to write.
 * \param[iotype] Numerical constant to describe the datatype of data. It is one
 *   of the TOM_IO_TYPE_XXX constants in tom_defines.h.
 * \param[in] sizex Size of the volume along x.
 * \param[in] sizey Size of the volume along y.
 * \param[in] sizez Size of the volume along z.
 * \param[in] stridex Stride in bytes along x. Set to 0 for default.
 * \param[in] stridey Stride in bytes along y. Set to 0 for default.
 * \param[in] stridez Stride in bytes along z. Set to 0 for default.
 * \param[in] swap_bytes Needs the bytes to be swapped for little endian vs. big endian.?
 * \param[in] fstridex Stride in the file in bytes along x. Set to 0 for default.
 * \param[in] fstridey Stride in the file in bytes along y. Set to 0 for default.
 * \param[in] fstridez Stride in the file in bytes along z. Set to 0 for default.
 * \returns TOM_ERR_OK in case of success. Otherwise an error status.
 *
 * The file is never closed.\n
 * By specifying the fstride parameters you can paste the volume into an
 * already existing file. In that case, the file must be large enough to hold the volume
 * (so that seeking does not fail). Otherwise the file will be enlarged if necessary.\n
 * When pasting the volume, it is up to the caller to ensure, where the volume is pasted
 * and that the preexisting type is the same as the new written one (\a iotype_write).
 **************************************************************************/
int tom_io_write_vol(       FILE *f, int iotype_write, const void *data, int iotype,
                            size_t sizex, size_t sizey, size_t sizez,
                            size_t stridex, size_t stridey, size_t stridez, int swap_bytes,
                            size_t fstridex, size_t fstridey, size_t fstridez) {

    int fseek_success;
    size_t numel;
    const int element_size = tom_io_iotype_datasize(iotype);
    const int element_size_write = tom_io_iotype_datasize(iotype_write);



    /* Check input parameters... */
    if (!f || ferror(f) || !data) {
        return TOM_ERR_WRONG_INPUT;
    }


    if (!tom_io_check_stride(element_size, sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        return TOM_ERR_WRONG_INPUT;
    }
    if (!tom_io_check_stride(element_size_write, sizex, sizey, sizez, &fstridex, &fstridey, &fstridez)) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (element_size<1 || element_size_write<1) {
        return TOM_ERR_IOTYPE_NOT_SUPPORTED;
    }

    fstridez -= sizey*fstridey;
    fstridey -= sizex*fstridex;
    fstridex -= element_size_write;

    numel = sizex * sizey * sizez;

    /* Set swap_bytes to false, if the data type is only one byte long. */
    swap_bytes = swap_bytes && iotype_write!=TOM_IO_TYPE_INT8;

    if (iotype_write==iotype && stridex==(size_t)element_size && stridey==stridex*sizex && stridez==sizey*stridey && !swap_bytes &&
        fstridex==0 && fstridey==0 && fstridez==0) {
        /* Write all the data at once. */
        if (fwrite(data, element_size, numel, f) != numel) {
            return TOM_ERR_WRITE_FILE;
        }
    } else {
        /* Write the data line per line. With either differing step width or swapping data. */
        uint32_t y, z;
        size_t i;

        void (*swap_fcn)(void *, size_t) = NULL;
        void (*copy_vector_fcn)(size_t, size_t, const char*, char*) = NULL;

        size_t copy_vector_var_n;
        size_t copy_vector_var_stridex;
        const char *copy_vector_var_src;
        char *copy_vector_var_dst;

        if (iotype!=iotype_write || swap_bytes || stridex!=(size_t)element_size) {
            switch (iotype_write) {
                case TOM_IO_TYPE_INT8:
                    switch (iotype) {
                        case TOM_IO_TYPE_INT8:     copy_vector_fcn = &copy_vector_int8_t__int8_t;           break;
                        default: return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                    }
                    break;
                case TOM_IO_TYPE_INT16:
                    if (swap_bytes) { swap_fcn = &tom_io_nswap02; }
                    switch (iotype) {
                        case TOM_IO_TYPE_INT8:      copy_vector_fcn = &copy_vector_int8_t__int16_t;         break;
                        case TOM_IO_TYPE_INT16:     copy_vector_fcn = &copy_vector_int16_t__int16_t;        break;
                        default: return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                    }
                    break;
                case TOM_IO_TYPE_INT32:
                    if (swap_bytes) { swap_fcn = &tom_io_nswap04; }
                    switch (iotype) {
                        case TOM_IO_TYPE_INT8:      copy_vector_fcn = &copy_vector_int8_t__int32_t;         break;
                        case TOM_IO_TYPE_INT16:     copy_vector_fcn = &copy_vector_int16_t__int32_t;        break;
                        case TOM_IO_TYPE_INT32:     copy_vector_fcn = &copy_vector_int32_t__int32_t;        break;
                        default: return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                    }
                    break;
                case TOM_IO_TYPE_FLOAT:
                    if (swap_bytes) { swap_fcn = &tom_io_nswap04; }
                    switch (iotype) {
                        case TOM_IO_TYPE_INT8:      copy_vector_fcn = &copy_vector_int8_t__float;           break;
                        case TOM_IO_TYPE_INT16:     copy_vector_fcn = &copy_vector_int16_t__float;          break;
                        case TOM_IO_TYPE_INT32:     copy_vector_fcn = &copy_vector_int32_t__float;          break;
                        case TOM_IO_TYPE_FLOAT:     copy_vector_fcn = &copy_vector_float__float;            break;
                        case TOM_IO_TYPE_DOUBLE:    copy_vector_fcn = &copy_vector_double__float;           break;
                        default: return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                    }
                    break;
                case TOM_IO_TYPE_DOUBLE:
                    if (swap_bytes) { swap_fcn = &tom_io_nswap08; }
                    switch (iotype) {
                        case TOM_IO_TYPE_INT8:      copy_vector_fcn = &copy_vector_int8_t__double;          break;
                        case TOM_IO_TYPE_INT16:     copy_vector_fcn = &copy_vector_int16_t__double;         break;
                        case TOM_IO_TYPE_INT32:     copy_vector_fcn = &copy_vector_int32_t__double;         break;
                        case TOM_IO_TYPE_FLOAT:     copy_vector_fcn = &copy_vector_float__double;           break;
                        case TOM_IO_TYPE_DOUBLE:    copy_vector_fcn = &copy_vector_double__double;          break;
                        default: return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                    }
                    break;
                case TOM_IO_TYPE_COMPLEX4:
                    if (swap_bytes) { swap_fcn = &tom_io_nswap_em_complex; }
                    switch (iotype) {
                        case TOM_IO_TYPE_INT8:      copy_vector_fcn = &copy_vector_int8_t__complex;         break;
                        case TOM_IO_TYPE_INT16:     copy_vector_fcn = &copy_vector_int16_t__complex;        break;
                        case TOM_IO_TYPE_INT32:     copy_vector_fcn = &copy_vector_int32_t__complex;        break;
                        case TOM_IO_TYPE_FLOAT:     copy_vector_fcn = copy_vector_float__complex;           break;
                        case TOM_IO_TYPE_COMPLEX4:  copy_vector_fcn = &copy_vector_complex__complex;        break;
                        case TOM_IO_TYPE_DOUBLE:    copy_vector_fcn = copy_vector_double__complex;          break;
                        default: return TOM_ERR_WRONG_IOTYPE_CONVERSION;
                    }
                    break;
                default:    /* Not expected*/
                    assert(1);
            }
            if (stridex==(size_t)element_size && iotype==iotype_write) {
                copy_vector_fcn = &copy_vector_memcpy; /* If there is only swapping needed, reset the copyfunction to memcpy. */
            }
            if (!(copy_vector_var_dst = (char *)malloc(sizex*element_size_write))) {
                return TOM_ERR_MALLOC;
            }
        } else {
            copy_vector_var_dst = (char *)data;
        }

        /* Set variable for the memcopy function. */
        copy_vector_var_n = sizex;
        if (copy_vector_fcn == &copy_vector_memcpy) {
            copy_vector_var_n *= element_size;
        }
        copy_vector_var_stridex = stridex;
        copy_vector_var_src = (char *)data;

        stridez -= sizey*stridey;

        for (z=0; z<sizez; z++) {
            for (y=0; y<sizey; y++) {
                if (copy_vector_fcn) {
                    /* Remind, that in case of copy_vector_fcn==&copy_vector_memcpy,
                       copy_vector_var_n is the number of bytes to copy.
                       Otherwise it is the number of elements per line and copy_vector_fcn
                       (using copy_vector_var_stridex) copies element out of the
                       data pointer into the contiguous buffer... */
                    copy_vector_fcn(copy_vector_var_n, copy_vector_var_stridex, copy_vector_var_src, copy_vector_var_dst);
                    if (swap_bytes) {
                        swap_fcn(copy_vector_var_dst, sizex);
                    }
                }
                if (fstridex) {
                    /* Write the first voxel of the line. */
                    if (fwrite(&copy_vector_var_dst[0*element_size_write], element_size_write, 1, f) != 1) {
                        if (copy_vector_fcn) { free(copy_vector_var_dst); }
                        return TOM_ERR_WRITE_FILE;
                    }
                    for (i=1; i<sizex; i++) {
                        ___fseek_cur(f, fstridex, fseek_success);
                        if (!fseek_success) {
                            if (copy_vector_fcn) { free(copy_vector_var_dst); }
                            return TOM_ERR_WRITE_FILE;
                        }
                        if (fwrite(&copy_vector_var_dst[i*element_size_write], element_size_write, 1, f) != 1) {
                            if (copy_vector_fcn) { free(copy_vector_var_dst); }
                            return TOM_ERR_WRITE_FILE;
                        }
                    }
                    if (y+1<sizey || z+1<sizez) {
                        ___fseek_cur(f, fstridex, fseek_success);
                        if (!fseek_success) {
                            if (copy_vector_fcn) { free(copy_vector_var_dst); }
                            return TOM_ERR_WRITE_FILE;
                        }
                    }
                } else {
                    if (fwrite(copy_vector_var_dst, element_size_write, sizex, f) != sizex) {
                        if (copy_vector_fcn) { free(copy_vector_var_dst); }
                        return TOM_ERR_WRITE_FILE;
                    }
                }
                if (copy_vector_fcn) {
                    copy_vector_var_src += stridey;
                } else {
                    copy_vector_var_dst += stridey;
                }
                if (fstridey && (y+1<sizey || z+1<sizez)) {
                    ___fseek_cur(f, fstridey, fseek_success);
                    if (!fseek_success) {
                        if (copy_vector_fcn) { free(copy_vector_var_dst); }
                        return TOM_ERR_WRITE_FILE;
                    }
                }
            }
            if (copy_vector_fcn) {
                copy_vector_var_src += stridez;
            } else {
                copy_vector_var_dst += stridez;
            }
            if (fstridez && z+1<sizez) {
                ___fseek_cur(f, fstridez, fseek_success);
                if (!fseek_success) {
                    if (copy_vector_fcn) { free(copy_vector_var_dst); }
                    return TOM_ERR_WRITE_FILE;
                }
            }
        }

        if (copy_vector_fcn) {
            free(copy_vector_var_dst);
        }
    }

    return TOM_ERR_OK;
}





/***********************************************************************//**
 * \brief Save a volume as em-file.
 *
 * Saves a volume to an em-file. Caution: existing files are overwritten.
 *
 * \param[in] filename The name of the file (existing files are overwritten)
 * \param[in] header   The header of the em-data. Must be initialized
 *   correctly. If the machine-code of the header says so, the data-elements
 *   are swaped before saving to file. Otherwise the header is saved directly.
 * \param[in] data     A pointer to the volume. The data-type and the dimension
 *   is saved in the header.
 * \param[in] stridex  How far two consecutive elements along x are distance from each
 *   other in bytes (0 defaults to a contiguous volume, thus the size of one data element
 *   where its type is specified by the header).
 * \param[in] stridey How far two consecutive elements along y are distance from each
 *   other in bytes (0 defaults to a contiguous volume, thus header->dims[0]*stridex).
 * \param[in] stridez How far two consecutive elements along z are distance from each
 *   other in bytes (0 defaults to a contiguous volume, thus header->dims[1]*stridey).
 * \return TOM_ERR_OK if no error occured. Otherwise an error-status is
 *   returned. In that case it is possible that at least a part of the file were written.
 *
 * If the entire volume should be saved, the stride parameters can be set to 0.
 * Otherwise they can be used if there is a gap between two consecutive data-elements,
 * lines or to consecutive "levels" (stridez). With them a subregion of the entire volume can be
 * saved to file. For that initialize the dimension of the header with the size of the subregion.
 * Pass as data not the pointer to the first element of the entire volume, but the pointer to
 * the first element of the subregion. The stride parameters can not be negative.
 **************************************************************************/
int tom_io_em_write(const char *filename, const tom_io_em_header *header, const void *data, int iotype, size_t stridex, size_t stridey, size_t stridez) {

    int swap_bytes, iotype_write, res;
    FILE *f;

    /* Check input parameters... */
    if (!filename || !header || !data) {
        return TOM_ERR_WRONG_INPUT;
    }
    if (!tom_io_em_is_valid_header(header, sizeof(*header))) {
        return TOM_ERR_WRONG_HEADER;
    }

    if (header->dims[0]<1 || header->dims[1]<1 || header->dims[2]<1) {
        return TOM_ERR_WRONG_HEADER;
    }
    iotype_write = tom_io_em_get_iotype(header);
    switch (iotype_write) {
        case TOM_IO_TYPE_INT8:
        case TOM_IO_TYPE_INT16:
        case TOM_IO_TYPE_INT32:
        case TOM_IO_TYPE_FLOAT:
        case TOM_IO_TYPE_DOUBLE:
        case TOM_IO_TYPE_COMPLEX4:
        break;
        default:
            return TOM_ERR_IOTYPE_NOT_SUPPORTED;
    }
    switch (iotype) {
        case TOM_IO_TYPE_INT8:
        case TOM_IO_TYPE_INT16:
        case TOM_IO_TYPE_INT32:
        case TOM_IO_TYPE_FLOAT:
        case TOM_IO_TYPE_DOUBLE:
        case TOM_IO_TYPE_COMPLEX4:
        break;
        default:
            return TOM_ERR_IOTYPE_NOT_SUPPORTED;
    }

    if (!tom_io_check_stride(tom_io_iotype_datasize(iotype), header->dims[0], header->dims[1], header->dims[2], &stridex, &stridey, &stridez)) {
        return TOM_ERR_WRONG_INPUT;
    }

    swap_bytes = tom_io_em_is_swaped(header);

    /* Open the file. */
    if (!(f = fopen(filename, "wb"))) {
        return TOM_ERR_OPEN_FILE;
    }

    {

        /* Write the tom_io_em_header to file. */
        tom_io_em_header header_write;
        const tom_io_em_header *pheader_write = header;

        assert(sizeof(tom_io_em_header) == 512);
        /* Does the compiler padd the struct wrong?
        * http://c-faq.com/struct/endpad.html
        * http://c-faq.com/struct/align.esr.html
        * http://c-faq.com/struct/padding.html
        * This should not happen since the largest type of tom_io_em_header is int32 with
        * 4 bytes, and all the fields of tom_io_em_header are nultiples of 4.
        * Otherwise the reading/writing of the header should be changed to
        * http://c-faq.com/stdio/extconform.html
        */

        if (swap_bytes) {
            header_write = *header;
            pheader_write = &header_write;
            tom_io_em_swap_header(&header_write);
        }
        if (fwrite(pheader_write, sizeof(*pheader_write), 1, f) != 1) {
            fclose(f);
            return TOM_ERR_WRITE_FILE;
        }
    }

    res = tom_io_write_vol(f, iotype_write, data, iotype, header->dims[0], header->dims[1], header->dims[2], stridex, stridey, stridez, swap_bytes, 0,0,0);

    fclose(f);
    return res;

}









/***********************************************************************//**
 * \brief Writes a particle stack to em-file by appending it to an existing one.
 *
 * Takes a volume and writes it to an existing em-file, by appending it
 * to the file along its last dimension (Z).
 * \param[in] filename The name of the existing em-file to which the
 *   particle should be appended.
 * \param[in] data Pointer to the volume.
 * \param[in] iotype The type of the data. It must be one of the numerical
 *   constants from tom_defines.h, named \c TOM_IO_TYPE_XXX. That is,
 *   the data-type of the existing em-file and the type of the data where
 *   the pointer points to, must be of the same type iotype.
 * \param[in] sizex The size along X of the volume of the data as well as
 *   of the volume in the em-file.
 * \param[in] sizey The size along Y of the volume of the data as well as
 *   of the volume in the em-file.
 * \param[in] sizez The size along Z of the volume. This can differ from the
 *   size in the already saved file.
 * \param[in] stridex The distance between two voxels along X in data.
 *   Measured in bytes (setting to 0 defaults to sizeof("iotype")).
 * \param[in] stridey The distance between two voxels along Y in data.
 *   Measured in bytes (setting to 0 defaults to stridex*sizex).
 * \param[in] stridey The distance between two voxels along Z in data.
 *   Measured in bytes (setting to 0 defaults to stridey*sizey).
 * \param[out] header_out Pointer to the em-header from the file AFTER
 *   a successfull write, i.e. if the function returns TOM_ERR_OK.
 *   If an error happened after successfully reading the file
 *   it contains the header of the existing em-file. In that
 *   case *header_read is set to 1. This happens for example if you try to
 *   append to an em-file with the wrong datatype or wrong size. In that case
 *   the em-file stays unchanged, the functions returns with an error,
 *   header_out contains the header of the em-file and *header_read is set
 *   to 1. If *header_read is 0, the content of header_out is undefined.\n
 *   Setting to \c NULL is ok, and means that you are not interested
 *   in the header.
 * \param[out] header_read Specifies whether header_out contains a valid header.
 *   if *header_read is 0, the content of header_out is undefined. If
 *   *header_read is 1, header_out either contains the header after writing
 *   successfully the file (if the function returns TOM_ERR_OK) or it contains
 *   the header of the file after opening it. This allows you in case of error,
 *   to trace down the reason. For example, was the iotype different, etc.?
 *   Can be set to \c NULL.
 * \param[in] allow_conversion If true, type conversions implemented in
 *   tom_io_write_vol are allowed. Otherwise if the datatype in the em-file
 *   differs from iotype, the function returns with error
 *   TOM_ERR_WRONG_IOTYPE_CONVERSION
 * \returns TOM_ERR_OK in case of success. Otherwise an error-code as defined
 *   in tom_defines.h
 *
 * In case of successfully writing, the only parameter changed in the em-header
 * on file is the size the last dimension Z.\n
 * If the file does not exist, it is not created. This is implemented that
 * way because creating a new em-file also needs an entire em-header, and it
 * is unclear how to fill it up. Thus it is suggested, in case of error, to
 * use tom_io_em_write to write a new file.\n
 * If an IO-error happens after beginning to append the volume, thats
 * a serious problem because the em-file may be damaged. However the IO-error
 * already indicates a problem...\n
 * An other serious problem is, if a segemention fault happens (because one
 * of the other input-pointers is invalid). In that case the program crashes,
 * leaving a damaged file. Its up to the caller to avoid this.
 **************************************************************************/
int tom_io_em_write_append_stack(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_em_header *header_out, int *header_read, int allow_conversion) {

    const int element_size = tom_io_iotype_datasize(iotype);
    tom_io_em_header header;
    uint32_t sizez_old;
    FILE *f;
    int res;
    int swap_bytes;
    int iotype_write;


    if (header_read) { *header_read = 0; }


    if (!filename || !data) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (!tom_io_em_set_iotype(&header, iotype)) {
        return TOM_ERR_IOTYPE_NOT_SUPPORTED;
    }

    if (!tom_io_check_stride(element_size, sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (!header_out) {
        header_out = &header;
    }

    if ((res = tom_io_em_read_header(filename, "r+b", header_out, &f)) != TOM_ERR_OK) {
        return res;
    }
    fseek(f, 0, SEEK_END);

    if (ferror(f)) {
        return TOM_ERR_READ_FILE;
    }
    if (header_read) { *header_read = 1; }

    iotype_write = tom_io_em_get_iotype(header_out);

    if (header_out->dims[0]!=sizex || header_out->dims[1]!=sizey) {
        return TOM_ERR_WRONG_DATA_SIZE;
    }
    sizez_old = header_out->dims[2];
    if ((size_t)sizez + (size_t)sizez_old > 0xFFFFFFFF) {
        return TOM_ERR_VOLUME_TOO_LARGE;
    }

    if (!allow_conversion && iotype_write!=iotype) {
        return TOM_ERR_WRONG_IOTYPE_CONVERSION;
    }

    swap_bytes = tom_io_em_is_swaped(header_out);

    res = tom_io_write_vol(f, iotype_write, data, iotype, sizex, sizey, sizez, stridex, stridey, stridez, swap_bytes, 0,0,0);

    if (res == TOM_ERR_OK && !ferror(f)) {
        fseek(f, 4*1 + 2*4, SEEK_SET); /* Seek to the position of the dimension Z: MAGIC(4) + DIMX(4) + DIMY(4) */

        sizez_old += sizez;
        if (swap_bytes) {
            tom_io_swap04(&sizez_old);
        }

        if (fwrite(&sizez_old, sizeof(sizez_old), 1, f) != 1) {
            res = TOM_ERR_WRITE_FILE;
        }
    }

    fclose(f);
    if (res == TOM_ERR_OK) {
        header_out->dims[2] += sizez;
    }
    return res;

}






/***********************************************************************//**
 * \param[in] filename The filename.
 * \param[in] data Pointer to the data.
 * \param[in] iotype
 * \param[in] sizex sizex - z of the volume to be written.
 * \param[in] stridex stridex - z The strides of the volume.
 * \param[in] header_out out If the header of the file on disk has been read, it is stored here.
 * \param[in] header_read True if header_out has been filled.
 * \param[in] allow_conversion
 * \param[in] first_voxel The first voxel which is overwritten by the volume.
 * \param[in] sampling
 **************************************************************************/
int tom_io_em_write_paste(  const char *filename,
                            const void *data,
                            int iotype,
                            uint32_t sizex, uint32_t sizey, uint32_t sizez,
                            size_t stridex, size_t stridey, size_t stridez,
                            tom_io_em_header *header_out, int *header_read,
                            int allow_conversion,
                            const uint32_t *first_voxel, const uint32_t *sampling) {

    const int element_size = tom_io_iotype_datasize(iotype);
    tom_io_em_header header;
    int iotype_write;
    int swap_bytes;
    FILE *f;
    int res;
    size_t i;
    uint32_t first_voxel_[3] = { 0,0,0 };
    uint32_t sampling_[3] = { 1,1,1 };
    size_t fstride[3];

    assert(sizeof(header) == 512);

    if (header_read) { *header_read = 0; }

    if (!filename || !data) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (element_size < 1) { return TOM_ERR_IOTYPE_NOT_SUPPORTED; }

    if (!tom_io_check_stride(element_size, sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        return TOM_ERR_WRONG_INPUT;
    }

    if (!header_out) { header_out = &header; }
    if (first_voxel) { for (i=0; i<3; i++) { first_voxel_[i] = first_voxel[i]                  ; } }
    if (sampling   ) { for (i=0; i<3; i++) { sampling_   [i] = sampling   [i] ? sampling[i] : 1; } }

    if ((res = tom_io_em_read_header(filename, "r+b", header_out, &f)) != TOM_ERR_OK) {
        return res;
    }

    if (ferror(f)) {
        return TOM_ERR_READ_FILE;
    }

    if (header_read) { *header_read = 1; }

    if (first_voxel_[0] + (sizex-1)*sampling_[0] + 1 > header_out->dims[0] ||
        first_voxel_[1] + (sizey-1)*sampling_[1] + 1 > header_out->dims[1] ||
        first_voxel_[2] + (sizez-1)*sampling_[2] + 1 > header_out->dims[2]) {
        return TOM_ERR_VOLUME_TOO_LARGE;
    }

    swap_bytes = tom_io_em_is_swaped(header_out);
    iotype_write = tom_io_em_get_iotype(header_out);

    assert(swap_bytes==1 || swap_bytes==0 && tom_io_iotype_datasize(iotype_write)>0);

    if (!allow_conversion && iotype_write!=iotype) {
        return TOM_ERR_WRONG_IOTYPE_CONVERSION;
    }

    fstride[0] = sampling_[0] * tom_io_iotype_datasize(iotype_write);
    fstride[1] = sampling_[1] * header_out->dims[0] * tom_io_iotype_datasize(iotype_write);
    fstride[2] = sampling_[2] * header_out->dims[1] * header_out->dims[0] * tom_io_iotype_datasize(iotype_write);

    {
        int fseek_success;
        ___fseek_set( f, sizeof(header) + ((first_voxel_[2]*header_out->dims[1] + first_voxel_[1])*header_out->dims[0] + first_voxel_[0] )*tom_io_iotype_datasize(iotype_write), fseek_success);
        if (!fseek_success) {
            return TOM_ERR_READ_FILE;
        }
    }


    res = tom_io_write_vol(f, iotype_write, data, iotype, sizex, sizey, sizez, stridex, stridey, stridez, swap_bytes, fstride[0], fstride[1], fstride[2]);

    fclose(f);
    return res;
}








