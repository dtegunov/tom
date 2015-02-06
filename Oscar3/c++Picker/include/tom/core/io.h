/****************************************************************************//**
 * \file io.h
 * \brief The header file of the io-functions in tom_io.c
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.11.2007
 *
 * Contains functions to read and write data from/to a file.
 * The used format is EM
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__IO_H__
#define ___INCLUDE_CORE__IO_H__

#define _FILE_OFFSET_BITS 64
#ifndef __LARGEFILE_SOURCE
    #define __LARGEFILE_SOURCE
#endif
#define _FILE_OFFSET_BITS 64


#ifdef __cplusplus
    extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#if defined __GNUC__ || defined __xlC__
    #include <stdint.h>
#else
    #include "system/pstdint.h"
#endif



#include "tom/core/defines.h"





/****************************************************************************//**
 * \brief Structure for the em-header.
 *
 * A binary compatible structure to the header of the em-file.
 * It contains the exactly same bits as stored in the file, exept the case
 * when the integers are swaped.
 *******************************************************************************/
typedef struct tom_io_em_header {
    int8_t machine;             /**< Byte 1: Machine Coding
                                            (OS-9;      0),
                                            (VAX;       1),
                                            (Convex;    2),
                                            (SGI;       3),
                                            (Mac;       5),
                                            (PC;        6). */
    int8_t byte2;               /**< General purpose. On OS-9 system: 0 old version 1 is new version. */
    int8_t byte3;               /**< Not used in standard EM-format, if this byte is 1 the header is abandoned. */
    int8_t type;                /**< Data Type Coding. */
    uint32_t dims[3];           /**< Three long integers (3x4 bytes) are image size in x, y, z Dimension. */
    int8_t comment[80];         /**< 80 Characters as comment. */
    int32_t emdata[40];         /**< 40 long integers (4 x 40 bytes) are user defined parameters. */
    int8_t  userdata[256];      /**< 256 Byte with userdata, i.e. the username. */
} tom_io_em_header;



/****************************************************************************//**
 * \brief Complex element corresponding to the complex em-type.
 *
 * In this format complex data is saved in em-format.
 *******************************************************************************/
typedef struct {
    float re; /**< real part */
    float im; /**< imaginary part */
} tom_io_em_complex;



int tom_io_iotype_datasize(int type);


#define TOM_IO_MACHINE_BYTEORDER_BIG_ENDIAN      0
#define TOM_IO_MACHINE_BYTEORDER_LITTLE_ENDIAN   1
int tom_io_machine_byteorder();
void tom_io_swap(void *x, size_t size);

int tom_io_em_read_header(const char *fname, const char *mode, tom_io_em_header *header, FILE **fhandle);
int tom_io_em_is_valid_header(const void *header, size_t size);
int tom_io_em_is_swaped(const tom_io_em_header *header);
int tom_io_em_swap_header(tom_io_em_header *header);
int tom_io_em_get_iotype(const tom_io_em_header *header);
int tom_io_em_set_iotype(tom_io_em_header *header, int iotype);
void tom_io_em_get_microscope_name(const tom_io_em_header *header, char *name, size_t bytes_reserved);

int tom_io_em_read( const char *fname, const uint32_t *subregion_, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header,
                    void **data, uint32_t *dims, int *restype,
                    void *(*fcn_malloc)(size_t size),
                    void  (*fcn_free)(void *ptr));
int tom_io_read_vol_sac(FILE *f, const uint32_t voldims[], int iotype, int swaped,
                    const uint32_t *subregion_, const uint32_t *sampling_, const uint32_t *binning_,
                    void **data, uint32_t *dims_, int *restype_,
                    void *(*fcn_malloc)(size_t size),
                    void  (*fcn_free)(void *ptr));
int tom_io_read_vol(FILE *f, const uint32_t voldims[], int iotype, int swaped,
                    const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning,
                    void *pdata, int restype);
int tom_io_calculate_sizes(const uint32_t *voldims,
                    const uint32_t *subregion_, uint32_t *subregion_out,
                    const uint32_t *sampling_, uint32_t *sampling_out,
                    const uint32_t *binning_, uint32_t *binning_out,
                    uint32_t *dims_sampled_out,
                    uint32_t *dims_out);

/*int tom_io_write_vol(FILE *f, int iotype_write, const void *data, int iotype, size_t sizex, size_t sizey, size_t sizez, size_t stridex, size_t stridey, size_t stridez, int swap_bytes);*/
int tom_io_write_vol(       FILE *f, int iotype_write, const void *data, int iotype,
                            size_t sizex, size_t sizey, size_t sizez,
                            size_t stridex, size_t stridey, size_t stridez, int swap_bytes,
                            size_t fstridex, size_t fstridey, size_t fstridez);

int tom_io_em_write(const char *filename, const tom_io_em_header *header, const void *data, int iotype, size_t stridex, size_t stridey, size_t stridez);
int tom_io_em_write_append_stack(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_em_header *header_out, int *header_read, int allow_conversion);
int tom_io_em_write_paste(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_em_header *header_out, int *header_read, int allow_conversion, const uint32_t *first_voxel, const uint32_t *sampling);

int tom_io_check_stride(size_t size_of_type, size_t sizex, size_t sizey, size_t sizez, size_t *stridex, size_t *stridey, size_t *stridez);



#ifdef __cplusplus
    } /* extern "C" */
#endif


#endif
