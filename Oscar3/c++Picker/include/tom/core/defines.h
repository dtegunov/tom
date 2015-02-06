/**************************************************************************//**
 * \file defines.h
 * Defines for several routines.
 ******************************************************************************/
#ifndef ___INCLUDE_CORE__DEFINES_HPP__
#define ___INCLUDE_CORE__DEFINES_HPP__



/* Error/returning codes by (some) functions. Make
 * these values nagtive!!!! */
#define TOM_ERR_OK                        ((int)( 0))
#define TOM_ERR_WRONG_INPUT               ((int)(-1))
#define TOM_ERR_WRONG_HEADER              ((int)(-2))
#define TOM_ERR_SUBREGION_OUT             ((int)(-3))
#define TOM_ERR_NOT_YET_IMPLEMENTED       ((int)(-4))
#define TOM_ERR_MALLOC                    ((int)(-5))
#define TOM_ERR_OPEN_FILE                 ((int)(-6))
#define TOM_ERR_FILESIZE_MISMATCH         ((int)(-7))
#define TOM_ERR_READ_FILE                 ((int)(-8))
#define TOM_ERR_BINNING_TOO_HIGH          ((int)(-9))
#define TOM_ERR_IOTYPE_NOT_SUPPORTED      ((int)(-10))
#define TOM_ERR_WRONG_IOTYPE_CONVERSION   ((int)(-11))
#define TOM_ERR_WRITE_FILE                ((int)(-12))
#define TOM_ERR_NO_COMPLEX_BINNING        ((int)(-13))
#define TOM_ERR_WRONG_DATA_SIZE           ((int)(-14))
#define TOM_ERR_VOLUME_TOO_LARGE          ((int)(-15))




/* These are integer numbers defining the raw datatype.
 * The numerical value corresponds to the parameter in the em-header.
 * This is arbitrary. Don't rely on the numerical values of the defines.
 * Only make them >= 1. See also iotype_datasize */
#define TOM_IO_TYPE_UNKNOWN                    ((int)0)
#define TOM_IO_TYPE_INT8                       ((int)1)
#define TOM_IO_TYPE_INT16                      ((int)2)
#define TOM_IO_TYPE_INT32                      ((int)4)
#define TOM_IO_TYPE_FLOAT                      ((int)5)
#define TOM_IO_TYPE_COMPLEX4                   ((int)8)
#define TOM_IO_TYPE_DOUBLE                     ((int)9)
#define TOM_IO_TYPE_INT64                     ((int)10)



#endif



