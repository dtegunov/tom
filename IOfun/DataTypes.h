#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

typedef unsigned int       	DWORD;
typedef int                 	BOOL;
typedef unsigned char       	BYTE;
typedef unsigned short      	WORD;
typedef float 			FLOAT;
typedef unsigned int		UINT;
typedef int			LONG;
typedef unsigned int		ULONG;
typedef unsigned short 		USHORT;

//typedef unsigned char			TCHAR:

typedef enum IMAGE_TYPE
{
	DT_UCHAR			= 1,
	DT_USHORT			= 2,
	DT_SHORT			= 3,
	DT_LONG				= 4,
	DT_FLOAT			= 5,
	DT_DOUBLE			= 6,
	DT_COMPLEX			= 7,
	DT_STRING			= 8,
	DT_BINARY			= 9,
	DT_RGB8				= 10,
	DT_RGB16			= 11,
	DT_EMVECTOR			= 12
};
