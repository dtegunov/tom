// FileTool.h: interface for the FileTool class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef _DATATYPES
	#include "DataTypes.h"
	#define _DATATYPES
#endif
#ifndef _EMVECT
	#include "EMVector.h"
	#define _EMVECT
#endif
#ifndef _TEMPLATE
	#include "var_list.cpp"
	#define _TEMPLATE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <iostream>


#define SWITCHENDIAN(d,s,n,b) { \
	if((b)) \
		for(int i=0,j=(n)-1;j>=0;i++,j--) \
			((BYTE *)(d))[i]=((BYTE *)(s))[j]; \
		else  \
			memcpy((d),(s),(n)); \
}


typedef int HRESULT;

class CFileTool  
{
public:
	CFileTool();
	virtual ~CFileTool();
	
	HRESULT DoReadTiff(TCHAR *tcpFilename, var_list<BYTE> *pDst, var_list<EMVECTOR> *pEmVector, var_list<BYTE> *pThumbNail);

private:
	
	void Swapper( char *cp, long lLen );
	void ConvertV1ToV2( void *p,void *TIFF_Data );
	void ApplyUnicodeString( char *cpDst, char *cpSrc );
	void GetDateInformation( long lSince1_1_1970, long *lpDay, long *lpMonth, long *lpYear );
	
	void ReadFromFile(FILE*, BYTE*, size_t, BOOL);
	void ReadFromFile(const BYTE *pBuffer, ULONG &pos, BYTE *pcDst, size_t size, BOOL bEndian);
	bool BufferReadTiff(const BYTE *pBuffer, var_list<BYTE> *pDst, var_list<EMVECTOR> *pEmVector, var_list<BYTE> *pThumbNail);

	struct TIFF_DATA {
		DWORD	ImageWidth,
				ImageLength,
				BitsPerSample,
				Compression,
				PhotoInterpret,
				StripOffsets,
				StripOffsetsCount,
				RowsPerStrip,
				StripByteCounts,
				XResolutionZ,
				XResolutionN,
				YResolutionZ,
				YResolutionN,
				ResolutionUnit;
	};

	struct TIFF_FIELD {
		WORD	tag;
		WORD	type;
		DWORD	count;
		DWORD	value;
	};

protected:
	
	void GetMinMaxMean( long lDataType, void *vpData, long lXDim, long lYDim, float *fpMin, float *fpMax, float *fpMean );
};

inline void CFileTool::ReadFromFile(const BYTE *pBuffer, ULONG &pos, BYTE *pcDst, size_t size, BOOL bEndian) 
{
	const BYTE *pcBuf = &pBuffer[pos];
	SWITCHENDIAN(pcDst, pcBuf, size, bEndian);
	pos+=size;
}

