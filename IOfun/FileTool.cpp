// FileTool.cpp: implementation of the FileTool class.
//
//////////////////////////////////////////////////////////////////////

#include "FileTool.h"
//#include "File.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

#define SIZEOF_EMVECTOR_IN_FILE	6312

/* TIFF FORMAT NUMBERS */
#define TYPE_BYTE		1
#define TYPE_ASCII		2
#define TYPE_SHORT		3
#define TYPE_LONG		4
#define TYPE_RATIONAL	5
#define TYPE_SBYTE		6
#define TYPE_UNDEFINED	7
#define TYPE_SSHORT		8
#define TYPE_SLONG		9
#define TYPE_SRATIONAL	10
#define TYPE_FLOAT		11
#define TYPE_DOUBLE		12

#define TAG_TEM_DATA	37706			// TVIPS tag 1
#define TAG_OVL_DATA	37707			// TVIPS tag 2
#define ROWS_PER_STRIP	1				// use always 1 row per TIFF strip


// In this macro I use this simple formula because 2000 is also a schaltjahr
#define IS_SCHALTJAHR(j)   ( (((j)/4)*4) == (j) )

#define RETURN_WITH_ERROR(txt, err) {printf(txt); return ( err ); }
#define RETURN_WITH_ERRORV(txt) {printf(txt); return; }


CFileTool::CFileTool()
{

}

CFileTool::~CFileTool()
{

}

/*
 *	CLU:
 *
 *	Read the data from a given TIFF file and fill - if the file contains the data - the
 *	EM-Vector
 */
void CFileTool::Swapper(char *cp, long lLen)
{
	char c;

	if( lLen == 2 )
	{
		c = cp[0];
		cp[0] = cp[1];
		cp[1] = c;
	}

	if( lLen == 4 )
	{
		c = cp[0];
		cp[0] = cp[3];
		cp[3] = c;
		c = cp[1];
		cp[1] = cp[2];
		cp[2] = c;
	}
}

void CFileTool::ConvertV1ToV2( void *p, void *p2 )
{
	//
	// This function evaluates the version 1 header of the TVIPS TIFF tag and converts it as good as possible
	// into the version 2 header.
	//

	TIFF_DATA	*pTIFF_Data;
	pTIFF_Data = (TIFF_DATA*)p2;

	EMVECTOR *pEMV;
	int i;

	pEMV = (EMVECTOR *)p;
	
	if( pTIFF_Data!= NULL )
		pEMV->ud_v2.lImgDataType = (pTIFF_Data->BitsPerSample==8?1:3);							// 8bit->DT_UCHAR, 16bit->DT_SHORT

	ApplyUnicodeString( (char*)&pEMV->ud_v2.szImgName[0],"");								// image name
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szImgFolder[0],"");								// path to image
	pEMV->ud_v2.lImgSizeX = (long)pEMV->ud_v1.lCCDPixelXY;									// images in EMMENU 3.0 are 2dim and square !!!
	pEMV->ud_v2.lImgSizeY = (long)pEMV->ud_v1.lCCDPixelXY;									// images in EMMENU 3.0 are 2dim and square !!!
	pEMV->ud_v2.lImgSizeZ = 0;
	pEMV->ud_v2.lImgSizeE = 0;
	long lDay, lMonth, lYear;
	GetDateInformation( pEMV->ud_v1.lDate, &lDay, &lMonth, &lYear );
	pEMV->ud_v2.lImgCreationDate = (65536*lYear + 256*lMonth + lDay);
	pEMV->ud_v2.lImgCreationTime = static_cast<long>(pEMV->ud_v1.lTime*.01);				// time ( from 1/100 s into s )
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szImgComment[0],&pEMV->ud_v1.strComment[0]);	//  image comment
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szImgHistory[0],"");							//  image history

	pEMV->ud_v2.lImgType = 0;																// 0=image
	pEMV->ud_v2.lImgDisplay = 0;															// 0=display as grey value image

	// The following formula contains some multiplication and division factors
	// by 1000 or 10000, respectively. They come from the different units, the 
	// individual structure elements use. For better readibility they remain in
	// the code:
	// ud_v1.lCCDSinglePixelSize * 1000 -> pixelsize in nm
	// ud_v1.lElOptMag*1000.f			-> real mag
	// (ud_v1.lPostMag*0.0001)			-> real postmag

	pEMV->ud_v2.fImgDistX = 
	pEMV->ud_v2.fImgDistY = ((pEMV->ud_v1.lCCDSinglePixelSize * 1000)*pEMV->ud_v1.lCCDBinningFactor)/(pEMV->ud_v1.lElOptMag*1000.f*(pEMV->ud_v1.lPostMag*0.0001));
	pEMV->ud_v2.fImgDistZ = 0.f;
	pEMV->ud_v2.fImgDistE = 0.f;
	
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szTemType[0],"");								// TEM name
	pEMV->ud_v2.fTemHT = pEMV->ud_v1.lHighTension*1000.f;									// HT
	pEMV->ud_v2.lTemMode = 0;																// bright field
	pEMV->ud_v2.fTemMagScr = pEMV->ud_v1.lElOptMag*1000.f;									// mag
	pEMV->ud_v2.fTemMagCor = 1.f;															// new factor in EM4, so set it to 1.0
	pEMV->ud_v2.fTemMagPst = pEMV->ud_v1.lPostMag*.0001f;									// postmag
	pEMV->ud_v2.lTemStgType = 0;															// 0 = manual
	pEMV->ud_v2.lTemShutter = 0;															// 0 = ???
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szCamType[0],"EMMENU 3.0 image");				// camera name
	pEMV->ud_v2.fCamPixel[0] = pEMV->ud_v1.lCCDSinglePixelSize * 1000;
	pEMV->ud_v2.fCamPixel[1] = pEMV->ud_v1.lCCDSinglePixelSize * 1000;
	pEMV->ud_v2.lCamOffX = static_cast<long>(pEMV->ud_v1.lCCDOffsetX);
	pEMV->ud_v2.lCamOffY = static_cast<long>(pEMV->ud_v1.lCCDOffsetY);
	pEMV->ud_v2.lCamBinX = static_cast<long>(pEMV->ud_v1.lCCDBinningFactor);
	pEMV->ud_v2.lCamBinY = static_cast<long>(pEMV->ud_v1.lCCDBinningFactor);
	pEMV->ud_v2.fCamExpTime = pEMV->ud_v1.lCCDExposureTime;
	pEMV->ud_v2.fCamGain = pEMV->ud_v1.lCCDGain;
	pEMV->ud_v2.fCamSpeed = pEMV->ud_v1.lCCDReadOutSpeed;
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szCamFlat[0],"");								// flat image link
	pEMV->ud_v2.fCamSense = pEMV->ud_v1.lCCDSensitivity;
	pEMV->ud_v2.fCamDose = 0.f;
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szAdaTietzSpecInfo[0],"");						// tecnai specific data 1
	ApplyUnicodeString( (char*)&pEMV->ud_v2.szAdaTietzMicInfo[0],"");						// tecnai specific data 2

	for( i=0; i<32; i++ )
	{
		pEMV->ud_v2.fImgMisc[i]			= 0;
		pEMV->ud_v2.fTemAberr[i]			= 0;
		pEMV->ud_v2.fTemEnergy[i]			= 0;
		pEMV->ud_v2.fTemMisc[i]			= 0;
		pEMV->ud_v2.fCamMisc[i]			= 0;

		if( i< 16 )
		{
			pEMV->ud_v2.fImgScaling[i]	= 0;
		}

		if( i < 7 )
		{
			pEMV->ud_v2.fTemTiling[i]		= 0;
		}

		if( i < 5 )
		{
			pEMV->ud_v2.fTemStgPos[i]		= 0;
		}

		if( i < 3 )
		{
			pEMV->ud_v2.fTemIllum[i]		= 0;
		}

		if( i < 2 )
		{
			pEMV->ud_v2.fTemImgShift[i]	= 0;
			pEMV->ud_v2.fTemBeamShift[i]	= 0;
			pEMV->ud_v2.fTemBeamTilt[i]	= 0;
		}
	}
}


void CFileTool::ApplyUnicodeString( char *cpDst, char *cpSrc )
{
	while( (*cpSrc) != 0 )
	{
		*cpDst = *cpSrc;
		cpDst++;
		*cpDst = 0;
		cpDst++;
		cpSrc++;
	}
	*cpDst = 0;
	cpDst++;
	*cpDst = 0;
	cpDst++;
}

void CFileTool::GetDateInformation( long lSince1_1_1970, long *lpDay, long *lpMonth, long *lpYear )
{
	long lDay, lMonth, lYear, lDaysLeft;
	lYear = 1970;
	lMonth = 1;
	lDay = 1;

	*lpDay = lDay;
	*lpMonth = lMonth;
	*lpYear = lYear;

	lDaysLeft = lSince1_1_1970+1;		// 1.1.1970 is day 1

	if( lDaysLeft >= 0 )
	{
		lYear = 1970;
		while( lDaysLeft > 0 )
		{
			if( IS_SCHALTJAHR(lYear-1) && lDaysLeft > 366 )
			{
				lYear++;
				lDaysLeft-= 366;
				continue;
			}
			if( lDaysLeft > 365 )
			{
				lYear++;
				lDaysLeft-= 365;
				continue;
			}

			*lpYear = lYear;
			break;
		}

		if( IS_SCHALTJAHR( lYear) )
		{
			if( lDaysLeft >= 336 ){	*lpMonth = 12;	*lpDay = lDaysLeft-335;		return; }
			if( lDaysLeft >= 306 ){	*lpMonth = 11;	*lpDay = lDaysLeft-305;		return; }
			if( lDaysLeft >= 275 ){	*lpMonth = 10;	*lpDay = lDaysLeft-274;		return; }
			if( lDaysLeft >= 245 ){	*lpMonth = 9;	*lpDay = lDaysLeft-244;		return; }
			if( lDaysLeft >= 214 ){	*lpMonth = 8;	*lpDay = lDaysLeft-213;		return; }
			if( lDaysLeft >= 183 ){	*lpMonth = 7;	*lpDay = lDaysLeft-182;		return; }
			if( lDaysLeft >= 153 ){	*lpMonth = 6;	*lpDay = lDaysLeft-152;		return; }
			if( lDaysLeft >= 122 ){	*lpMonth = 5;	*lpDay = lDaysLeft-121;		return; }
			if( lDaysLeft >=  92 ){	*lpMonth = 4;	*lpDay = lDaysLeft-91;		return; }
			if( lDaysLeft >=  61 ){	*lpMonth = 3;	*lpDay = lDaysLeft-60;		return; }
			if( lDaysLeft >=  32 ){	*lpMonth = 2;	*lpDay = lDaysLeft-32;		return; }
			if( lDaysLeft >=   1 ){	*lpMonth = 1;	*lpDay = lDaysLeft;			return; }
		}
		else
		{
			if( lDaysLeft >= 335 ){	*lpMonth = 12;	*lpDay = lDaysLeft-334;		return; }
			if( lDaysLeft >= 305 ){	*lpMonth = 11;	*lpDay = lDaysLeft-304;		return; }
			if( lDaysLeft >= 274 ){	*lpMonth = 10;	*lpDay = lDaysLeft-273;		return; }
			if( lDaysLeft >= 244 ){	*lpMonth = 9;	*lpDay = lDaysLeft-243;		return; }
			if( lDaysLeft >= 213 ){	*lpMonth = 8;	*lpDay = lDaysLeft-212;		return; }
			if( lDaysLeft >= 182 ){	*lpMonth = 7;	*lpDay = lDaysLeft-181;		return; }
			if( lDaysLeft >= 152 ){	*lpMonth = 6;	*lpDay = lDaysLeft-151;		return; }
			if( lDaysLeft >= 121 ){	*lpMonth = 5;	*lpDay = lDaysLeft-120;		return; }
			if( lDaysLeft >=  91 ){	*lpMonth = 4;	*lpDay = lDaysLeft-90;		return; }
			if( lDaysLeft >=  60 ){	*lpMonth = 3;	*lpDay = lDaysLeft-59;		return; }
			if( lDaysLeft >=  32 ){	*lpMonth = 2;	*lpDay = lDaysLeft-31;		return; }
			if( lDaysLeft >=   1 ){	*lpMonth = 1;	*lpDay = lDaysLeft;			return; }
		}
	}

	return;
}


/*
 * CLU: reads data (byte stream) from a file into a destination buffer and does little - big
 * endian conversion.
 *
 * Parameters:
 * FILE		*pFile:		pointer to FILE
 * BYTE		*pcDst:		pointer to destination buffer
 * size_t	size:		number of bytes to be copied
 * BOOL		bEndian:	TRUE: reverse byte order
 *						FALSE: just copy
 *
 * Note: the calling routine must ensure the destination buffer size is sufficient
 *
 */
void CFileTool::ReadFromFile(FILE *pFile, BYTE *pcDst, size_t size, BOOL bEndian) {

	BYTE*	pcBuf;																// CLU: read buffer

	if (!(pcBuf = (BYTE *)malloc(size))) {
		RETURN_WITH_ERRORV("Out of Memory!");
	}

	if (fread(pcBuf, sizeof(BYTE), size, pFile) != size) {						// CLU: corrupt file
		free(pcBuf);

		RETURN_WITH_ERRORV("Cannot read from file!");
	}

	SWITCHENDIAN(pcDst, pcBuf, size, bEndian);												// CLU: interpret data
		
	free(pcBuf);

	//return(S_OK);
}


// K:B 
bool CFileTool::BufferReadTiff(const BYTE *pBuffer, var_list<BYTE> *pDst, var_list<EMVECTOR> *pEmVector, var_list<BYTE> *pThumbNail) 
{
	UINT		nStripOffsetsSize,
			nStripBytesSize;
	LONG		lFilePosition;
	DWORD		*pdwStripOffsets	= NULL,
			*pdwBytesPerStrip	= NULL;
	double		dXResolution 		= 72,
			dYResolution 		= 72;

//	HRESULT		hr;
	ULONG		pos = 0L;

	TIFF_DATA	TIFF_Data;
	TIFF_FIELD	TIFF_Field;

	
	if (!pDst && !pEmVector && !pThumbNail) {
		RETURN_WITH_ERROR("Illegal parameter (pDst and pEmVector and pThumbNail are NULL", 1);
	}

	try {

		/*
		 * CLU: initiate the TIFF data fields with -1
		 */
		TIFF_Data.ImageWidth		=
		TIFF_Data.ImageLength		=
		TIFF_Data.BitsPerSample		=
		TIFF_Data.Compression		=
		TIFF_Data.PhotoInterpret	=
		TIFF_Data.StripOffsets		=
		TIFF_Data.StripOffsetsCount	=
		TIFF_Data.RowsPerStrip		=
		TIFF_Data.StripByteCounts	=
		TIFF_Data.XResolutionZ		=
		TIFF_Data.XResolutionN		=
		TIFF_Data.YResolutionZ		=
		TIFF_Data.YResolutionN		=
		TIFF_Data.ResolutionUnit	=	-1;

		/*
		 * CLU: check byte order and TIFF signature
		 */
		BYTE	cBuf_4[4];															// CLU: read buffer

		ReadFromFile(pBuffer, pos, cBuf_4, 4, FALSE);

		if ((*(int *)cBuf_4 != 0x002a4949 && *(int *)cBuf_4 != 0x2a004d4d)) {		// CLU: no tiff file
			RETURN_WITH_ERROR("No TIFF file detected", 1);
		}

		BOOL bBigEndian = (*(short *)cBuf_4 == 0x4d4d);								// CLU: the byte order

		/*
		 * CLU: get first IFD
		 */
		ULONG IFD_Offset;

		ReadFromFile(pBuffer, pos, (BYTE *)&IFD_Offset, sizeof(ULONG), bBigEndian);

		if (IFD_Offset == 0) {														// CLU: image contains no data
			RETURN_WITH_ERROR("Image contains no data", 1);
		}

		// now we are checking whether there is special TVIPS thumbnail information in this image
		if( IFD_Offset >= 0xc010 ) {
			char ca[16];

			ReadFromFile(pBuffer, pos, (BYTE *)&ca[0], 8, bBigEndian);

			if( strncmp( ca, "THUMB128", 8 ) == 0 )	{		// if 0, then we found it
				if( pThumbNail != NULL && pThumbNail->vpData != NULL &&
					pThumbNail->lDimensions == 2 &&
					pThumbNail->lDim[0] == 128 &&
					pThumbNail->lDim[1] == 128 &&
					pThumbNail->lType == DT_RGB8 ) {
					ReadFromFile(pBuffer, pos, (BYTE *)pThumbNail->vpData, 128*128*3, bBigEndian);
				}
			}
		}

		/*
		 * CLU: read all IFD's (Image File Directories)
		 */
		while (IFD_Offset != 0) {
			pos = IFD_Offset;
			/*
			 * CLU: get number of fields
			 */
			WORD	wFields;		
			WORD	wTmp, *wpTmp;

			ReadFromFile(pBuffer, pos, (BYTE *)&wFields, sizeof(WORD), bBigEndian); // num(Directory Entries)
			/*
			 * CLU: read all fields
			 */
			for (int i = 0; i < wFields; i++ ) {

				/*
				 * CLU: get field tag
				 */
				ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Field.tag, sizeof(WORD), bBigEndian);

				/*
				 * CLU: get field type
				 */
				ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Field.type, sizeof(WORD), bBigEndian);

				/*
				 * CLU: get field count
				 */
				ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Field.count, sizeof(DWORD), bBigEndian);

				/*
				 * CLU: get field value
				 */
				ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Field.value, sizeof(DWORD), bBigEndian);

				/*
				 * CLU: decode TIFF fields
				 */
				switch (TIFF_Field.tag) {
					case 256:								// CLU: ImageWidth
						switch (TIFF_Field.type) {
							case TYPE_SHORT:
								if (bBigEndian) {
									wpTmp = (WORD*)&TIFF_Field.value;
									wTmp = wpTmp[1];
									TIFF_Data.ImageWidth = (DWORD)wTmp;
								}
								else {
									wTmp = *( (WORD*)(&TIFF_Field.value)  );
									TIFF_Data.ImageWidth = (DWORD)wTmp;
								}
								break;
							case TYPE_LONG:
								TIFF_Data.ImageWidth = TIFF_Field.value;
								break;

							default:
								RETURN_WITH_ERROR("Cannot decode TIFF data", 1);
						}
						break;

					case 257:								// CLU: ImageLength
						switch (TIFF_Field.type) {
							case TYPE_SHORT:
								if (bBigEndian )
								{
									wpTmp = (WORD*)&TIFF_Field.value;
									wTmp = wpTmp[1];
									TIFF_Data.ImageLength = (DWORD)wTmp;
								}
								else
								{
									wTmp = *( (WORD*)(&TIFF_Field.value) );
									TIFF_Data.ImageLength = (DWORD)wTmp;
								}
								break;
							case TYPE_LONG:
								TIFF_Data.ImageLength = TIFF_Field.value;
								break;

							default:
								RETURN_WITH_ERROR("Cannot decode TIFF data", 1);
						}
						break;

					case 258:														// CLU: BitsPerSample
						if (TIFF_Field.type != TYPE_SHORT) {
							RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}

						if( TRUE == bBigEndian )
						{
							wpTmp = (WORD*)&TIFF_Field.value;
							wTmp = wpTmp[1];
							TIFF_Data.BitsPerSample = (DWORD)wTmp;
						}
						else
							TIFF_Data.BitsPerSample = TIFF_Field.value;


						switch (TIFF_Data.BitsPerSample) {
							case  8:
							case 16:
								break;
							default:
								RETURN_WITH_ERROR("Cannot decode TIFF data", 1);
						}
						break;

					case 259:														// CLU: Compression
						if (TIFF_Field.type != TYPE_SHORT) {
							RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}

						if( bBigEndian )
						{
							wpTmp = (WORD*)&TIFF_Field.value;
							wTmp = wpTmp[1];
							TIFF_Data.Compression = (DWORD)wTmp;
						}
						else
							TIFF_Data.Compression = TIFF_Field.value;
						
						if (TIFF_Data.Compression != 1) {
							RETURN_WITH_ERROR("Cannot decode TIFF data", 1);
						}
						break;

					case 262:														// CLU: PhotometricInterpretation
						if (TIFF_Field.type != TYPE_SHORT) {
							RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}

						if( TRUE == bBigEndian )
						{
							wpTmp = (WORD*)&TIFF_Field.value;
							wTmp = wpTmp[1];
							TIFF_Data.PhotoInterpret= (DWORD)wTmp;
						}
						else
							TIFF_Data.PhotoInterpret = TIFF_Field.value;
							
						switch (TIFF_Data.PhotoInterpret) {
							case 0:
							case 1:
								break;
							default:
								RETURN_WITH_ERROR("Cannot decode TIFF data", 1);
						}
						break;

					case 273:														// CLU: StripOffsets
						switch (TIFF_Field.type) {
							case TYPE_SHORT:
							case TYPE_LONG:
								TIFF_Field.type == TYPE_SHORT ? nStripOffsetsSize = 2 : nStripOffsetsSize = 4;

								TIFF_Data.StripOffsetsCount = TIFF_Field.count;

								if (!(pdwStripOffsets = (DWORD*)malloc(sizeof(DWORD)*TIFF_Data.StripOffsetsCount))) {
									RETURN_WITH_ERROR("No memory", 1);
								}
								
								if (TIFF_Data.StripOffsetsCount == 1) *pdwStripOffsets = (DWORD)TIFF_Field.value;	
								else {
									lFilePosition = pos;
									pos = TIFF_Field.value;

									for (int j = 0; j < (int)TIFF_Data.StripOffsetsCount; j++) {
										ReadFromFile(pBuffer, pos, (BYTE *)&pdwStripOffsets[j], nStripOffsetsSize, bBigEndian);
									}
									pos = lFilePosition;
								}
								break;

							default:
								RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}
						break;

					case 278:														// CLU: RowsPerStrip
						switch (TIFF_Field.type) {
							case TYPE_SHORT:
								if( TRUE == bBigEndian )
								{
									wpTmp = (WORD*)&TIFF_Field.value;
									wTmp = wpTmp[1];
									TIFF_Data.RowsPerStrip = (DWORD)wTmp;
								}
								else
								{
									wTmp = *( (WORD*)(&TIFF_Field.value) );
									TIFF_Data.RowsPerStrip = (DWORD)wTmp;
								}
								break;
							case TYPE_LONG:
								TIFF_Data.RowsPerStrip = TIFF_Field.value;
								break;

							default:
								RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}
						break;

					case 279:														// CLU: StripByteCount
						switch (TIFF_Field.type) {
							case TYPE_SHORT:
							case TYPE_LONG:
								TIFF_Field.type == TYPE_SHORT ? nStripBytesSize = 2 : nStripBytesSize = 4;

								TIFF_Data.StripByteCounts = TIFF_Field.count;

								if (!(pdwBytesPerStrip = (DWORD*)malloc(sizeof(DWORD)*TIFF_Data.StripByteCounts))) {
									RETURN_WITH_ERROR("No memory", 1);
								}

								if (TIFF_Data.StripByteCounts == 1) *pdwBytesPerStrip = TIFF_Field.value;
								else {
									lFilePosition = pos;
									pos = TIFF_Field.value;
									for (int j = 0; j < (int)TIFF_Data.StripByteCounts; j++) {
										ReadFromFile(pBuffer, pos, (BYTE *)&pdwBytesPerStrip[j], nStripBytesSize, bBigEndian);
									}
									pos = lFilePosition;
								} // else
								break;
							default:
								RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}
						break;

					case 282:														// CLU: XResolution
						if (TIFF_Field.type != TYPE_RATIONAL) {
							RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}

						lFilePosition = pos;
						pos = TIFF_Field.value;

						ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Data.XResolutionZ, sizeof(DWORD), bBigEndian);
						ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Data.XResolutionN, sizeof(DWORD), bBigEndian);

						dXResolution = (double)TIFF_Data.XResolutionZ/TIFF_Data.XResolutionN;

						pos = lFilePosition;
						break;

					case 283:														// CLU: YResolution
						if (TIFF_Field.type != TYPE_RATIONAL) {
							RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}

						lFilePosition = pos;
						pos = TIFF_Field.value;

						ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Data.YResolutionZ, sizeof(DWORD), bBigEndian);
						ReadFromFile(pBuffer, pos, (BYTE *)&TIFF_Data.YResolutionN, sizeof(DWORD), bBigEndian);

						dYResolution = (double)TIFF_Data.YResolutionZ/TIFF_Data.YResolutionN;
						pos = lFilePosition;
						break;

					case 296:														// CLU: ResolutionUnit
						if (TIFF_Field.type != TYPE_SHORT) {
							RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}

						TIFF_Data.ResolutionUnit = TIFF_Field.value;

						switch (TIFF_Data.ResolutionUnit) {
							case 1:
							case 2:
							case 3: break;
							default:
								RETURN_WITH_ERROR("Wrong TIFF data type", 1);
						}
						break;

					case TAG_OVL_DATA:												// TVIPS specific Overlay tag
						{
							int lOverlaySize;
							BYTE *cpOverlayData;
							lOverlaySize = (int)TIFF_Field.count;

							if( lOverlaySize == 0 )									// Oops, ignore this tag silently
								break;

							cpOverlayData = new BYTE[lOverlaySize];
							if( cpOverlayData == NULL )								// Oops, no buffer, then overlay data are ignored
								break;

							int lOverlayOffset;
							lOverlayOffset = TIFF_Field.value;

							// store current file position
							lFilePosition = pos;

							pos = lOverlayOffset;
							// seek file position of overlay data

							// read overlay data
							ReadFromFile(pBuffer, pos, (BYTE *)cpOverlayData, lOverlaySize, FALSE /*binary, no byte reversal*/);

							// restore file position
							pos = lFilePosition;

							//g_DataContainer.WriteOverlayData(lImgHandle, cpOverlayData, lOverlaySize);
							delete[] cpOverlayData;
						}
						break;
					case TAG_TEM_DATA:												// TVIPS specific TemData tag
						if (pEmVector) {
							if (TIFF_Field.type != TYPE_LONG) {						// CLU: read only if needed
								RETURN_WITH_ERROR("Wrong TIFF data type", 1);
							}

							int EMVectorVersion;									// CLU: version of EM vector
														
							lFilePosition = pos;;
							pos = TIFF_Field.value;

							ReadFromFile(pBuffer, pos, (BYTE *)&EMVectorVersion, sizeof(DWORD), bBigEndian);

							if (EMVectorVersion >= 1) {
								switch (EMVectorVersion) {
									case 1:											// CLU: em_vector  V == 1 starts here
										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v1.strComment, 80, bBigEndian);

										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lHighTension = (*(int *)cBuf_4)/10000;

										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lSphericalAberration = (*(int *)cBuf_4)/10000;

										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lIllumAperture = (*(int *)cBuf_4)/10000;

										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lElOptMag = (*(int *)cBuf_4)/10000;

										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
//EMVECTOR *pEMV;
//pEMV = (EMVECTOR *)pEmVector->vpData;
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lPostMag = (*(int *)cBuf_4);
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);

										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lFocalLength = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDefocus = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);

										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lAstigmatismNM = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lAstigmatismMRAD = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lBiprismTension = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lSpecimenTiltAngle = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lSpecimenTiltDirection = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lIllumTiltAngle = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lIllumTiltDirection = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lMode = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lEnergySpread = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lChromaticalAberration = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lShutterType = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDefocusSpread = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDNumber = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDPixelXY = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDOffsetX = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDOffsetY = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDSinglePixelSize = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDBinningFactor = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDGain = (*(int *)cBuf_4)/10000;
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDReadOutSpeed = (*(int *)cBuf_4)/10000;
										
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDSensitivity = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lCCDExposureTime = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lFlatfieldCorrection = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDeadPixelCorrection = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lMeanValue = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lStandardDeviation = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDisplacementX = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDisplacementY = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
//										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDate = 0.0001f**(float *)cBuf_4;
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lDate = (float)(*(int *)cBuf_4);
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
//										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lTime = 0.0001f**(float *)cBuf_4;
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lTime = (float)(*(int *)cBuf_4);
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lMinimum = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lMaximum = (*(int *)cBuf_4)/10000;
										
										ReadFromFile(pBuffer, pos, cBuf_4, sizeof(int), bBigEndian);
										(*(EMVECTOR *)pEmVector->vpData).ud_v1.lQualityFactor = (*(int *)cBuf_4)/10000;

										ConvertV1ToV2( (void *)pEmVector->vpData, (void*)&TIFF_Data );
										break;	// case 1:

									case 2:					// CLU: em_vector  V == 2 starts here
										pos += 240;

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szImgName,
															   80*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szImgFolder,
															   80*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgSizeX,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgSizeY,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgSizeZ,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgSizeE,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgDataType,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgCreationDate,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgCreationTime,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szImgComment,
															   512*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szImgHistory,
															   512*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fImgScaling,
															   16*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.cmplxImgStat,
															   16*sizeof(_complex), bBigEndian);
										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgType,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lImgDisplay,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fImgDistX,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fImgDistY,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fImgDistZ,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fImgDistE,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fImgMisc,
															   32*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szTemType,
															   80*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemHT,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemAberr,
															   32*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemEnergy,
															   32*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lTemMode,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemMagScr,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemMagCor,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemMagPst,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lTemStgType,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemStgPos,
															   5*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemImgShift,
															   2*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemBeamShift,
															   2*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemBeamTilt,
															   2*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemTiling,
															   7*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemIllum,
															   3*sizeof(float), bBigEndian);


										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lTemShutter,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fTemMisc,
															   32*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szCamType,
															   80*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamPixel,
															   2*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lCamOffX,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lCamOffY,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lCamBinX,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.lCamBinY,
															   sizeof(int), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamExpTime,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamGain,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamSpeed,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szCamFlat,
															   80*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamSense,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamDose,
															   sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.fCamMisc,
															   32*sizeof(float), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szAdaTietzMicInfo,
															   512*sizeof(USHORT), bBigEndian);

										ReadFromFile(pBuffer, pos, (BYTE *)&(*(EMVECTOR *)pEmVector->vpData).ud_v2.szAdaTietzSpecInfo,
															   512*sizeof(USHORT), bBigEndian);

										int lMagic;								// PSP: Added magic number $aaaaaaaa at the end of the structure
										ReadFromFile(pBuffer, pos, (BYTE *)&lMagic, 4, bBigEndian);

										if( lMagic != 0xaaaaaaaa )
										{
											RETURN_WITH_ERROR("Magic number in EMVector missing", 1);
										}

										break;	// case 2:

											/*
									case 3:  * CLU: add future versions (V == 3) here, and (V > 3) below
									case 4:  * 
									case 5:	 */ 

									default:
										RETURN_WITH_ERROR("Unknown EM vector", 1);

								} // switch (EMVectorVersion)
							}  // if (EMVectorVersion >= 1)

							pos = lFilePosition;
						} // if (pEmVector)
						break; // TAG_TEM_DATA

					default:														// CLU: ignore unknown or not supported tags
						break;

				} // switch (TIFF_Field.tag)
			} // for (int i = 0; i < wFields; i++ )

			if (pDst) {																// CLU: read data if required
				/*
				 * Check the required elements of the structure TIFF_Data against -1.
				 * If there are required elements equal to -1 then return FALSE, if all is OK go on.
				 * Use some assumptions if you can - instead of error message.
				 */
				if (TIFF_Data.ImageWidth	== -1 ||
					TIFF_Data.ImageLength	== -1 ||
					TIFF_Data.BitsPerSample	== -1 ||
					TIFF_Data.Compression	== -1 ||
					TIFF_Data.RowsPerStrip	== -1) {

					free(pdwStripOffsets);											// CLU: free memory
					free(pdwBytesPerStrip);											// CLU: free memory

					RETURN_WITH_ERROR("Required entry in TIFF image missing", 1);
				}

				if (TIFF_Data.PhotoInterpret == -1)  TIFF_Data.PhotoInterpret = 1;

				lFilePosition = pos;

				DWORD 	  nSizeX	= TIFF_Data.ImageWidth,
					  nSizeY	= TIFF_Data.ImageLength,
					  nDataSize	= (TIFF_Data.BitsPerSample / 8);

				/*
				 * CLU: adapt the parameters of the destination variable
				 */
				pDst->lDimensions	= 2L;
				pDst->lDim[0]		= (int)nSizeX;
				pDst->lDim[1]		= (int)nSizeY;
				pDst->lDim[2]		=
				pDst->lDim[3]		= 0L;
				pDst->lType			= nDataSize == 1 ? DT_UCHAR : DT_SHORT;
				
				if ((pDst->vpData = (BYTE *)realloc(pDst->vpData, nSizeX*nSizeY*nDataSize)) == NULL) {
					free(pdwStripOffsets);											// CLU: free memory
					free(pdwBytesPerStrip);											// CLU: free memory

					RETURN_WITH_ERROR("No memory", 1);
				}

				BYTE *pbB = NULL,  *pB = (BYTE *)pDst->vpData;

				for (int j = 0; j < (int)TIFF_Data.StripOffsetsCount; j++) {

					pos = (int)pdwStripOffsets[j];

					if ((pbB = (BYTE *)realloc(pbB, pdwBytesPerStrip[j])) == NULL) {
						free(pdwStripOffsets);										// CLU: free memory
						free(pdwBytesPerStrip);										// CLU: free memory
						RETURN_WITH_ERROR("No memory", 1);
						free(pdwStripOffsets);										// CLU: free memory
					}

					ReadFromFile(pBuffer, pos, pbB, pdwBytesPerStrip[j], FALSE);

					if (nDataSize == 1)
						memcpy(pB, pbB, (int)(pdwBytesPerStrip[j]));

					if (nDataSize > 1) 
					for (int k = 0; k < (int)(pdwBytesPerStrip[j]); k += nDataSize)
						SWITCHENDIAN(pB + k, pbB + k, nDataSize, bBigEndian);		// CLU: interpret data

					pB += pdwBytesPerStrip[j];

				} // for int j = 0; j < (int)TIFF_Data.StripOffsetsCount; j++)
				free(pbB);															// CLU: free memory
			} // if (pDst)

			free(pdwStripOffsets);													// CLU: free memory
			free(pdwBytesPerStrip);													// CLU: free memory

			pos = lFilePosition;

			/*
			 * CLU: check next IFD (if any)
			 *
			 * The next offset should equals to 0, since only 1 image per file is
			 * supported (till 2day (= 21.10.2002))
			 */
			ReadFromFile(pBuffer, pos, (BYTE *)&IFD_Offset, sizeof(ULONG), bBigEndian);

			if (IFD_Offset != 0) {													// CLU: file contains more image(s)
				IFD_Offset = 0;
			}

		} // while (IFD_Offset != 0)
	}
	catch (...) {
		return false;
	}
	return(true);
}


HRESULT CFileTool::DoReadTiff(TCHAR *ptcFilename, struct var_list<BYTE> *pDst, struct var_list<EMVECTOR> *pEmVector, struct var_list<BYTE> *pThumbNail) 
{
	bool hr;
	struct stat FileInfo;
	stat(ptcFilename, &FileInfo);

	int Fsize;
	Fsize = FileInfo.st_size;
	void *map = NULL;

	int fd;
	fd = open(ptcFilename, O_RDONLY);
	if (fd == -1) {
		perror("Error opening file for reading");
		exit(EXIT_FAILURE);
	}
	
	map = mmap(0, Fsize, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED) {
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	}
	
	hr = BufferReadTiff((BYTE *)map, pDst, pEmVector, pThumbNail);
	
	if (munmap(map, Fsize ) == -1) {
		perror("Error un-mmapping the file");
    	}
    	close(fd);
	return 0L;

	/*
	HRESULT hr = S_OK;
	WIN32_FILE_ATTRIBUTE_DATA wfad;
	if (!GetFileAttributesEx(ptcFilename,  GetFileExInfoStandard, &wfad)) {
		RETURN_WITH_ERROR("File don't exist");
	}
	HANDLE hFile = INVALID_HANDLE_VALUE;
	HANDLE hMapFile = NULL;
	void *lpFileData = NULL;
	__try {
		hFile = CreateFile(ptcFilename, GENERIC_READ, 0, NULL, OPEN_EXISTING, 0, NULL);
		if (hFile == INVALID_HANDLE_VALUE) {
			RETURN_WITH_ERROR("Failed to open file");
		}
		hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, wfad.nFileSizeLow, NULL);
		if (hMapFile == NULL) {
			RETURN_WITH_ERROR("Failed to open file");
		}

		lpFileData = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);
		if (lpFileData == NULL) {
			RETURN_WITH_ERROR("Failed to open file");
		}
		hr = BufferReadTiff((BYTE *)lpFileData, pDst, pEmVector, pThumbNail);
	}
	__finally {
		if (lpFileData) 
			UnmapViewOfFile(lpFileData);
		if (hMapFile)
			CloseHandle(hMapFile);
		if (hFile != INVALID_HANDLE_VALUE)
			CloseHandle(hFile);
	}
	return hr;
*/
}
