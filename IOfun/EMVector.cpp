#include "EMVector.h"
#include "stdlib.h"
#include "errno.h"

#define _tcstol strtol
#define _tcstod strtod

//typedef char CHAR;
typedef char *LPSTR;
typedef LPSTR LPTSTR;

CEMVector::CEMVector()
{
	Init(); 
}

CEMVector::~CEMVector()
{
}

CEMVector::CEMVector(const EMVECTOR *pEMV)
{
	Copy(pEMV);
}

CEMVector& CEMVector::Copy(const EMVECTOR *pEMV)
{
	if(pEMV != 0)
	{
		memcpy(&m_EMVector, pEMV, sizeof(struct tagEMVECTOR));
	}

	return *this;
}

CEMVector::CEMVector(CEMVector& other) //cctor
{
	*this = other;
}

CEMVector& CEMVector::operator=(CEMVector& src)
{
	return Copy(&src.m_EMVector);
}

CEMVector& CEMVector::operator=(const EMVECTOR *pEMV)
{
	return Copy(pEMV);
}

EMVECTOR& CEMVector::GetStructEMVECTOR()
{
	return (this->m_EMVector);
}

void CEMVector::ConvertFromCharToUnicode(const char *pSrc, unsigned short *pDest)
{
	char *cp;
	cp = (char*)pDest;

	for(int n=0; n < (int)strlen(pSrc); n++ )
	{
		*cp = pSrc[n];
		cp++;
		*cp=0;
		cp++;
	}
	*cp=0;
	cp++;
	*cp=0;
	cp++;
}

void CEMVector::ConvertFromUnicodeToChar(const unsigned short *pSrc, char *pDest)
{
	char *cp=NULL;

	cp = (char*)pSrc;
	int i=0;
	while( *cp )
	{
		pDest[i] = *cp;
		cp++;
		cp++;
		i++;
	}
}

bool CEMVector::IsEmpty()
{
	bool bEmpty = true;

	bEmpty = (m_EMVector.ud_v2.lImgSizeX				== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgSizeY				== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgSizeZ				== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgSizeE				== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgDataType			== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgCreationDate		== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgCreationTime		== 0);

	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 0]		== 0.0f);		// Scale mode 0 = (contrast/brightness)
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 1]		== 0.0f);		// Contrast ( default )
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 2]		== 0.0f);		// Brightness ( default )
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 3]		== 0.0f);		// Left value
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 4]		== 0.0f);		// Right value
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 5]		== 1.0f);		// Power left scaling
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 6]		== 1.0f);		// Power right scaling
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 7]		== 0.0f);		// Absolute minimum
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 8]		== 0.0f);		// 0.7% Minimum
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[ 9]		== 0.0f);		// 0.7% Maximum
	bEmpty &= (m_EMVector.ud_v2.fImgScaling[10]		== 0.0f);		// Absolute maximum

	int i;
	for( i=11; i<16; i++ )
		bEmpty &= (m_EMVector.ud_v2.fImgScaling[i]		== 0.f);

	bEmpty &= (m_EMVector.ud_v2.lImgType				== 0);
	bEmpty &= (m_EMVector.ud_v2.lImgDisplay			== 0);

	bEmpty &= (m_EMVector.ud_v2.fImgDistX				== 0.f);
	bEmpty &= (m_EMVector.ud_v2.fImgDistY				== 0.f);
	bEmpty &= (m_EMVector.ud_v2.fImgDistZ				== 0.f);
	bEmpty &= (m_EMVector.ud_v2.fImgDistE				== 0.f);

	for( i=0; i<32; i++ )
		bEmpty &= (m_EMVector.ud_v2.fImgMisc[i]			== 0.f);

	bEmpty &= (m_EMVector.ud_v2.fTemHT					== 0.f);
	for( i=0; i<32; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemAberr[i]		== 0.f);
	for( i=0; i<32; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemEnergy[i]		== 0.f);


	bEmpty &= (m_EMVector.ud_v2.lTemMode				== 0);

	bEmpty &= (m_EMVector.ud_v2.fTemMagScr				== 100.f);
	bEmpty &= (m_EMVector.ud_v2.fTemMagCor				==   1.f);
	bEmpty &= (m_EMVector.ud_v2.fTemMagPst				==   1.f);

	bEmpty &= (m_EMVector.ud_v2.lTemStgType				== 0);

	for( i=0; i<5; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemStgPos[i]		== 0.f);
	for( i=0; i<2; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemImgShift[i]		== 0.f);
	for( i=0; i<2; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemBeamShift[i]	== 0.f);
	for( i=0; i<2; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemBeamTilt[i]		== 0.f);
	for( i=0; i<7; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemTiling[i]		== 0.f);
	for( i=0; i<3; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemIllum[i]		== 0.f);

	bEmpty &= (m_EMVector.ud_v2.lTemShutter				== 0);

	for( i=0; i<32; i++ )
		bEmpty &= (m_EMVector.ud_v2.fTemMisc[i]			== 0.f);

//	for( i=0; i<80;i++)
//		m_EMVector.ud_v2.szCamType[i]		=	(unsigned short)0;		// Empty UNICODE string

	for( i=0; i<2; i++ )
		bEmpty &= (m_EMVector.ud_v2.fCamPixel[i]		== 0.f);

	bEmpty &= (m_EMVector.ud_v2.lCamOffX				== 0);
	bEmpty &= (m_EMVector.ud_v2.lCamOffY				== 0);
	bEmpty &= (m_EMVector.ud_v2.lCamBinX				== 0);
	bEmpty &= (m_EMVector.ud_v2.lCamBinY				== 0);

	bEmpty &= (m_EMVector.ud_v2.fCamExpTime				== 0.f);
	bEmpty &= (m_EMVector.ud_v2.fCamGain				== 0.f);
	bEmpty &= (m_EMVector.ud_v2.fCamSpeed				== 0.f);

//	for( i=0; i<80;i++)
//		m_EMVector.ud_v2.szCamFlat[i]		=	(unsigned short)0;		// Empty UNICODE string

	bEmpty &= (m_EMVector.ud_v2.fCamSense				== 0.f);
	bEmpty &= (m_EMVector.ud_v2.fCamDose				== 0.f);

	for( i=0; i<32; i++ )
		bEmpty &= (m_EMVector.ud_v2.fCamMisc[i]			== 0.f);

	return bEmpty;
}

void CEMVector::Init()
{
	m_EMVector.version = EV_VERSION;

	strcpy( m_EMVector.ud_v1.strComment,"empty");
	
	m_EMVector.ud_v1.lHighTension			=	0;
	m_EMVector.ud_v1.lSphericalAberration	=	0;
	m_EMVector.ud_v1.lIllumAperture			=	0;
	m_EMVector.ud_v1.lElOptMag				=	0;
	m_EMVector.ud_v1.lPostMag				=	0;
	m_EMVector.ud_v1.lFocalLength			=	0;
	m_EMVector.ud_v1.lDefocus				=	0;
	m_EMVector.ud_v1.lAstigmatismNM			=	0;
	m_EMVector.ud_v1.lAstigmatismMRAD		=	0;
	m_EMVector.ud_v1.lBiprismTension		=	0;
	m_EMVector.ud_v1.lSpecimenTiltAngle		=	0;
	m_EMVector.ud_v1.lSpecimenTiltDirection	=	0;
	m_EMVector.ud_v1.lIllumTiltDirection	=	0;
	m_EMVector.ud_v1.lIllumTiltAngle		=	0;
	m_EMVector.ud_v1.lMode					=	0;
	m_EMVector.ud_v1.lEnergySpread			=	0;
	m_EMVector.ud_v1.lChromaticalAberration	=	0;
	m_EMVector.ud_v1.lShutterType			=	0;
	m_EMVector.ud_v1.lDefocusSpread			=	0;
	m_EMVector.ud_v1.lCCDNumber				=	0;
	m_EMVector.ud_v1.lCCDPixelXY			=	0;
	m_EMVector.ud_v1.lCCDOffsetX			=	0;
	m_EMVector.ud_v1.lCCDOffsetY			=	0;
	m_EMVector.ud_v1.lCCDSinglePixelSize	=	0;
	m_EMVector.ud_v1.lCCDBinningFactor		=	0;
	m_EMVector.ud_v1.lCCDReadOutSpeed		=	0;
	m_EMVector.ud_v1.lCCDGain				=	0;
	m_EMVector.ud_v1.lCCDSensitivity		=	0;
	m_EMVector.ud_v1.lCCDExposureTime		=	0;
	m_EMVector.ud_v1.lFlatfieldCorrection	=	0;
	m_EMVector.ud_v1.lDeadPixelCorrection	=	0;
	m_EMVector.ud_v1.lMeanValue				=	0;
	m_EMVector.ud_v1.lStandardDeviation		=	0;
	m_EMVector.ud_v1.lDisplacementX			=	0;
	m_EMVector.ud_v1.lDisplacementY			=	0;
	m_EMVector.ud_v1.lDate					=	0;
	m_EMVector.ud_v1.lTime					=	0;
	m_EMVector.ud_v1.lMinimum				=	0;
	m_EMVector.ud_v1.lMaximum				=	0;
	m_EMVector.ud_v1.lQualityFactor			=	0;

	int i;
	for( i=0; i<80;i++)
		m_EMVector.ud_v2.szImgName[i]		=	(unsigned short)0;		// Empty UNICODE string
	for( i=0; i<80;i++)
		m_EMVector.ud_v2.szImgFolder[i]		=	(unsigned short)0;		// Empty UNICODE string

	m_EMVector.ud_v2.lImgSizeX				= 0;
	m_EMVector.ud_v2.lImgSizeY				= 0;
	m_EMVector.ud_v2.lImgSizeZ				= 0;
	m_EMVector.ud_v2.lImgSizeE				= 0;
	m_EMVector.ud_v2.lImgDataType			= 0;
	m_EMVector.ud_v2.lImgCreationDate		= 0;
	m_EMVector.ud_v2.lImgCreationTime		= 0;

	for( i=0; i<512;i++)
	m_EMVector.ud_v2.szImgComment[i]			=	(unsigned short)0;		// Empty UNICODE string
	for( i=0; i<512;i++)
	m_EMVector.ud_v2.szImgHistory[i]			=	(unsigned short)0;		// Empty UNICODE string

	m_EMVector.ud_v2.fImgScaling[ 0]		= 0.0f;		// Scale mode 0 = (contrast/brightness)
	m_EMVector.ud_v2.fImgScaling[ 1]		= 0.0f;		// Contrast ( default )
	m_EMVector.ud_v2.fImgScaling[ 2]		= 0.0f;		// Brightness ( default )
	m_EMVector.ud_v2.fImgScaling[ 3]		= 0.0f;		// Left value
	m_EMVector.ud_v2.fImgScaling[ 4]		= 0.0f;		// Right value
	m_EMVector.ud_v2.fImgScaling[ 5]		= 1.0f;		// Power left scaling
	m_EMVector.ud_v2.fImgScaling[ 6]		= 1.0f;		// Power right scaling
	m_EMVector.ud_v2.fImgScaling[ 7]		= 0.0f;		// Absolute minimum
	m_EMVector.ud_v2.fImgScaling[ 8]		= 0.0f;		// 0.7% Minimum
	m_EMVector.ud_v2.fImgScaling[ 9]		= 0.0f;		// 0.7% Maximum
	m_EMVector.ud_v2.fImgScaling[10]		= 0.0f;		// Absolute maximum

	for( i=11; i<16; i++ )
		m_EMVector.ud_v2.fImgScaling[i]		= 0.f;

	for( i=0; i<16; i++ )
	{
		m_EMVector.ud_v2.cmplxImgStat[i].x		= 0.f;
		m_EMVector.ud_v2.cmplxImgStat[i].y		= 0.f;
	}

	m_EMVector.ud_v2.lImgType				= 0;
	m_EMVector.ud_v2.lImgDisplay				= 0;

	m_EMVector.ud_v2.fImgDistX				= 0.f;
	m_EMVector.ud_v2.fImgDistY				= 0.f;
	m_EMVector.ud_v2.fImgDistZ				= 0.f;
	m_EMVector.ud_v2.fImgDistE				= 0.f;

	for( i=0; i<32; i++ )
		m_EMVector.ud_v2.fImgMisc[i]			= 0.f;

	m_EMVector.ud_v2.fImgMisc[0] = 0.7f;
	m_EMVector.ud_v2.fImgMisc[1] = 0.7f;

	for( i=0; i<80;i++)
		m_EMVector.ud_v2.szTemType[i]		=	(unsigned short)0;		// Empty UNICODE string

	m_EMVector.ud_v2.fTemHT					= 0.f;
	for( i=0; i<32; i++ )
		m_EMVector.ud_v2.fTemAberr[i]		= 0.f;
	for( i=0; i<32; i++ )
		m_EMVector.ud_v2.fTemEnergy[i]		= 0.f;


	m_EMVector.ud_v2.lTemMode				= 0;

	m_EMVector.ud_v2.fTemMagScr				= 100.f;
	m_EMVector.ud_v2.fTemMagCor				=   1.f;
	m_EMVector.ud_v2.fTemMagPst				=   1.f;

	m_EMVector.ud_v2.lTemStgType				= 0;

	for( i=0; i<5; i++ )
		m_EMVector.ud_v2.fTemStgPos[i]		= 0.f;
	for( i=0; i<2; i++ )
		m_EMVector.ud_v2.fTemImgShift[i]		= 0.f;
	for( i=0; i<2; i++ )
		m_EMVector.ud_v2.fTemBeamShift[i]	= 0.f;
	for( i=0; i<2; i++ )
		m_EMVector.ud_v2.fTemBeamTilt[i]		= 0.f;
	for( i=0; i<7; i++ )
		m_EMVector.ud_v2.fTemTiling[i]		= 0.f;
	for( i=0; i<3; i++ )
		m_EMVector.ud_v2.fTemIllum[i]		= 0.f;

	m_EMVector.ud_v2.lTemShutter				= 0;

	for( i=0; i<32; i++ )
		m_EMVector.ud_v2.fTemMisc[i]			= 0.f;

	for( i=0; i<80;i++)
		m_EMVector.ud_v2.szCamType[i]		=	(unsigned short)0;		// Empty UNICODE string

	for( i=0; i<2; i++ )
		m_EMVector.ud_v2.fCamPixel[i]		= 0.f;

	m_EMVector.ud_v2.lCamOffX				= 0;
	m_EMVector.ud_v2.lCamOffY				= 0;
	m_EMVector.ud_v2.lCamBinX				= 0;
	m_EMVector.ud_v2.lCamBinY				= 0;

	m_EMVector.ud_v2.fCamExpTime				= 0.f;
	m_EMVector.ud_v2.fCamGain				= 0.f;
	m_EMVector.ud_v2.fCamSpeed				= 0.f;

	for( i=0; i<80;i++)
		m_EMVector.ud_v2.szCamFlat[i]		=	(unsigned short)0;		// Empty UNICODE string

	m_EMVector.ud_v2.fCamSense				= 0.f;
	m_EMVector.ud_v2.fCamDose				= 0.f;

	for( i=0; i<32; i++ )
		m_EMVector.ud_v2.fCamMisc[i]			= 0.f;

	for( i=0; i<512;i++)
		m_EMVector.ud_v2.szAdaTietzMicInfo[i]	=	(unsigned short)0;		// Empty UNICODE string
	for( i=0; i<512;i++)
		m_EMVector.ud_v2.szAdaTietzSpecInfo[i]	=	(unsigned short)0;		// Empty UNICODE string
}

bool CEMVector::operator==(const CEMVector &other)
{
	bool bEqual;

	bEqual = (m_EMVector.version == other.m_EMVector.version);

	//strcpy( m_EMVector.ud_v1.strComment,"empty");
	
	//ud_v1 cannot be changed, it is used only for compatibility with 3.0 => make no sense to check it
	/*
	bEqual &= (m_EMVector.ud_v1.lHighTension			== other.m_EMVector.ud_v1.lHighTension );
	bEqual &= (m_EMVector.ud_v1.lSphericalAberration	== other.m_EMVector.ud_v1.lSphericalAberration);
	bEqual &= (m_EMVector.ud_v1.lIllumAperture			== other.m_EMVector.ud_v1.lIllumAperture);
	bEqual &= (m_EMVector.ud_v1.lElOptMag				== other.m_EMVector.ud_v1.lElOptMag);
	bEqual &= (m_EMVector.ud_v1.lPostMag				== other.m_EMVector.ud_v1.lPostMag);
	bEqual &= (m_EMVector.ud_v1.lFocalLength			== other.m_EMVector.ud_v1.lFocalLength);
	bEqual &= (m_EMVector.ud_v1.lDefocus				== other.m_EMVector.ud_v1.lDefocus);
	bEqual &= (m_EMVector.ud_v1.lAstigmatismNM			== other.m_EMVector.ud_v1.lAstigmatismNM);
	bEqual &= (m_EMVector.ud_v1.lAstigmatismMRAD		== other.m_EMVector.ud_v1.lAstigmatismMRAD);
	bEqual &= (m_EMVector.ud_v1.lBiprismTension			== other.m_EMVector.ud_v1.lBiprismTension);
	bEqual &= (m_EMVector.ud_v1.lSpecimenTiltAngle		== other.m_EMVector.ud_v1.lSpecimenTiltAngle);
	bEqual &= (m_EMVector.ud_v1.lSpecimenTiltDirection	== other.m_EMVector.ud_v1.lSpecimenTiltDirection);
	bEqual &= (m_EMVector.ud_v1.lIllumTiltDirection		== other.m_EMVector.ud_v1.lIllumTiltDirection);
	bEqual &= (m_EMVector.ud_v1.lIllumTiltAngle			== other.m_EMVector.ud_v1.lIllumTiltAngle);
	bEqual &= (m_EMVector.ud_v1.lMode					== other.m_EMVector.ud_v1.lMode);
	bEqual &= (m_EMVector.ud_v1.lEnergySpread			== other.m_EMVector.ud_v1.lEnergySpread);
	bEqual &= (m_EMVector.ud_v1.lChromaticalAberration	== other.m_EMVector.ud_v1.lChromaticalAberration);
	bEqual &= (m_EMVector.ud_v1.lShutterType			== other.m_EMVector.ud_v1.lShutterType);
	bEqual &= (m_EMVector.ud_v1.lDefocusSpread			== other.m_EMVector.ud_v1.lDefocusSpread);
	bEqual &= (m_EMVector.ud_v1.lCCDNumber				== other.m_EMVector.ud_v1.lCCDNumber);
	bEqual &= (m_EMVector.ud_v1.lCCDPixelXY				== other.m_EMVector.ud_v1.lCCDPixelXY);
	bEqual &= (m_EMVector.ud_v1.lCCDOffsetX				== other.m_EMVector.ud_v1.lCCDOffsetX);
	bEqual &= (m_EMVector.ud_v1.lCCDOffsetY				== other.m_EMVector.ud_v1.lCCDOffsetY);
	bEqual &= (m_EMVector.ud_v1.lCCDSinglePixelSize		== other.m_EMVector.ud_v1.lCCDSinglePixelSize);
	bEqual &= (m_EMVector.ud_v1.lCCDBinningFactor		== other.m_EMVector.ud_v1.lCCDBinningFactor);
	bEqual &= (m_EMVector.ud_v1.lCCDReadOutSpeed		== other.m_EMVector.ud_v1.lCCDReadOutSpeed);
	bEqual &= (m_EMVector.ud_v1.lCCDGain				== other.m_EMVector.ud_v1.lCCDGain);
	bEqual &= (m_EMVector.ud_v1.lCCDSensitivity			== other.m_EMVector.ud_v1.lCCDSensitivity);
	bEqual &= (m_EMVector.ud_v1.lCCDExposureTime		== other.m_EMVector.ud_v1.lCCDExposureTime);
	bEqual &= (m_EMVector.ud_v1.lFlatfieldCorrection	== other.m_EMVector.ud_v1.lFlatfieldCorrection);
	bEqual &= (m_EMVector.ud_v1.lDeadPixelCorrection	== other.m_EMVector.ud_v1.lDeadPixelCorrection);
	bEqual &= (m_EMVector.ud_v1.lMeanValue				== other.m_EMVector.ud_v1.lMeanValue);
	bEqual &= (m_EMVector.ud_v1.lStandardDeviation		== other.m_EMVector.ud_v1.lStandardDeviation);
	bEqual &= (m_EMVector.ud_v1.lDisplacementX			== other.m_EMVector.ud_v1.lDisplacementX);
	bEqual &= (m_EMVector.ud_v1.lDisplacementY			== other.m_EMVector.ud_v1.lDisplacementY);
	bEqual &= (m_EMVector.ud_v1.lDate					== other.m_EMVector.ud_v1.lDate);
	bEqual &= (m_EMVector.ud_v1.lTime					== other.m_EMVector.ud_v1.lTime);
	bEqual &= (m_EMVector.ud_v1.lMinimum				== other.m_EMVector.ud_v1.lMinimum);
	bEqual &= (m_EMVector.ud_v1.lMaximum				== other.m_EMVector.ud_v1.lMaximum);
	bEqual &= (m_EMVector.ud_v1.lQualityFactor			== other.m_EMVector.ud_v1.lQualityFactor);
	*/

	int i;
	//for( i=0; i<80;i++)
	//	bEqual &= (m_EMVector.ud_v2.szImgName[i]		== other.m_EMVector.ud_v2.szImgName[i]);		// UNICODE string
	//for( i=0; i<80;i++)
	//	bEqual &= (m_EMVector.ud_v2.szImgFolder[i]		== other.m_EMVector.ud_v2.szImgFolder[i]);		// UNICODE string

	bEqual &= (m_EMVector.ud_v2.lImgSizeX				== other.m_EMVector.ud_v2.lImgSizeX);
	bEqual &= (m_EMVector.ud_v2.lImgSizeY				== other.m_EMVector.ud_v2.lImgSizeY);
	bEqual &= (m_EMVector.ud_v2.lImgSizeZ				== other.m_EMVector.ud_v2.lImgSizeZ);
	bEqual &= (m_EMVector.ud_v2.lImgSizeE				== other.m_EMVector.ud_v2.lImgSizeE);
	bEqual &= (m_EMVector.ud_v2.lImgDataType			== other.m_EMVector.ud_v2.lImgDataType);
	bEqual &= (m_EMVector.ud_v2.lImgCreationDate		== other.m_EMVector.ud_v2.lImgCreationDate);
	bEqual &= (m_EMVector.ud_v2.lImgCreationTime		== other.m_EMVector.ud_v2.lImgCreationTime);

/*	for( i=0; i<512;i++)
		bEqual &= m_EMVector.ud_v2.szImgComment[i]		== other.m_EMVector;		// Empty UNICODE string
	for( i=0; i<512;i++)
		bEqual &= m_EMVector.ud_v2.szImgHistory[i]		== other.m_EMVector;		// Empty UNICODE string
*/
	for(i=0; i<16; i++)
		bEqual &= (m_EMVector.ud_v2.fImgScaling[i] == other.m_EMVector.ud_v2.fImgScaling[i]);

	//for( i=0; i<16; i++ )
	//{
	//	bEqual &= (m_EMVector.ud_v2.cmplxImgStat[i].x == other.m_EMVector.ud_v2.cmplxImgStat[i].x);
	//	bEqual &= (m_EMVector.ud_v2.cmplxImgStat[i].y == other.m_EMVector.ud_v2.cmplxImgStat[i].y);
	//}

	bEqual &= (m_EMVector.ud_v2.lImgType			== other.m_EMVector.ud_v2.lImgType);
	bEqual &= (m_EMVector.ud_v2.lImgDisplay			== other.m_EMVector.ud_v2.lImgDisplay);

	bEqual &= (m_EMVector.ud_v2.fImgDistX == other.m_EMVector.ud_v2.fImgDistX);
	bEqual &= (m_EMVector.ud_v2.fImgDistY == other.m_EMVector.ud_v2.fImgDistY);
	bEqual &= (m_EMVector.ud_v2.fImgDistZ == other.m_EMVector.ud_v2.fImgDistZ);
	bEqual &= (m_EMVector.ud_v2.fImgDistE == other.m_EMVector.ud_v2.fImgDistE);

	for( i=0; i<32; i++ )
		bEqual &= (m_EMVector.ud_v2.fImgMisc[i] == other.m_EMVector.ud_v2.fImgMisc[i]);
/*
	for( i=0; i<80;i++)
		bEqual &= (m_EMVector.ud_v2.szTemType[i] == other.m_EMVector.ud_v2.szTemType[i]);		// UNICODE string
*/
	bEqual &= (m_EMVector.ud_v2.fTemHT == other.m_EMVector.ud_v2.fTemHT);

	for( i=0; i<32; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemAberr[i]		== other.m_EMVector.ud_v2.fTemAberr[i]);
	for( i=0; i<32; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemEnergy[i]		== other.m_EMVector.ud_v2.fTemEnergy[i]);

	bEqual &= (m_EMVector.ud_v2.lTemMode				== other.m_EMVector.ud_v2.lTemMode);

	bEqual &= (m_EMVector.ud_v2.fTemMagScr				== other.m_EMVector.ud_v2.fTemMagScr);
	bEqual &= (m_EMVector.ud_v2.fTemMagCor				== other.m_EMVector.ud_v2.fTemMagCor);
	bEqual &= (m_EMVector.ud_v2.fTemMagPst				== other.m_EMVector.ud_v2.fTemMagPst);

	bEqual &= (m_EMVector.ud_v2.lTemStgType			== other.m_EMVector.ud_v2.lTemStgType);

	for( i=0; i<5; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemStgPos[i]		== other.m_EMVector.ud_v2.fTemStgPos[i]);
	for( i=0; i<2; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemImgShift[i]	== other.m_EMVector.ud_v2.fTemImgShift[i]);
	for( i=0; i<2; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemBeamShift[i]	== other.m_EMVector.ud_v2.fTemBeamShift[i]);
	for( i=0; i<2; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemBeamTilt[i]	== other.m_EMVector.ud_v2.fTemBeamTilt[i]);
	for( i=0; i<7; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemTiling[i]		== other.m_EMVector.ud_v2.fTemTiling[i]);
	for( i=0; i<3; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemIllum[i]		== other.m_EMVector.ud_v2.fTemIllum[i]);

	bEqual &= (m_EMVector.ud_v2.lTemShutter			== other.m_EMVector.ud_v2.lTemShutter);

	for( i=0; i<32; i++ )
		bEqual &= (m_EMVector.ud_v2.fTemMisc[i]		== other.m_EMVector.ud_v2.fTemMisc[i]);

/*	for( i=0; i<80;i++)
		bEqual &= (m_EMVector.ud_v2.szCamType[i]	== other.m_EMVector.ud_v2.szCamType[i]);		// UNICODE string
*/
	for( i=0; i<2; i++ )
		bEqual &= (m_EMVector.ud_v2.fCamPixel[i]	== other.m_EMVector.ud_v2.fCamPixel[i]);

	bEqual &= (m_EMVector.ud_v2.lCamOffX				== other.m_EMVector.ud_v2.lCamOffX);
	bEqual &= (m_EMVector.ud_v2.lCamOffY				== other.m_EMVector.ud_v2.lCamOffY);
	bEqual &= (m_EMVector.ud_v2.lCamBinX				== other.m_EMVector.ud_v2.lCamBinX);
	bEqual &= (m_EMVector.ud_v2.lCamBinY				== other.m_EMVector.ud_v2.lCamBinY);

	bEqual &= (m_EMVector.ud_v2.fCamExpTime			== other.m_EMVector.ud_v2.fCamExpTime);
	bEqual &= (m_EMVector.ud_v2.fCamGain			== other.m_EMVector.ud_v2.fCamGain);
	bEqual &= (m_EMVector.ud_v2.fCamSpeed			== other.m_EMVector.ud_v2.fCamSpeed);

/*	for( i=0; i<80;i++)
		bEqual &= (m_EMVector.ud_v2.szCamFlat[i]		== other.m_EMVector.ud_v2.szCamFlat[i]);		// UNICODE string
*/
	bEqual &= (m_EMVector.ud_v2.fCamSense				== other.m_EMVector.ud_v2.fCamSense);
	bEqual &= (m_EMVector.ud_v2.fCamDose				== other.m_EMVector.ud_v2.fCamDose);

	for( i=0; i<32; i++ )
		bEqual &= (m_EMVector.ud_v2.fCamMisc[i]		== other.m_EMVector.ud_v2.fCamMisc[i]);

/*	for( i=0; i<512;i++)
		m_EMVector.ud_v2.szAdaTietzMicInfo[i]== other.m_EMVector;		// Empty UNICODE string
	for( i=0; i<512;i++)
		m_EMVector.ud_v2.szAdaTietzSpecInfo[i]== other.m_EMVector;		// Empty UNICODE string
*/
	return bEqual;
}

void CEMVector::SetImageName(const string &strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy((TCHAR*)m_EMVector.ud_v2.szImgName, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szImgName);
#endif
}

string CEMVector::GetImageName()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szImgName)/sizeof(TCHAR), '\0');

	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szImgName);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szImgName, pVal);
#endif

	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void CEMVector::SetCameraFlat(const string& strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy(m_EMVector.ud_v2.szCamFlat, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szCamFlat);
#endif
}

string CEMVector::GetCameraFlat()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szCamFlat)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szCamFlat);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szCamFlat, pVal);
#endif
	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void CEMVector::SetImageFolder(const string& strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy(m_EMVector.ud_v2.szImgFolder, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szImgFolder);
#endif
	//InformObservers(false);
}

string CEMVector::GetImageFolder()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szImgFolder)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szImgFolder);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szImgFolder, pVal);
#endif
	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void CEMVector::SetImageSizeX(int lVal)
{
	m_EMVector.ud_v2.lImgSizeX = lVal;
	//InformObservers(false);
}

int CEMVector::GetImageSizeX() const
{
	return m_EMVector.ud_v2.lImgSizeX;
}

void		CEMVector::SetImageSizeY(int lVal)
{
	m_EMVector.ud_v2.lImgSizeY = lVal;
	//InformObservers(false);
}

int		CEMVector::GetImageSizeY() const
{
	return m_EMVector.ud_v2.lImgSizeY;
}

void		CEMVector::SetImageSizeZ(int lVal)
{
	m_EMVector.ud_v2.lImgSizeZ = lVal;
	//InformObservers(false);
}

int		CEMVector::GetImageSizeZ() const
{
	return m_EMVector.ud_v2.lImgSizeZ;
}

void		CEMVector::SetImageSizeE(int lVal)
{
	m_EMVector.ud_v2.lImgSizeE = lVal;
	//InformObservers(false);
}

int		CEMVector::GetImageSizeE() const
{
	return m_EMVector.ud_v2.lImgSizeE;
}

void		CEMVector::SetImageDataType(IMAGE_TYPE imgType)
{
	m_EMVector.ud_v2.lImgDataType = (int)imgType;
	//InformObservers(false);
}

IMAGE_TYPE	CEMVector::GetImageDataType() const
{
	return (IMAGE_TYPE)m_EMVector.ud_v2.lImgDataType;
}

void		CEMVector::SetImageCreationDate(int lVal)
{
	m_EMVector.ud_v2.lImgCreationDate = lVal;
	//InformObservers(false);
}

int		CEMVector::GetImageCreationDate() const
{
	return m_EMVector.ud_v2.lImgCreationDate;
}

void		CEMVector::SetImageCreationTime(int lVal)
{
	m_EMVector.ud_v2.lImgCreationTime = lVal;
	//InformObservers(false);
}

int		CEMVector::GetImageCreationTime() const
{
	return m_EMVector.ud_v2.lImgCreationTime;
}

void		CEMVector::SetImageComment(const string &strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy(m_EMVector.ud_v2.szImgComment, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szImgComment);
#endif
	//InformObservers(false);
}

string		CEMVector::GetImageComment()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szImgComment)/sizeof(TCHAR), '\0');

	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szImgComment);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szImgComment, pVal);
#endif

	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void		CEMVector::SetImageHistory(const string &strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy(m_EMVector.ud_v2.szImgHistory, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szImgHistory);
#endif
	//InformObservers(false);
}

string		CEMVector::GetImageHistory()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szImgHistory)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szImgHistory);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szImgHistory, pVal);
#endif

	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void		CEMVector::SetScaleMode(float fVal)
{
	if(fVal < 0)
		return;

	m_EMVector.ud_v2.fImgScaling[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetScaleMode() const // Scale mode 0 = (contrast/brightness)
{
	return m_EMVector.ud_v2.fImgScaling[0];
}

bool		CEMVector::GetScaleModes(float *buf, int nNums) const
{
	if(nNums <0 || nNums > 16)
		return false;

	for(int i=0; i<nNums; i++)
		buf[i] = m_EMVector.ud_v2.fImgScaling[i];

	return true;
}

void		CEMVector::SetContrast(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetContrast() const
{
	return m_EMVector.ud_v2.fImgScaling[1];
}

void		CEMVector::SetBrightness(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[2] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetBrightness() const
{
	return m_EMVector.ud_v2.fImgScaling[2];
}

void		CEMVector::SetOverflow(float fVal)
{
	m_EMVector.ud_v2.fImgMisc[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetOverflow() const
{
	return m_EMVector.ud_v2.fImgMisc[0];
}

void		CEMVector::SetUnderflow(float fVal)
{
	m_EMVector.ud_v2.fImgMisc[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetUnderflow() const
{
	return m_EMVector.ud_v2.fImgMisc[1];
}

void		CEMVector::SetMinScalingValue(float fVal)
{
	m_EMVector.ud_v2.fImgMisc[2] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetMinScalingValue() const
{
	return m_EMVector.ud_v2.fImgMisc[2];
}

void		CEMVector::SetMaxScalingValue(float fVal)
{
	m_EMVector.ud_v2.fImgMisc[3] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetMaxScalingValue() const
{
	return m_EMVector.ud_v2.fImgMisc[3];
}

void		CEMVector::SetNavigatorPositionX(float fVal)
{
	m_EMVector.ud_v2.fImgMisc[4] = fVal;
}

float		CEMVector::GetNavigatorPositionX() const
{
	return m_EMVector.ud_v2.fImgMisc[4];
}

void		CEMVector::SetNavigatorPositionY(float fVal)
{
	m_EMVector.ud_v2.fImgMisc[5] = fVal;
}

float		CEMVector::GetNavigatorPositionY() const
{
	return m_EMVector.ud_v2.fImgMisc[5];
}

void		CEMVector::SetLeftScalingValue(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[3] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetLeftScalingValue() const
{
	return m_EMVector.ud_v2.fImgScaling[3];
}

void		CEMVector::SetRightScalingValue(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[4] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetRightScalingValue() const
{
	return m_EMVector.ud_v2.fImgScaling[4];
}

void		CEMVector::SetLeftPowerScalingValue(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[5] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetLeftPowerScalingValue() const
{
	return m_EMVector.ud_v2.fImgScaling[5];
}

void		CEMVector::SetRightPowerScalingValue(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[6] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetRightPowerScalingValue() const
{
	return m_EMVector.ud_v2.fImgScaling[6];
}

void		CEMVector::SetScalingAbsoluteMinimum(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[7] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetScalingAbsoluteMinimum() const
{
	return m_EMVector.ud_v2.fImgScaling[7];
}

void		CEMVector::SetScalingMinimum(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[8] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetScalingMinimum() const
{
	return m_EMVector.ud_v2.fImgScaling[8];
}

void		CEMVector::SetScalingMaximum(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[9] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetScalingMaximum() const
{
	return m_EMVector.ud_v2.fImgScaling[9];
}

void		CEMVector::SetScalingAbsoluteMaximum(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[10] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetScalingAbsoluteMaximum() const
{
	return m_EMVector.ud_v2.fImgScaling[10];
}

void		CEMVector::SetImageType(int nVal)
{
	m_EMVector.ud_v2.lImgType = (int)nVal;
	//InformObservers(false);
}

int			CEMVector::GetImageType() const
{
	return (int)m_EMVector.ud_v2.lImgType;
}

void		CEMVector::SetDisplayType(int lVal)
{
	m_EMVector.ud_v2.lImgDisplay = lVal;
	//InformObservers(false);
}

int		CEMVector::GetDisplayType() const
{
	return m_EMVector.ud_v2.lImgDisplay;
}

void		CEMVector::SetImageDistX(float fVal)
{
	m_EMVector.ud_v2.fImgDistX = fVal;

	//InformObservers(true);
}

float		CEMVector::GetImageDistX() const
{
	return m_EMVector.ud_v2.fImgDistX;
}

void		CEMVector::SetImageDistY(float fVal)
{
	m_EMVector.ud_v2.fImgDistY = fVal;
	//InformObservers(true);
}

float		CEMVector::GetImageDistY() const
{
	return m_EMVector.ud_v2.fImgDistY;
}

void		CEMVector::SetImageDistZ(float fVal)
{
	m_EMVector.ud_v2.fImgDistZ= fVal;
	//InformObservers(false);
}

float		CEMVector::GetImageDistZ() const
{
	return m_EMVector.ud_v2.fImgDistZ;
}

void		CEMVector::SetImageDistE(float fVal)
{
	m_EMVector.ud_v2.fImgDistE = fVal;
	//InformObservers(false);
}

float		CEMVector::GetImageDistE() const
{
	return m_EMVector.ud_v2.fImgDistE;
}

void		CEMVector::SetTEMType(const string &strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy(m_EMVector.ud_v2.szTemType, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szTemType);
#endif
	//InformObservers(false);
}

string		CEMVector::GetTEMType()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szTemType)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szTemType);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szTemType, pVal);
#endif

	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void		CEMVector::SetTEMHighTension(float fVal)
{
	if(fVal < 0)
		return;

	m_EMVector.ud_v2.fTemHT = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMHighTension() const
{
	return m_EMVector.ud_v2.fTemHT;
}



void		CEMVector::SetTEMMode(int lTemMode)
{
	m_EMVector.ud_v2.lTemMode = lTemMode;
	//InformObservers(false);
}

int		CEMVector::GetTEMMode() const
{
	return m_EMVector.ud_v2.lTemMode;
}

void		CEMVector::SetCurrentMagnificationRelativeToScreen(float fVal)
{
	if(fVal <= 0)
		return;

	m_EMVector.ud_v2.fTemMagScr = fVal;

	//InformObservers(true);

}

float		CEMVector::GetCurrentMagnificationRelativeToScreen() const
{
	return m_EMVector.ud_v2.fTemMagScr;
}

void		CEMVector::SetCorrectedMagnification(float fVal)
{
	if(fVal <= 0)
		return;

	m_EMVector.ud_v2.fTemMagCor = fVal;
	//InformObservers(true);
}

float		CEMVector::GetCorrectedMagnification() const
{
	return m_EMVector.ud_v2.fTemMagCor;
}

void		CEMVector::SetPostMagnification(float fVal)
{
	if(fVal <= 0)
		return;

	m_EMVector.ud_v2.fTemMagPst = fVal;
	//InformObservers(true);
}

float		CEMVector::GetPostMagnification() const
{
	return m_EMVector.ud_v2.fTemMagPst;
}

void		CEMVector::SetStageType(int lVal)
{
	m_EMVector.ud_v2.lTemStgType = lVal;
	//InformObservers(false);
}

int		CEMVector::GetStageType() const
{
	return m_EMVector.ud_v2.lTemStgType;
}

bool		CEMVector::GetStgPos(float *pBuf, int nNums)
{
	if(nNums < 0 || nNums > 5)
		return false;

	for(int i=0; i < nNums; i++)
		pBuf[i] = m_EMVector.ud_v2.fTemStgPos[i];

	return true;
}

void		CEMVector::SetGonioX(float fVal)
{
	m_EMVector.ud_v2.fTemStgPos[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetGonioX() const
{
	return m_EMVector.ud_v2.fTemStgPos[0];
}

void		CEMVector::SetGonioY(float fVal)
{
	m_EMVector.ud_v2.fTemStgPos[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetGonioY() const
{
	return m_EMVector.ud_v2.fTemStgPos[1];
}

void		CEMVector::SetGonioZ(float fVal)
{
	m_EMVector.ud_v2.fTemStgPos[2] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetGonioZ() const
{
	return m_EMVector.ud_v2.fTemStgPos[2];
}

void		CEMVector::SetGonioTiltA(float fVal)
{
	if(fVal < -360 || fVal > 360)
		return;

	m_EMVector.ud_v2.fTemStgPos[3] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetGonioTiltA() const
{
	return m_EMVector.ud_v2.fTemStgPos[3];
}

void		CEMVector::SetGonioTiltB(float fVal)
{
	if(fVal < -360 || fVal > 360)
		return;

	m_EMVector.ud_v2.fTemStgPos[4] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetGonioTiltB() const
{
	return m_EMVector.ud_v2.fTemStgPos[4];
}

void		CEMVector::SetTEMTilingNrImgX(float fVal)
{
	if(fVal < 1)
		return;

	m_EMVector.ud_v2.fTemTiling[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMTilingNrImgX() const
{
	return m_EMVector.ud_v2.fTemTiling[1];
}

void		CEMVector::SetTEMTilingNrImgY(float fVal)
{
	if(fVal < 1)
		return;

	m_EMVector.ud_v2.fTemTiling[2] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMTilingNrImgY() const
{
	return m_EMVector.ud_v2.fTemTiling[2];
}

void		CEMVector::SetTEMTilingOverlapX(float fVal)
{
	m_EMVector.ud_v2.fTemTiling[5] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMTilingOverlapX() const
{
	return m_EMVector.ud_v2.fTemTiling[5];
}

void		CEMVector::SetTEMTilingOverlapY(float fVal)
{
	m_EMVector.ud_v2.fTemTiling[6] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMTilingOverlapY() const
{
	return m_EMVector.ud_v2.fTemTiling[6];
}

void		CEMVector::SetTEMTilingType(int nVal)
{
	if(nVal <0 || nVal >1)
		return;

	m_EMVector.ud_v2.fTemTiling[0] = (float)nVal;
	//InformObservers(false);
}

int			CEMVector::GetTEMTilingType() const
{
	return (int)m_EMVector.ud_v2.fTemTiling[0];
}

void		CEMVector::SetTEMTilingMaxNumberX(float fVal)
{
	m_EMVector.ud_v2.fTemTiling[3] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMTilingMaxNumberX() const
{
	return m_EMVector.ud_v2.fTemTiling[3];
}

void		CEMVector::SetTEMTilingMaxNumberY(float fVal)
{
	m_EMVector.ud_v2.fTemTiling[4] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMTilingMaxNumberY() const
{
	return m_EMVector.ud_v2.fTemTiling[4];
}

void		CEMVector::SetTEMIlluminationSpotSize(float fVal)
{
	m_EMVector.ud_v2.fTemIllum[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMIlluminationSpotSize() const
{
	return m_EMVector.ud_v2.fTemIllum[0];
}

void		CEMVector::SetTEMIlluminationIntensity(float fVal)
{
	m_EMVector.ud_v2.fTemIllum[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMIlluminationIntensity() const
{
	return m_EMVector.ud_v2.fTemIllum[1];
}

void		CEMVector::SetTEMShutter(int lVal)
{
	m_EMVector.ud_v2.lTemShutter = lVal;
	//InformObservers(false);
}

int		CEMVector::GetTEMShutter() const
{
	return m_EMVector.ud_v2.lTemShutter;
}

void		CEMVector::SetCameraType(const string &strVal)
{
	if(strVal.size() == 0)
		return;

	CStdString csValue(strVal);

#ifdef _UNICODE
	//wcscpy(m_EMVector.ud_v2.szCamType, csValue.GetBuffer(csValue.GetLength()));
#else
	ConvertFromCharToUnicode(csValue.GetBuffer(csValue.GetLength()), m_EMVector.ud_v2.szCamType);
#endif
	//InformObservers(false);
}

string		CEMVector::GetCameraType()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szCamType)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szCamType);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szCamType, pVal);
#endif
	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

void		CEMVector::SetCameraPixelSizeX(float fVal)
{
	m_EMVector.ud_v2.fCamPixel[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraPixelSizeX() const
{
	return m_EMVector.ud_v2.fCamPixel[0];
}

void		CEMVector::SetCameraPixelSizeY(float fVal)
{
	m_EMVector.ud_v2.fCamPixel[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraPixelSizeY() const
{
	return m_EMVector.ud_v2.fCamPixel[1];
}

void		CEMVector::SetCameraOffsetX(int lVal)
{
	m_EMVector.ud_v2.lCamOffX = lVal;
	//InformObservers(false);
}

int		CEMVector::GetCameraOffsetX() const
{
	return m_EMVector.ud_v2.lCamOffX;
}

void		CEMVector::SetCameraOffsetY(int lVal)
{
	m_EMVector.ud_v2.lCamOffY = lVal;
	//InformObservers(false);
}

int		CEMVector::GetCameraOffsetY() const
{
	return m_EMVector.ud_v2.lCamOffY;
}

void		CEMVector::SetCameraBinningX(int lVal)
{
	m_EMVector.ud_v2.lCamBinX = lVal;
	//InformObservers(true);
}

int		CEMVector::GetCameraBinningX() const
{
	return m_EMVector.ud_v2.lCamBinX;
}

void		CEMVector::SetCameraBinningY(int lVal)
{
	m_EMVector.ud_v2.lCamBinY = lVal;
	//InformObservers(false);
}

int		CEMVector::GetCameraBinningY() const
{
	return m_EMVector.ud_v2.lCamBinY;
}

void		CEMVector::SetCameraExposureTime(float fVal)
{
	m_EMVector.ud_v2.fCamExpTime = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraExposureTime() const
{
	return m_EMVector.ud_v2.fCamExpTime;
}

void		CEMVector::SetCameraGainFactor(float fVal)
{
	m_EMVector.ud_v2.fCamGain = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraGainFactor() const
{
	return m_EMVector.ud_v2.fCamGain;
}

void		CEMVector::SetCameraReadoutRate(float fVal)
{
	m_EMVector.ud_v2.fCamSpeed = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraReadoutRate() const
{
	return m_EMVector.ud_v2.fCamSpeed;
}

void		CEMVector::SetCameraSense(float fVal)
{
	if(fVal < 0)
		return;

	m_EMVector.ud_v2.fCamSense = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraSense() const
{
	return m_EMVector.ud_v2.fCamSense;
}

void		CEMVector::SetCameraElectronDose(float fVal)
{
	if(fVal < 0)
		return;

	m_EMVector.ud_v2.fCamDose = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraElectronDose() const
{
	return m_EMVector.ud_v2.fCamDose;
}

void		CEMVector::SetCameraGroup(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraGroup() const
{
	return m_EMVector.ud_v2.fCamMisc[0];
}

void		CEMVector::SetGainIndex(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetGainIndex() const
{
	return m_EMVector.ud_v2.fCamMisc[1];
}

void		CEMVector::SetSpeedIndex(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[2] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetSpeedIndex() const
{
	return m_EMVector.ud_v2.fCamMisc[2];
}

void		CEMVector::SetCorrectionNrImages(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[3] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCorrectionNrImages() const
{
	return m_EMVector.ud_v2.fCamMisc[3];
}

void		CEMVector::SetCameraDynamic(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[4] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetCameraDynamic() const
{
	return m_EMVector.ud_v2.fCamMisc[4];
}

void		CEMVector::SetVersion(int lVal)
{
	m_EMVector.version = lVal;
	//InformObservers(false);
}

int		CEMVector::GetVersion() const
{
	return m_EMVector.version;
}

bool		CEMVector::GetImgStatus(_complex *pbuf, int nNum)
{
	if(nNum < 0 || nNum > 16)
		return false;

	for(int i=0; i<nNum; i++)
		pbuf[i] = m_EMVector.ud_v2.cmplxImgStat[i];

	return true;
}

void		CEMVector::SetImgStatMini(double dVal)
{
	m_EMVector.ud_v2.cmplxImgStat[0].x = dVal;
	//InformObservers(false);
}

double		CEMVector::GetImgStatMini() const
{
	return m_EMVector.ud_v2.cmplxImgStat[0].x;
}

void		CEMVector::SetImgStatMax(double dVal)
{
	m_EMVector.ud_v2.cmplxImgStat[1].x = dVal;
	//InformObservers(false);
}

double		CEMVector::GetImgStatMax() const
{
	return m_EMVector.ud_v2.cmplxImgStat[1].x;
}

void		CEMVector::SetImgStatMean(double dVal)
{
	m_EMVector.ud_v2.cmplxImgStat[2].x = dVal;
	//InformObservers(false);
}

double		CEMVector::GetImgStatMean() const
{
	return m_EMVector.ud_v2.cmplxImgStat[2].x;
}

void		CEMVector::SetImgStatDev(double dVal)
{
	m_EMVector.ud_v2.cmplxImgStat[3].x = dVal;
	//InformObservers(false);
}

double		CEMVector::GetImgStatDev() const
{
	return m_EMVector.ud_v2.cmplxImgStat[3].x;
}

void		CEMVector::SetBorder(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[5] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetBorder() const
{
	return m_EMVector.ud_v2.fCamMisc[5];
}

void		CEMVector::SetClippedRadiusPercent(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[11] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetClippedRadiusPercent() const
{
	return m_EMVector.ud_v2.fImgScaling[11];
}

void		CEMVector::SetUseLogarithm(bool bVal)
{
	 m_EMVector.ud_v2.fImgScaling[12] = bVal ? 1.f : 0.f;
	 //InformObservers(false);
}

bool		CEMVector::GetUseLogarithm() const
{
	bool bVal;
	bVal = (m_EMVector.ud_v2.fImgScaling[12] == 1.f);
	return bVal;
}

void		CEMVector::SetTEMabberationsCS(float fVal)
{
	m_EMVector.ud_v2.fTemAberr[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMabberationsCS() const
{
	return m_EMVector.ud_v2.fTemAberr[0];
}

float		CEMVector::GetTEMabberationsCC() const
{
	return m_EMVector.ud_v2.fTemAberr[5];
}

void		CEMVector::SetTEMabberationsCC(float fVal)
{
	m_EMVector.ud_v2.fTemAberr[5] = fVal;
	//InformObservers(false);
}

void		CEMVector::SetTEMenergy(float fVal)
{
	if(fVal < 0)
		return;

	m_EMVector.ud_v2.fTemEnergy[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMenergy() const
{
	return m_EMVector.ud_v2.fTemEnergy[1];
}

float		CEMVector::GetTEMImageShiftX() const
{
	return m_EMVector.ud_v2.fTemImgShift[0];
}

void		CEMVector::SetTEMImageShiftX(float fVal)
{
	m_EMVector.ud_v2.fTemImgShift[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMImageShiftY() const
{
	return m_EMVector.ud_v2.fTemImgShift[1];
}

void		CEMVector::SetTEMImageShiftY(float fVal)
{
	m_EMVector.ud_v2.fTemImgShift[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMbeamshiftX() const
{
	return m_EMVector.ud_v2.fTemBeamShift[0];
}

void		CEMVector::SetTEMbeamshiftX(float fVal)
{
	m_EMVector.ud_v2.fTemBeamShift[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMbeamshiftY() const
{
	return m_EMVector.ud_v2.fTemBeamShift[1];
}

void		CEMVector::SetTEMbeamshiftY(float fVal)
{
	m_EMVector.ud_v2.fTemBeamShift[1] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMbeamtilt() const
{
	return m_EMVector.ud_v2.fTemBeamTilt[0];
}

void		CEMVector::SetTEMbeamtilt(float fVal)
{
	m_EMVector.ud_v2.fTemBeamTilt[0] = fVal;
	//InformObservers(false);
}

void		CEMVector::SetTEMMisc(float fVal)
{
	if(fVal < -360 || fVal > 360)
		return;

	m_EMVector.ud_v2.fTemMisc[0] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetTEMMisc() const
{
	return m_EMVector.ud_v2.fTemMisc[0];
}

string		CEMVector::GetTecnaiMicInfo()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szAdaTietzMicInfo)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szAdaTietzMicInfo);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szAdaTietzMicInfo, pVal);
#endif

	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

string		CEMVector::GetTecnaiSpecInfo()
{
	CStdString strVal(sizeof(m_EMVector.ud_v2.szAdaTietzSpecInfo)/sizeof(TCHAR), '\0');
	TCHAR *pVal = strVal.GetBuffer(strVal.GetLength());

#ifdef _UNICODE
	//strVal = string(m_EMVector.ud_v2.szAdaTietzSpecInfo);
#else
	ConvertFromUnicodeToChar(m_EMVector.ud_v2.szAdaTietzSpecInfo, pVal);
#endif

	strVal.Trim();
	string st(strVal.begin(), strVal.end());
	return st;
}

/*void		CEMVector::SetInformObservers(bool bVal)
{
	m_bInformObservers = bVal;
}*/

/*bool		CEMVector::GetInformObservers() const
{
	return m_bInformObservers;
}*/

/*void		CEMVector::InformObservers(bool bScalebar)
{
	if(!m_bInformObservers)
		return;

	itEMVObs itObs;
	for(itObs = m_vectObservers.begin(); itObs != m_vectObservers.end(); itObs++)
	{
		(*itObs)->Update();
		if(bScalebar)
			(*itObs)->UpdateScalebar();
	}
}*/

string		CEMVector::GetValue(EMVECTOR_PARAMS emvParam)
{
	CStdString strVal;

	int	lVal = 0;
	float	fVal = 0;
	float	fValues[32];

	BOOL bOK = FALSE;

	switch(emvParam)
	{
		case IMAGE_NAME:
			{
				strVal = GetImageName( );
				bOK = TRUE;
			}break;
		case IMAGE_FOLDER:
			{
				strVal = GetImageFolder();
				bOK = TRUE;
			}break;
		case IMAGE_SIZE_X:
			{
				strVal.Format( _T("%ld"), GetImageSizeX() );
				bOK = TRUE;
			}break;
		case IMAGE_SIZE_Y:
			{
				strVal.Format( _T("%ld"), GetImageSizeY() );
				bOK = TRUE;
			}break;
		case IMAGE_SIZE_Z:
			{
				strVal.Format( _T("%ld"), GetImageSizeZ() );
				bOK = TRUE;
			}break;
		case IMAGE_SIZE_E:
			{
				strVal.Format( _T("%ld"), GetImageSizeE() );
				bOK = TRUE;
			}break;
		case IMAGE_DATA_TYPE:
			{
				strVal = GetDataType((int)GetImageDataType());
				bOK = TRUE;
			}break;
		case IMAGE_CREATION_DATE:
			{
				int	lTmp1, lTmp2, lTmp3;

				lTmp1 = GetImageCreationDate() / 65536;
				lTmp2 = GetImageCreationDate() - ( lTmp1 * 65536 );
				lTmp2 /= 256;
				lTmp3 = GetImageCreationDate() - ( lTmp1 * 65536 + lTmp2 * 256);
				strVal.Format(_T("%ld-%02ld-%02ld"), lTmp1, lTmp2, lTmp3 );
				bOK = TRUE;
			}break;
		case IMAGE_CREATION_TIME:
			{
				int	lTmp1, lTmp2, lTmp3;

				lTmp1 = GetImageCreationTime() / 3600;
				lTmp2 = GetImageCreationTime() - ( lTmp1 * 3600 );
				lTmp2 /= 60;
				lTmp3 = GetImageCreationTime() - ( lTmp1 * 3600 + lTmp2 * 60);
				strVal.Format(_T("%02ld:%02ld:%02ld"), lTmp1, lTmp2, lTmp3 );
				bOK = TRUE;
			}break;
		case IMAGE_COMMENT:
			{
				strVal = GetImageComment();
				bOK = TRUE;
			}break;
		case IMAGE_SCALING_LOWER:
			{
				GetScaleModes(&fValues[0], 16);	
				float fLeftVal = fValues[3];
				strVal.Format(_T("%.0f"), fLeftVal);
				bOK = TRUE;
			}break;
		case IMAGE_SCALING_UPPER:
			{
				GetScaleModes(&fValues[0], 16);	
				fVal = fValues[4];
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_MIN:
			{
				fVal = (float)GetImgStatMini();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_MAX:
			{
				fVal = (float)GetImgStatMax();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_MEAN:
			{
				fVal = (float)GetImgStatMean();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_STD_DEVIATION:
			{
				fVal = (float)GetImgStatDev();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_TYPE_PARAM:
			{
				int nType = GetImageType();
				switch( nType )
				{
				case 0:		strVal = _T("(Image)");			break;
				case 1:		strVal = _T("(Dark image)");	break;
				case 2:		strVal = _T("(Gain image)");	break;
				case 3:		strVal = _T("(Power)");			break;
				case 4:		strVal = _T("(Diffraction)");	break;
				default:	strVal = _T("(unknown)");		break;
				}
				bOK = TRUE;
			}break;
		case IMAGE_PIXEL_SIZE_X:
			{
				fVal = GetImageDistX();			
				strVal.Format(_T("%.7f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_PIXEL_SIZE_Y:
			{
				fVal = GetImageDistY();			
				strVal.Format(_T("%.7f"), fVal);
				bOK = TRUE;
			}break;
		case IMAGE_PIXEL_SIZE_Z:
			{
				fVal = GetImageDistZ();			
				strVal.Format(_T("%.3f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_TYPE:
			{
				strVal = GetTEMType();						
				bOK = TRUE;
			}break;
		case TEM_HT:
			{
				fVal = GetTEMHighTension();			
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_ABBERATIONS_CS:
			{
				fVal = GetTEMabberationsCS();	
				strVal.Format(_T("%.1f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_ABBERATIONS_CC:
			{
				fVal = GetTEMabberationsCC();	
				strVal.Format(_T("%.1f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_ENERGY:
			{
				fVal = GetTEMenergy();
				strVal.Format(_T("%.1f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_MAGNIFICATION:
			{
				fVal = GetCurrentMagnificationRelativeToScreen();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_MAG_COR:
			{
				fVal = GetCorrectedMagnification();	
				strVal.Format(_T("%.3f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_POSTMAG:
			{
				fVal = GetPostMagnification();		
				strVal.Format(_T("%.3f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_GONIO_POS_X:
			{
				fVal = GetGonioX();
				strVal.Format(_T("%.3f"), fVal*1e6f);	// MCL works in [m], outout is in [?m]
				bOK = TRUE;
			}break;
		case TEM_GONIO_POS_Y:
			{
				fVal = GetGonioY();
				strVal.Format(_T("%.3f"), fVal*1e6f);	// MCL works in [m], outout is in [?m]
				bOK = TRUE;
			}break;
		case TEM_GONIO_POS_Z:
			{
				fVal = GetGonioZ();
				strVal.Format(_T("%.3f"), fVal*1e6);	// MCL works in [m], outout is in [?m]
				bOK = TRUE;
			}break;
		case TEM_GONIO_TILT_A:
			{
				fVal = GetGonioTiltA();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_GONIO_TILT_B:
			{
				fVal = GetGonioTiltB();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_IMAGE_SHIFT_X:
			{
				fVal = GetTEMImageShiftX();		
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_IMAGE_SHIFT_Y:
			{
				fVal = GetTEMImageShiftY();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_BEAM_SHIFT_X:
			{
				fVal = GetTEMbeamshiftX();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_BEAM_SHIFT_Y:
			{
				fVal = GetTEMbeamshiftY();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_BEAM_TILT:
			{
				fVal = GetTEMbeamtilt();		
				strVal.Format(_T("%.1f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_TILING:
			{
				fVal = (float)GetTEMTilingType();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_TILING_NR_IMG_X:
			{
				fVal = GetTEMTilingNrImgX();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_TILING_NR_IMG_Y:
			{
				fVal = GetTEMTilingNrImgY();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_TILING_OVERLAP_X:
			{
				fVal = GetTEMTilingOverlapX();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_TILING_OVERLAP_Y:
			{
				fVal = GetTEMTilingOverlapY();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_SPOT_SIZE:
			{
				fVal = GetTEMIlluminationSpotSize();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_INTENSITY:
			{
				fVal = GetTEMIlluminationIntensity();
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_SHUTTER_TYPE:
			{
				lVal = GetTEMShutter();				
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_MISC:
			{
				fVal = GetTEMMisc();	
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_TYPE:
			{
				strVal = GetCameraType();
				bOK = TRUE;
			}break;
		case CAMERA_PIXELSIZE_X:
			{
				fVal = GetCameraPixelSizeX();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_PIXELSIZE_Y:
			{
				fVal = GetCameraPixelSizeY();
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_CCD_OFFSET_X:
			{
				lVal = GetCameraOffsetX();			
				strVal.Format( _T("%ld"), lVal );
				bOK = TRUE;
			}break;
		case CAMERA_CCD_OFFSET_Y:
			{
				lVal = GetCameraOffsetY();			
				strVal.Format( _T("%ld"), lVal );
				bOK = TRUE;
			}break;
		case CAMERA_BINNING_X:
			{
				lVal = GetCameraBinningX();			
				strVal.Format( _T("%ld"), lVal );
				bOK = TRUE;
			}break;
		case CAMERA_BINNING_Y:
			{
				lVal = GetCameraBinningY();			
				strVal.Format( _T("%ld"), lVal );
				bOK = TRUE;
			}break;
		case CAMERA_EXP_TIME:
			{
				fVal = GetCameraExposureTime();		
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_GAIN:
			{
				fVal = GetCameraGainFactor();			
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_SPEED:
			{
				fVal = GetCameraReadoutRate();		
				strVal.Format(_T("%.0f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_SENSITIVITY:
			{
				fVal = GetCameraSense();				
				strVal.Format(_T("%.2f"), fVal);
				bOK = TRUE;
			}break;
		case CAMERA_DOSE:
			{
				fVal = GetCameraElectronDose();		
				strVal.Format(_T("%.1f"), fVal);
				bOK = TRUE;
			}break;
		case TEM_MODE:
			{
				int nType = GetTEMMode();
				switch( nType )
				{
				case 0:		strVal = _T("(Brightfield)");		break;
				case 1:		strVal = _T("(DarkField Carth.)");	break;
				case 2:		strVal = _T("(DarkField Conical)");	break;
				case 3:		strVal = _T("(Diffraction)");		break;
				default:	strVal = _T("(unknown)");			break;
				}
				bOK = TRUE;
			}break;
		case CAM_SCX_AMPLIFIER:
			{
				float fVal;
				fVal = GetSCXAmplifier();

				int lVal = (int)(fVal+.1);
				switch( lVal )
				{
				case 0:		strVal = _T("<undefined>");		break;
				case 1:		strVal = _T("(low noise)");		break;
				case 2:		strVal = _T("(high capacity)");	break;
				case 3:		strVal = _T("[undefined]");		break;
				default:	strVal = _T("-");				break;
				}
				bOK = TRUE;
			}
			break;
		case CAM_NUM_PORTS:
			{
				float fVal;
				fVal = GetNumberOfPorts();

				strVal.Format(_T("%ld"), (int)(fVal+.1f));
				bOK = TRUE;
			}
			break;
		case CAM_READOUT_MODE:
			{
				float fVal;
				fVal = GetReadoutMode();

				int lVal = (int)(fVal+.1);
				switch( lVal )
				{
				case 0:
				case 1:		strVal = _T("(normal)");	break;
				case 2:		strVal = _T("(frame transfer)");	break;
				case 3:		strVal = _T("(rolling shutter)");	break;
				default:	strVal = _T("-");		break;
				}
				bOK = TRUE;
			}
			break;
		case CAM_GEOMETRY:
			{
				float fVal;
				fVal = GetReadoutGeometry();

				int lVal = (int)(fVal+.1);

				strVal = _T("(normal)");

				if( lVal > 0 )
				{
					BOOL bFound = FALSE;

					strVal = _T("(");

					if( (lVal &  1) ==  1 )	{ strVal += "mirror h-axis"; bFound = TRUE; }
					if( (lVal &  2) ==  2 )	{ if( bFound ) strVal += _T(" + "); strVal += "mirror v-axis"; bFound = TRUE; }
					if( (lVal &  4) ==  4 )	{ if( bFound ) strVal += _T(" + "); strVal += "rotation 90?";  bFound = TRUE; }
					if( (lVal &  8) ==  8 )	{ if( bFound ) strVal += _T(" + "); strVal += "rotation 180?"; bFound = TRUE; }
					if( (lVal & 16) == 16 )	{ if( bFound ) strVal += _T(" + "); strVal += "rotation -90?"; bFound = TRUE; }
					strVal += _T(")");
				}
				bOK = TRUE;
			}
			break;
		case IMAGE_FF_NUMBER:
			{
				float fVal;
				fVal = GetFlatfieldNumber();

				strVal.Format(_T("%ld"), (int)(fVal+.1f));
				bOK = TRUE;
			}
			break;
	}
	strVal.Trim();

	if(!bOK)
		strVal.empty();
	string st(strVal.begin(), strVal.end());
	return st;
}

BOOL		CEMVector::SetValue(EMVECTOR_PARAMS emvParam, const string &strVal)
{
	int	lVal = 0;
	float	fVal = 0;
	BOOL bOK = FALSE;
	LPTSTR aEndPtr;

	CStdString csValue(strVal);


	switch(emvParam)
	{
		case IMAGE_NAME:
			{
				SetImageName( strVal );
				bOK = TRUE;
			}break;
		case IMAGE_FOLDER:
			{
				string st(csValue.begin(), csValue.end());
				SetImageFolder(st);
				//SetImageFolder(csValue.GetBuffer(csValue.GetLength()));
				bOK = TRUE;
			}break;
		case IMAGE_SIZE_X:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetImageSizeX(lVal);
					bOK = TRUE;
				}
			}break;
		case IMAGE_SIZE_Y:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetImageSizeY(lVal);
					bOK = TRUE;
				}
			}break;
		case IMAGE_SIZE_Z:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetImageSizeZ(lVal);
					bOK = TRUE;
				}
			}break;
		case IMAGE_SIZE_E:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetImageSizeE(lVal);		
					bOK = TRUE;
				}
			}break;
			/*
		case IMAGE_DATA_TYPE:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetImageDataType((IMAGE_TYPE)lVal);
					bOK = TRUE;
				}
			}break;
		case IMAGE_CREATION_DATE:
			{
				CString csValue;
				long	lTmp1, lTmp2, lTmp3;

				lTmp1 = GetImageCreationDate() / 65536;
				lTmp2 = GetImageCreationDate() - ( lTmp1 * 65536 );
				lTmp2 /= 256;
				lTmp3 = GetImageCreationDate() - ( lTmp1 * 65536 + lTmp2 * 256);
				csValue.Format(_T("%ld-%02ld-%02ld"), lTmp1, lTmp2, lTmp3 );
				bAdd ? AddUserDataString( lLine, strDescription, csValue, emvParam) : UpdateUserDataString(	lLine, csValue, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_CREATION_TIME:
			{
				CString csValue;
				long	lTmp1, lTmp2, lTmp3;

				lTmp1 = GetImageCreationTime() / 3600;
				lTmp2 = GetImageCreationTime() - ( lTmp1 * 3600 );
				lTmp2 /= 60;
				lTmp3 = GetImageCreationTime() - ( lTmp1 * 3600 + lTmp2 * 60);
				csValue.Format(_T("%02ld:%02ld:%02ld"), lTmp1, lTmp2, lTmp3 );
				bAdd ? AddUserDataString( lLine, strDescription, csValue, emvParam) : UpdateUserDataString(	lLine, csValue, emvParam );
				bOK = TRUE;
			}break;
			*/
		case IMAGE_COMMENT:
			{
				if(csValue.GetLength() > 512)
					csValue = csValue.Mid(0, 512);
				string st(csValue.begin(), csValue.end());
				SetImageFolder(st);
				//SetImageComment(csValue.GetBuffer(csValue.GetLength()));
				bOK = TRUE;
			}break;
		case IMAGE_SCALING_LOWER:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetLeftScalingValue(fVal);	
					bOK = TRUE;
				}
			}break;
		case IMAGE_SCALING_UPPER:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetRightScalingValue(fVal);
					bOK = TRUE;
				}
			}break;
			/*
		case IMAGE_MIN:
			{
				fVal = (float)GetImgStatMini();
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat( lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_MAX:
			{
				fVal = (float)GetImgStatMax();
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat( lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_MEAN:
			{
				fVal = (float)GetImgStatMean();
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat( lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_STD_DEVIATION:
			{
				fVal = (float)GetImgStatDev();
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat( lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_TYPE_PARAM:
			{
				CString csValue;
				int nType = GetImageType();
				csValue = (nType==0) ? _T("(Image)") : (nType==1) ? _T("(Dark image)") : (nType==2) ? _T("(Gain image)") : _T("(unknwon)");
				bAdd ? AddUserDataString( lLine, strDescription, csValue, emvParam) : UpdateUserDataString(	lLine, csValue, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_PIXEL_SIZE_X:
			{
				fVal = GetImageDistX();			
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat( lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_PIXEL_SIZE_Y:
			{
				fVal = GetImageDistY();			
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat(	lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
		case IMAGE_PIXEL_SIZE_Z:
			{
				fVal = GetImageDistZ();			
				bAdd ? AddUserDataFloat( lLine, strDescription, fVal, emvParam ) : UpdateUserDataFloat( lLine, fVal, emvParam );
				bOK = TRUE;
			}break;
			*/
		case TEM_TYPE:
			{
				string st(csValue.begin(), csValue.end());
				SetImageFolder(st);
				//SetTEMType(csValue.GetBuffer(csValue.GetLength()));						
				bOK = TRUE;
			}break;
		case TEM_HT:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMHighTension(fVal);			
					bOK = TRUE;
				}
			}break;
		case TEM_ABBERATIONS_CS:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMabberationsCS(fVal);	
					bOK = TRUE;
				}
			}break;
		case TEM_ABBERATIONS_CC:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMabberationsCC(fVal);	
					bOK = TRUE;
				}
			}break;
		case TEM_ENERGY:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMenergy(fVal);
					if(fVal >=0 )
						bOK = TRUE;
				}
			}break;
		case TEM_MAGNIFICATION:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCurrentMagnificationRelativeToScreen(fVal);
					if(fVal > 0)
					{
						CalculateImageDist();
						//UpdateInfo(FALSE);
						bOK = TRUE;
					}
				}
			}break;
		case TEM_MAG_COR:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCorrectedMagnification(fVal);
					if(fVal > 0)
					{
						CalculateImageDist();
						bOK = TRUE;
					}
				}
			}break;
		case TEM_POSTMAG:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetPostMagnification(fVal);
					if(fVal > 0)
					{
						CalculateImageDist();
						bOK = TRUE;
					}
				}
			}break;
		case TEM_GONIO_POS_X:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetGonioX(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_GONIO_POS_Y:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetGonioY(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_GONIO_POS_Z:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetGonioZ(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_GONIO_TILT_A:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetGonioTiltA(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_GONIO_TILT_B:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetGonioTiltB(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_IMAGE_SHIFT_X:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMImageShiftX(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_IMAGE_SHIFT_Y:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMImageShiftY(fVal);		
					bOK = TRUE;
				}
			}break;
		case TEM_BEAM_SHIFT_X:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMbeamshiftX(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_BEAM_SHIFT_Y:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMbeamshiftY(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_BEAM_TILT:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMbeamtilt(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_TILING:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMTilingType((int)fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_TILING_NR_IMG_X:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMTilingNrImgX(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_TILING_NR_IMG_Y:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMTilingNrImgY(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_TILING_OVERLAP_X:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMTilingOverlapX(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_TILING_OVERLAP_Y:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMTilingOverlapY(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_SPOT_SIZE:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMIlluminationSpotSize(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_INTENSITY:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMIlluminationIntensity(fVal);
					bOK = TRUE;
				}
			}break;
		case TEM_SHUTTER_TYPE:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMShutter(lVal);				
					bOK = TRUE;
				}
			}break;
		case TEM_MISC:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetTEMMisc(fVal);	
					bOK = TRUE;
				}
			}break;
		case CAMERA_TYPE:
			{
				string st(csValue.begin(), csValue.end());
				SetImageFolder(st);
				//SetCameraType(csValue.GetBuffer(csValue.GetLength()));
				bOK = TRUE;
			}break;
		case CAMERA_PIXELSIZE_X:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraPixelSizeX(fVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_PIXELSIZE_Y:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraPixelSizeY(fVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_CCD_OFFSET_X:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraOffsetX(lVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_CCD_OFFSET_Y:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraOffsetY(lVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_BINNING_X:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraBinningX(lVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_BINNING_Y:
			{
				lVal = _tcstol(csValue, &aEndPtr, 10);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraBinningY(lVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_EXP_TIME:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraExposureTime(fVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_GAIN:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraGainFactor(fVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_SPEED:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraReadoutRate(fVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_SENSITIVITY:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraSense(fVal);
					bOK = TRUE;
				}
			}break;
		case CAMERA_DOSE:
			{
				fVal = (float)_tcstod(csValue, &aEndPtr);
				if (! (*aEndPtr) && errno != ERANGE)
				{
					SetCameraElectronDose(fVal);
					bOK = TRUE;
				}
			}break;
	}

	return bOK;
}

string		CEMVector::GetDataType(int lVal)
{
	CStdString csType;
	IMAGE_TYPE imgType= (IMAGE_TYPE)lVal;

	switch(imgType)
	{
		case DT_UCHAR:				/* unsigned char ( 8 bit )       */
			{
				csType = _T("unsigned char ( 8 bit )");
			} break;
		case DT_USHORT:				/* unsigned short ( 16 bit )     */
			{
				csType = _T("unsigned short ( 16 bit )");
			} break;
		case DT_SHORT:				/* signed short ( 16 bit )       */
			{
				csType = _T("signed short ( 16 bit )");
			} break;
		case DT_LONG:					/* signed long( 32 bit )         */
			{
				csType = _T("signed long( 32 bit )");
			} break;
		case DT_FLOAT:				/* float ( 32 bit )              */
			{
				csType = _T("float ( 32 bit )");
			} break;
		case DT_DOUBLE:				/* double ( 64 bit )             */
			{
				csType = _T("double ( 64 bit )");
			} break;
		case DT_COMPLEX:				/* complex ( 2x32 bit )          */
			{
				csType = _T("complex ( 2x32 bit )");
			} break;
		case DT_STRING:				/* string                        */
			{
				csType = _T("string");
			} break;
		case DT_BINARY:				/* 8bit, but just 0 and 1 are OK */
			{
				csType = _T("8bit");
			} break;
		case DT_RGB8:					/* 3* 8bit unsigned char         */
			{
				csType = _T("3* 8bit unsigned char");
			} break;
		case DT_RGB16:				/* 3* 16bit unsigned char        */
			{
				csType = _T("3* 16bit unsigned char");
			} break;
		default:
			{
				csType = _T("unknown type");
			}
	}
	string csT(csType.begin(), csType.end());
	return csT;
}

void		CEMVector::CalculateImageDist()
{
	float fCamPixelX = 0, fCamPixelY = 0;

	fCamPixelX = GetCameraPixelSizeX();
	fCamPixelY = GetCameraPixelSizeY();

	fCamPixelX *= GetCameraBinningX();	
	fCamPixelY *= GetCameraBinningY();	

	float fMagScr = GetCurrentMagnificationRelativeToScreen();
	float fMagCor = GetCorrectedMagnification();	
	float fMagPost = GetPostMagnification();		

	//_ASSERTE(fMagScr != 0 && fMagCor != 0 && fMagPost != 0);

	fCamPixelX /= (fMagScr * fMagCor * fMagPost);
	fCamPixelY /= (fMagScr * fMagCor * fMagPost);
	
	SetImageDistX(fCamPixelX);
	SetImageDistY(fCamPixelY);
}

void		CEMVector::SetSCXAmplifier(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[6] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetSCXAmplifier() const
{
	return m_EMVector.ud_v2.fCamMisc[6];
}

void		CEMVector::SetNumberOfPorts(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[7] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetNumberOfPorts() const
{
	return m_EMVector.ud_v2.fCamMisc[7];
}

void		CEMVector::SetReadoutMode(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[8] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetReadoutMode() const
{
	return m_EMVector.ud_v2.fCamMisc[8];
}

void		CEMVector::SetReadoutGeometry(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[9] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetReadoutGeometry() const
{
	return m_EMVector.ud_v2.fCamMisc[9];
}

void		CEMVector::SetFlatfieldNumber(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[10] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetFlatfieldNumber() const
{
	return m_EMVector.ud_v2.fCamMisc[10];
}

void		CEMVector::SetLSDiffractionCenterX(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[11] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetLSDiffractionCenterX() const
{
	return m_EMVector.ud_v2.fCamMisc[11];
}

void		CEMVector::SetLSDiffractionCenterY(float fVal)
{
	m_EMVector.ud_v2.fCamMisc[12] = fVal;
	//InformObservers(false);
}
float		CEMVector::GetLSDiffractionCenterY() const
{
	return m_EMVector.ud_v2.fCamMisc[12];
}

//-------------------------------

void		CEMVector::SetHistogramRangeMode(float fVal)
{
	if( ( fVal < 1 ) || ( fVal > 3 ) )
		return;

	m_EMVector.ud_v2.fImgScaling[13] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetHistogramRangeMode() const
{
	return m_EMVector.ud_v2.fImgScaling[13];
}

void		CEMVector::SetHistogramRangeLow(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[14] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetHistogramRangeLow() const
{
	return m_EMVector.ud_v2.fImgScaling[14];
}

void		CEMVector::SetHistogramRangeHigh(float fVal)
{
	m_EMVector.ud_v2.fImgScaling[15] = fVal;
	//InformObservers(false);
}

float		CEMVector::GetHistogramRangeHigh() const
{
	return m_EMVector.ud_v2.fImgScaling[15];
}

