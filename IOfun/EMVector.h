//#pragma once

#ifndef _DATATYPES
	#include "DataTypes.h"
	#define _DATATYPES
#endif
#include "EMVectorDefines.h"
#include "string.h"
#include "StdString.h"
#include "stdio.h"
#include "math.h"

using namespace std;

class CEMVector
{
public:
	CEMVector(); //ctr
	CEMVector(const EMVECTOR *pEMV); //2nd ctr
	virtual ~CEMVector();

	CEMVector(CEMVector& other); //cctor
	CEMVector& operator=(CEMVector& src);
	CEMVector& operator=(const EMVECTOR *pEMV); 

	bool operator==(const CEMVector &other);

	//GetStructEMVECTOR() should be used only for get,put EMVECTOR to COM methods
	//const EMVECTOR& GetStructEMVECTOR();
	EMVECTOR&	GetStructEMVECTOR();

	bool	IsEmpty();

	string		GetValue(EMVECTOR_PARAMS emvParam);
	BOOL		SetValue(EMVECTOR_PARAMS emvParam, const string &strVal);

	////properties //////////////////////////////////////////////////////////////////////

	void		SetImageDataType(IMAGE_TYPE imgType);
	IMAGE_TYPE	GetImageDataType() const;
	
	PROPERTY_STRING( ImageName )
	PROPERTY_STRING(ImageFolder)
	PROPERTY_LONG(Version)
	PROPERTY_LONG(ImageSizeX)
	PROPERTY_LONG(ImageSizeY)
	PROPERTY_LONG(ImageSizeZ)
	PROPERTY_LONG(ImageSizeE)
	PROPERTY_LONG(ImageCreationDate)
	PROPERTY_LONG(ImageCreationTime)
	PROPERTY_STRING(ImageComment)
	PROPERTY_STRING(ImageHistory)
	PROPERTY_FLOAT(Contrast)
	PROPERTY_FLOAT(Brightness)
	PROPERTY_FLOAT(Overflow)
	PROPERTY_FLOAT(Underflow)
	PROPERTY_FLOAT(MinScalingValue)
	PROPERTY_FLOAT(MaxScalingValue)
	PROPERTY_FLOAT(LeftScalingValue)
	PROPERTY_FLOAT(RightScalingValue)
	PROPERTY_FLOAT(LeftPowerScalingValue)
	PROPERTY_FLOAT(RightPowerScalingValue)
	PROPERTY_FLOAT(ScalingAbsoluteMinimum)
	PROPERTY_FLOAT(ScalingMinimum)
	PROPERTY_FLOAT(ScalingMaximum)
	PROPERTY_FLOAT(ScalingAbsoluteMaximum)
	PROPERTY_INT(ImageType)
	PROPERTY_LONG(DisplayType)
	PROPERTY_FLOAT(ImageDistX)
	PROPERTY_FLOAT(ImageDistY)
	PROPERTY_FLOAT(ImageDistZ)
	PROPERTY_FLOAT(ImageDistE)
	PROPERTY_FLOAT(NavigatorPositionX)
	PROPERTY_FLOAT(NavigatorPositionY)
	PROPERTY_FLOAT(CameraPixelSizeX)
	PROPERTY_FLOAT(CameraPixelSizeY)
	PROPERTY_LONG(CameraOffsetX)
	PROPERTY_LONG(CameraOffsetY)
	PROPERTY_LONG(CameraBinningX)
	PROPERTY_LONG(CameraBinningY)
	PROPERTY_FLOAT(CameraExposureTime)
	PROPERTY_FLOAT(CameraGainFactor)
	PROPERTY_FLOAT(CameraReadoutRate)
	PROPERTY_FLOAT(CameraSense)
	PROPERTY_FLOAT(CameraElectronDose)
	PROPERTY_FLOAT(CameraDynamic)
	PROPERTY_DOUBLE(ImgStatMini)
	PROPERTY_DOUBLE(ImgStatMax)
	PROPERTY_DOUBLE(ImgStatMean)
	PROPERTY_DOUBLE(ImgStatDev)
	PROPERTY_FLOAT(Border)
	PROPERTY_FLOAT(ClippedRadiusPercent)
	PROPERTY_BOOL(UseLogarithm)
	PROPERTY_FLOAT(ScaleMode) // Scale mode 0 = (contrast/brightness)
	PROPERTY_FLOAT(HistogramRangeMode) // 0=undefined, 1=dynamic, 2=full, 3=user
	PROPERTY_FLOAT(HistogramRangeLow)
	PROPERTY_FLOAT(HistogramRangeHigh)
	PROPERTY_STRING(CameraType)
	PROPERTY_FLOAT(CameraGroup)
	PROPERTY_FLOAT(GainIndex)
	PROPERTY_FLOAT(SpeedIndex)
	PROPERTY_FLOAT(CorrectionNrImages)
	PROPERTY_STRING(CameraFlat)
	PROPERTY_FLOAT(SCXAmplifier)
	PROPERTY_FLOAT(NumberOfPorts)
	PROPERTY_FLOAT(ReadoutMode)
	PROPERTY_FLOAT(ReadoutGeometry)
	PROPERTY_FLOAT(FlatfieldNumber)
	PROPERTY_FLOAT(LSDiffractionCenterX)
	PROPERTY_FLOAT(LSDiffractionCenterY)

	bool		GetScaleModes(float *buf, int nNums) const;
	//bool		GetCamMisc(float *pBuf, int nNums);

	bool		GetImgStatus(_complex *pbuf, int nNum);

	////TEM
	PROPERTY_STRING(TEMType)
	PROPERTY_FLOAT(TEMHighTension)
	PROPERTY_LONG(TEMMode)
	PROPERTY_FLOAT(CurrentMagnificationRelativeToScreen)
	PROPERTY_FLOAT(CorrectedMagnification)
	PROPERTY_FLOAT(PostMagnification)
	PROPERTY_LONG(StageType)
	PROPERTY_FLOAT(GonioX)	//fTemStgPos[0]
	PROPERTY_FLOAT(GonioY)	//fTemStgPos[1]
	PROPERTY_FLOAT(GonioZ)	//fTemStgPos[2]
	PROPERTY_FLOAT(GonioTiltA)
	PROPERTY_FLOAT(GonioTiltB)
	PROPERTY_INT(TEMTilingType)
	PROPERTY_FLOAT(TEMTilingMaxNumberX)
	PROPERTY_FLOAT(TEMTilingMaxNumberY)
	PROPERTY_FLOAT(TEMIlluminationSpotSize)
	PROPERTY_FLOAT(TEMIlluminationIntensity)
	PROPERTY_LONG(TEMShutter)
	PROPERTY_FLOAT(TEMabberationsCS)
	PROPERTY_FLOAT(TEMabberationsCC)
	PROPERTY_FLOAT(TEMenergy)
	PROPERTY_FLOAT(TEMImageShiftX)
	PROPERTY_FLOAT(TEMImageShiftY)
	PROPERTY_FLOAT(TEMbeamshiftX)
	PROPERTY_FLOAT(TEMbeamshiftY)
	PROPERTY_FLOAT(TEMbeamtilt)
	PROPERTY_FLOAT(TEMTilingNrImgX)
	PROPERTY_FLOAT(TEMTilingNrImgY)
	PROPERTY_FLOAT(TEMTilingOverlapX)
	PROPERTY_FLOAT(TEMTilingOverlapY)
	PROPERTY_FLOAT(TEMMisc)

	string		GetTecnaiMicInfo();
	string		GetTecnaiSpecInfo();
	bool		GetStgPos(float *pBuf, int nNums);
	
	/////

	///////////////////////////////////////////////////////////////////////////////////

private:
	CEMVector& Copy(const EMVECTOR *pEMV);

	void Init();
	void ConvertFromCharToUnicode(const char *pSrc, unsigned short *pDest);
	void ConvertFromUnicodeToChar(const unsigned short *pSrc, char *pDest);

	std::string GetDataType(int lVal);
	void	CalculateImageDist();


private:
	EMVECTOR m_EMVector;
	
};
