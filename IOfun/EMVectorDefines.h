#if !defined(EM_VECT)
#define EM_VECT

#define _T(x) x

#define PROPERTY_BOOL( _name_ )			\
	void Set##_name_ (bool bVal);		\
	bool Get##_name_ () const;

#define PROPERTY_INT( _name_ )			\
	void Set##_name_ (int nVal);		\
	int Get##_name_ () const;
	
#define PROPERTY_LONG( _name_ )			\
	void Set##_name_ (int lVal);		\
	int Get##_name_ () const;

#define PROPERTY_FLOAT( _name_ )		\
	void Set##_name_ (float fVal);		\
	float Get##_name_ () const;

#define PROPERTY_DOUBLE( _name_ ) 		\
	void Set##_name_ (double dVal);		\
	double Get##_name_ () const;

#define PROPERTY_STRING( _name_ )		\
	void Set##_name_ (const string &sstrVal);	\
	string Get##_name_ ();


typedef enum EMVECTOR_PARAMS		//Extended			Editable
{
	UNDEFINED_PARAM		= -1,
	IMAGE_NAME			= 0,		//						x
	IMAGE_FOLDER		= 1,		//						x
	IMAGE_SIZE_X		= 2,		//
	IMAGE_SIZE_Y		= 3,		//
	IMAGE_SIZE_Z		= 4,		//	x
	IMAGE_SIZE_E		= 5,		//	x
	IMAGE_DATA_TYPE		= 6,		//	x
	IMAGE_CREATION_DATE	= 7,		//
	IMAGE_CREATION_TIME	= 8,		//
	IMAGE_COMMENT		= 9,		//						x
	IMAGE_SCALING_LOWER	= 10,		//						x
	IMAGE_SCALING_UPPER	= 11,		//						x
	IMAGE_MIN			= 12,		//
	IMAGE_MAX			= 13,		//
	IMAGE_MEAN			= 14,		//
	IMAGE_STD_DEVIATION	= 15,		//
	IMAGE_TYPE_PARAM	= 16,		//	x
	IMAGE_PIXEL_SIZE_X	= 17,		//
	IMAGE_PIXEL_SIZE_Y	= 18,		//
	IMAGE_PIXEL_SIZE_Z	= 19,		//
	TEM_TYPE			= 20,		//	x					x
	TEM_HT				= 21,		//						x
	TEM_ABBERATIONS_CS	= 22,		//	x					x
	TEM_ABBERATIONS_CC	= 23,		//	x					x
	TEM_ENERGY			= 24,		//	x					x
	TEM_MAGNIFICATION	= 25,		//						x
	TEM_MAG_COR			= 26,		//						x
	TEM_POSTMAG			= 27,		//						x
	TEM_GONIO_POS_X		= 28,		//	x					x
	TEM_GONIO_POS_Y		= 29,		//	x					x
	TEM_GONIO_POS_Z		= 30,		//	x					x
	TEM_GONIO_TILT_A	= 31,		//						x
	TEM_GONIO_TILT_B	= 32,		//	x					x
	TEM_IMAGE_SHIFT_X	= 33,		//	x					x
	TEM_IMAGE_SHIFT_Y	= 34,		//	x					x
	TEM_BEAM_SHIFT_X	= 35,		//	x					x
	TEM_BEAM_SHIFT_Y	= 36,		//	x					x
	TEM_BEAM_TILT		= 37,		//	x					x
	TEM_TILING			= 38,		//	x					x
	TEM_TILING_NR_IMG_X	= 39,		//	x					x
	TEM_TILING_NR_IMG_Y	= 40,		//	x					x
	TEM_TILING_OVERLAP_X= 41,		//	x					x
	TEM_TILING_OVERLAP_Y= 42,		//	x					x
	TEM_SPOT_SIZE		= 43,		//	x					x
	TEM_INTENSITY		= 44,		//	x					x
	TEM_SHUTTER_TYPE	= 45,		//	x
	TEM_MISC			= 46,		//						x
	CAMERA_TYPE			= 47,		//
	CAMERA_PIXELSIZE_X	= 48,		//
	CAMERA_PIXELSIZE_Y	= 49,		//
	CAMERA_CCD_OFFSET_X	= 50,		//
	CAMERA_CCD_OFFSET_Y	= 51,		//
	CAMERA_BINNING_X	= 52,		//
	CAMERA_BINNING_Y	= 53,		//
	CAMERA_EXP_TIME		= 54,		//
	CAMERA_GAIN			= 55,		//
	CAMERA_SPEED		= 56,		//	x
	CAMERA_SENSITIVITY	= 57,		//						x
	CAMERA_DOSE			= 58,		//						x
	TEM_MODE			= 59,		//
	CAM_SCX_AMPLIFIER	= 60,		//
	CAM_NUM_PORTS		= 61,		//
	CAM_READOUT_MODE	= 62,		//
	CAM_GEOMETRY		= 63,		//
	IMAGE_FF_NUMBER		= 64		//
};

#include <math.h>																// CLU: needed for data type _complex
#define EV_VERSION	0x00000002													// CLU: version number of EMVECTOR

struct _complex {
	double x,y;
};

/*
 * CLU:
 * typedefs for EMVECTOR, which collects the historicaly grown, so called "user data" of
 * TEM data format. Till today (17.10.2002) there are two versions of user data, represented
 * by USER_DATA_V1 and USER_DATA_V2.
 */
typedef struct tagUSER_DATA_V1 {
	char	strComment[80];				// CLU: comment
	long	lHighTension;				// CLU: high tension [kV]
	long	lSphericalAberration;		// CLU: spherical abberation [mm]
	long	lIllumAperture;				// CLU: illum. aperture [mrad]
	long	lElOptMag;					// CLU: electron optical magnification [x 1000]
	long	lPostMag;					// CLU: post magnification [x 1]
	long	lFocalLength;				// CLU: focal length [mm]
	long	lDefocus;					// CLU: defocus [nm]
	long	lAstigmatismNM;				// CLU: astigmatism [nm]
	long	lAstigmatismMRAD;			// CLU: astigmatism [mrad]
	long	lBiprismTension;			// CLU: voltage of Biprism [V]
	long	lSpecimenTiltAngle;			// CLU: specimen tilt angle [mrad]
	long	lSpecimenTiltDirection;		// CLU: specimen tilt direction [mrad]
	long	lIllumTiltDirection;		// CLU: illum. tilt direction [mrad]
	long	lIllumTiltAngle;			// CLU: illum. tilt angle [mrad]
	long	lMode;						// CLU: Mode: 0  image, 1  diffraction
	long	lEnergySpread;				// CLU: energy spread [eV]
	long	lChromaticalAberration;		// CLU: chromatically abberation [mm]
	long	lShutterType;				// CLU: 0 = beamblanker, 1 = shutter (or vice versa ?)
	long	lDefocusSpread;				// CLU: spread of defocus [nm]
	long	lCCDNumber;					// CLU: 0 = CCD 0, 10000 = CCD 1, ...
	long	lCCDPixelXY;				// CLU: CCD pixels in x and y direction
	long	lCCDOffsetX;				// CLU: CCD offset in x direction
	long	lCCDOffsetY;				// CLU: CCD offset in y direction
	long	lCCDSinglePixelSize;		// CLU: CCD single pixel size [m m]
	long	lCCDBinningFactor;			// CLU: CCD binning factor
	long	lCCDReadOutSpeed;			// CLU: CCD readout speed [kHz]
	long	lCCDGain;					// CLU: CCD gain ( 0  low, 1  high )
	long	lCCDSensitivity;			// CLU: CCD sensitivity [ADU/primary electron]
	long	lCCDExposureTime;			// CLU: CCD exposure time [ms]
	long	lFlatfieldCorrection;		// CLU: flatfield correction done ( 0  no, 1  yes )
	long	lDeadPixelCorrection;		// CLU: dead pixel correction done ( 0  no, 1  yes )
	long	lMeanValue;					// CLU: mean value [ADU]
	long	lStandardDeviation;			// CLU: standard deviation [ADU]
	long	lDisplacementX;				// CLU: x displacement within series [pixel]
	long	lDisplacementY;				// CLU: y displacement within series [pixel]
	long	lDate;						// CLU: days since 1970-01-01
	long	lTime;						// CLU:	1/100 seconds since midnight
	long	lMinimum;					// CLU: Minimum [ADU]
	long	lMaximum;					// CLU: Maximum [ADU]
	long	lQualityFactor;				// CLU: Quality factor for statistic calculation [%]
} USER_DATA_V1;

typedef struct tagUSER_DATA_V2 {
	unsigned short		szImgName[80];			// CLU: image name
	unsigned short		szImgFolder[80];		// CLU: path to image
	long				lImgSizeX;				// CLU: size x
	long				lImgSizeY;				// CLU: size y
	long				lImgSizeZ;				// CLU: size z
	long				lImgSizeE;				// CLU: size e
	long				lImgDataType;			// CLU: data type
	long				lImgCreationDate;		// CLU: creation date
	long				lImgCreationTime;		// CLU: creation time
	unsigned short		szImgComment[512];		// CLU: image comment
	unsigned short		szImgHistory[512];		// CLU: image history
	float				fImgScaling[16];		// CLU: scaling
	_complex			cmplxImgStat[16];		// CLU: statistics
	long				lImgType;				// CLU: image type
	long				lImgDisplay;			// CLU: display type
	float				fImgDistX;				// CLU: size x
	float				fImgDistY;				// CLU: size y
	float				fImgDistZ;				// CLU: size z
	float				fImgDistE;				// CLU: size e
	float				fImgMisc[32];			// CLU: reserved
	unsigned short		szTemType[80];			// CLU: TEM type
	float				fTemHT;					// CLU: high tension
	float				fTemAberr[32];			// CLU: aberrations, defocus, Cs, Cc, ...
	float				fTemEnergy[32];			// CLU: TEM energy parameters
	long				lTemMode;				// CLU: TEM mode
	float				fTemMagScr;				// CLU: current magnification relative to screen
	float				fTemMagCor;				// CLU: corrected magnification
	float				fTemMagPst;				// CLU: post-magnification
	long				lTemStgType;			// CLU: stage type
	float				fTemStgPos[5];			// CLU: stage position
	float				fTemImgShift[2];		// CLU: image shift
	float				fTemBeamShift[2];		// CLU: beam shift
	float				fTemBeamTilt[2];		// CLU: beam tilt
	float				fTemTiling[7];			// CLU: series:	[0]:	type (0-none, 1-tiling)
												// CLU: 		[1]:	actual number x
												// CLU: 		[2]:	actual number y
												// CLU: 		[3]:	max number x
												// CLU: 		[4]:	max number y
												// CLU: 		[5]:	overlap x
												// CLU: 		[6]:	overlap y
	float				fTemIllum[3];			// CLU: TEM illumination: 0 spotsize
												// CLU: TEM illumination: 1 intensity
	long				lTemShutter;			// CLU: 
	float				fTemMisc[32];			// CLU: reserved
	unsigned short		szCamType[80];			// CLU: camera type
	float				fCamPixel[2];			// CLU: pixelsize
	long				lCamOffX;				// CLU: offset x
	long				lCamOffY;				// CLU: offset y
	long				lCamBinX;				// CLU: binning x
	long				lCamBinY;				// CLU: binning y
	float				fCamExpTime;			// CLU: exposure time
	float				fCamGain;				// CLU: gain factor
	float				fCamSpeed;				// CLU: readout rate
	unsigned short		szCamFlat[80];			// CLU: 
	float				fCamSense;				// CLU: counts per primary electron
	float				fCamDose;				// CLU: electron dose
	float				fCamMisc[32];			// CLU: reserved
	unsigned short		szAdaTietzMicInfo[512];	// CLU: microscope info from FEI'a ADATietz, clipped to 512 characters
	unsigned short		szAdaTietzSpecInfo[512];// CLU: specimen info from FEI'a ADATietz, clipped to 512 characters
} USER_DATA_V2;

typedef struct tagEMVECTOR {
	unsigned long		version;
	USER_DATA_V1		ud_v1;
	USER_DATA_V2		ud_v2;
} EMVECTOR;

#endif
