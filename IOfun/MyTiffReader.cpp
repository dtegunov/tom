// MyTiffReader.cpp : Defines the entry point for the console application.
//
#include "mex.h"
#include "matrix.h"
#ifndef _DATATYPES
	#include "DataTypes.h"
	#define _DATATYPES
#endif
#ifndef _EMVECT
	#include "EMVector.h"
	#define _EMVECT
#endif
#include <iostream>
#include "FileTool.h"

#ifndef _TEMPLATE
	#include "var_list.cpp"
	#define _TEMPLATE
#endif



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
{
    
    char            *Filename;
    double          *outImage;
    Filename                    = mxArrayToString(prhs[0]);
    //Filename                    = "../../Distortion/Data/8kMessung/Shear8kV_001.tif";
    
	BYTE			*pDst		= (BYTE*)mxMalloc(sizeof(BYTE));
	BYTE			*pThumb		= (BYTE*)mxMalloc(128*128*3*sizeof(BYTE));
	EMVECTOR		*EMV		= new EMVECTOR();
	
	
	var_list<BYTE>		*lpDst		= new var_list<BYTE>(pDst, 0, 0, 4);
	var_list<BYTE>		*lpThumb	= new var_list<BYTE>(pThumb, 0, 0, 2);
	var_list<EMVECTOR> 	*lpEMVector 	= new var_list<EMVECTOR>(EMV, 0, 0, 0);

	CFileTool		*CFT		= new CFileTool();

	CFT->DoReadTiff(Filename , lpDst, lpEMVector, lpThumb);                 // Read Image from file 
	CEMVector		*pEMVector	= new CEMVector(EMV);

    plhs[0]       = mxCreateDoubleMatrix(lpDst->lDim[1], lpDst->lDim[0], mxREAL);
    outImage      = mxGetPr(plhs[0]);
    
    
    //******************************** EM Header Information ********************************************
    int dims[2] = {1, 1};
    mxArray *struct_array_ptr, *field_value;
    double *pr;
    
	// 77 elements
    const char *field_names[] = {"ImageName", "ImageFolder", "ImageSizeX", "ImageSizeY", "ImageSizeZ", "ImageSizeE", "DataType",
                                 "CreationDate", "CreationTime", "ImageComment", "ImageHistory", "ImageScaling", "StatMin","StatMax", "StatMean", "StatDevi",
                                 "ImageType", "ImageDisplyType", "ImageDistanceX", "ImageDistanceY", "ImageDistanceZ", "ImageDistanceE",
                                 "TEMType", "HighTension", "AberrationCs","AberrationCc", "EnergyParameters", "TEMMode", "CurMagRelPlate", "CorrMag",
                                 "PostMag", "StageType", "GonioX", "GonioY", "GonioZ", "GonioTiltA", "GonioTiltB", "ImageShiftX",
                                 "ImageShiftY", "BeamShiftX", "BeamShiftY", "BeamTilt", "TilingType","TilingNrImgX","TilingNrImgY", "TilingOverlapX", "TilingOverlapY",
                                 "TEMSpotSize", "TEMIntensity", "TEMShutterType", "CameraType", "CameraPixelSizeX", "CameraPixelSizeY", 
                                 "CameraOffsetX", "CameraOffsetY", "CameraBinningX", "CameraBinningY",
                                 "CameraExposureTime", "CameraGain", "CameraReadoutRate", "CameraSensitivity", "CameraDose", "CamGroup", "CamGainIndex", "CamSpeedIndex",
                                 "CamCorrectionNrImages", "CamDynamic", "CamBorderSize", "CamSCXAmplifier", "CamNumPorts", "CamReadoutMode", "CamReadoutGeometry", 
                                 "CamFlatfieldNumber", "CamLSDiffractionCenterX", "CamLSDiffractionCenterY", "MicroscopeInfo", "SpecimenInfo"};
        
    struct_array_ptr       = mxCreateStructArray(2, dims, 77, field_names);
	
    
    
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_NAME).c_str());
    mxSetField(struct_array_ptr, 0, "ImageName", field_value); 
    
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_FOLDER).c_str());
    mxSetField(struct_array_ptr, 0, "ImageFolder", field_value); 
    
    mxSetField(struct_array_ptr, 0, "ImageSizeX", mxCreateDoubleScalar(pEMVector->GetImageSizeX()));
    
    mxSetField(struct_array_ptr, 0, "ImageSizeY", mxCreateDoubleScalar(pEMVector->GetImageSizeY()));
       
    mxSetField(struct_array_ptr, 0, "ImageSizeZ", mxCreateDoubleScalar(pEMVector->GetImageSizeZ()));
    
    mxSetField(struct_array_ptr, 0, "ImageSizeE", mxCreateDoubleScalar(pEMVector->GetImageSizeE()));
      
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_DATA_TYPE).c_str());
    mxSetField(struct_array_ptr, 0, "DataType", field_value);
    
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_CREATION_DATE).c_str());
    mxSetField(struct_array_ptr, 0, "CreationDate", field_value); 
    
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_CREATION_TIME).c_str());
    mxSetField(struct_array_ptr, 0, "CreationTime", field_value); 
    
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_COMMENT).c_str());
    mxSetField(struct_array_ptr, 0, "ImageComment", field_value); 
    
    field_value = mxCreateString(pEMVector->GetImageHistory().c_str());
    mxSetField(struct_array_ptr, 0, "ImageHistory", field_value); 
    
    field_value = mxCreateDoubleMatrix(16, 1, mxREAL);
    pr = mxGetPr(field_value);
    int i;
    for (i = 0; i < 16; i++)
        pr[i] = EMV->ud_v2.fImgScaling[i];
    mxSetField(struct_array_ptr, 0, "ImageScaling", field_value);
    
    mxSetField(struct_array_ptr, 0, "StatMin", mxCreateDoubleScalar(pEMVector->GetImgStatMini()));
    
    mxSetField(struct_array_ptr, 0, "StatMax", mxCreateDoubleScalar(pEMVector->GetImgStatMax()));
    
    mxSetField(struct_array_ptr, 0, "StatMean", mxCreateDoubleScalar(pEMVector->GetImgStatMean()));
    
    mxSetField(struct_array_ptr, 0, "StatDevi", mxCreateDoubleScalar(pEMVector->GetImgStatDev()));
    
    field_value = mxCreateString(pEMVector->GetValue(IMAGE_TYPE_PARAM).c_str());
    mxSetField(struct_array_ptr, 0, "ImageType", field_value);
    
    mxSetField(struct_array_ptr, 0, "ImageDisplyType", mxCreateDoubleScalar(pEMVector->GetDisplayType()));
    
    mxSetField(struct_array_ptr, 0, "ImageDistanceX", mxCreateDoubleScalar(pEMVector->GetImageDistX()));
    
    mxSetField(struct_array_ptr, 0, "ImageDistanceY", mxCreateDoubleScalar(pEMVector->GetImageDistY()));
    
    mxSetField(struct_array_ptr, 0, "ImageDistanceZ", mxCreateDoubleScalar(pEMVector->GetImageDistZ()));
    
    mxSetField(struct_array_ptr, 0, "ImageDistanceE", mxCreateDoubleScalar(pEMVector->GetImageDistE()));
    
    field_value = mxCreateString(pEMVector->GetValue(TEM_TYPE).c_str());
    mxSetField(struct_array_ptr, 0, "TEMType", field_value); 
    
    mxSetField(struct_array_ptr, 0, "HighTension", mxCreateDoubleScalar(pEMVector->GetTEMHighTension()));
    
    mxSetField(struct_array_ptr, 0, "AberrationCs", mxCreateDoubleScalar(pEMVector->GetTEMabberationsCS()));
    
    mxSetField(struct_array_ptr, 0, "AberrationCc", mxCreateDoubleScalar(pEMVector->GetTEMabberationsCC()));
    
    mxSetField(struct_array_ptr, 0, "EnergyParameters", mxCreateDoubleScalar(pEMVector->GetTEMenergy()));
    
    field_value = mxCreateString(pEMVector->GetValue(TEM_MODE).c_str());
    mxSetField(struct_array_ptr, 0, "TEMMode", field_value);
    
    mxSetField(struct_array_ptr, 0, "CurMagRelPlate", mxCreateDoubleScalar(pEMVector->GetCurrentMagnificationRelativeToScreen()));
    
    mxSetField(struct_array_ptr, 0, "CorrMag", mxCreateDoubleScalar(pEMVector->GetCorrectedMagnification()));
    
    mxSetField(struct_array_ptr, 0, "PostMag", mxCreateDoubleScalar(pEMVector->GetPostMagnification()));
    
    mxSetField(struct_array_ptr, 0, "StageType", mxCreateDoubleScalar(pEMVector->GetStageType()));
    
    mxSetField(struct_array_ptr, 0, "GonioX", mxCreateDoubleScalar(pEMVector->GetGonioX()));
    
    mxSetField(struct_array_ptr, 0, "GonioY", mxCreateDoubleScalar(pEMVector->GetGonioY()));
    
    mxSetField(struct_array_ptr, 0, "GonioZ", mxCreateDoubleScalar(pEMVector->GetGonioZ()));
    
    mxSetField(struct_array_ptr, 0, "GonioTiltA", mxCreateDoubleScalar(pEMVector->GetGonioTiltA()));
    
    mxSetField(struct_array_ptr, 0, "GonioTiltB", mxCreateDoubleScalar(pEMVector->GetGonioTiltB()));
    
    mxSetField(struct_array_ptr, 0, "ImageShiftX", mxCreateDoubleScalar(pEMVector->GetTEMImageShiftX()));
    
    mxSetField(struct_array_ptr, 0, "ImageShiftY", mxCreateDoubleScalar(pEMVector->GetTEMImageShiftY()));
    
    mxSetField(struct_array_ptr, 0, "BeamShiftX", mxCreateDoubleScalar(pEMVector->GetTEMbeamshiftX()));
    
    mxSetField(struct_array_ptr, 0, "BeamShiftY", mxCreateDoubleScalar(pEMVector->GetTEMbeamshiftY()));
    
    mxSetField(struct_array_ptr, 0, "BeamTilt", mxCreateDoubleScalar(pEMVector->GetTEMbeamtilt()));
    
    mxSetField(struct_array_ptr, 0, "TilingType", mxCreateDoubleScalar(pEMVector->GetTEMTilingType()));
    
    mxSetField(struct_array_ptr, 0, "TilingNrImgX", mxCreateDoubleScalar(pEMVector->GetTEMTilingNrImgX()));
    
    mxSetField(struct_array_ptr, 0, "TilingNrImgY", mxCreateDoubleScalar(pEMVector->GetTEMTilingNrImgY()));
    
    mxSetField(struct_array_ptr, 0, "TilingOverlapX", mxCreateDoubleScalar(pEMVector->GetTEMTilingOverlapX()));
    
    mxSetField(struct_array_ptr, 0, "TilingOverlapY", mxCreateDoubleScalar(pEMVector->GetTEMTilingOverlapY()));
    
    mxSetField(struct_array_ptr, 0, "TEMSpotSize", mxCreateDoubleScalar(pEMVector->GetTEMIlluminationSpotSize()));
    
    mxSetField(struct_array_ptr, 0, "TEMIntensity", mxCreateDoubleScalar(pEMVector->GetTEMIlluminationIntensity()));
    
    mxSetField(struct_array_ptr, 0, "TEMShutterType", mxCreateDoubleScalar(pEMVector->GetTEMShutter()));
    
    field_value = mxCreateString(pEMVector->GetValue(CAMERA_TYPE).c_str());
    mxSetField(struct_array_ptr, 0, "CameraType", field_value); 
    
    mxSetField(struct_array_ptr, 0, "CameraPixelSizeX", mxCreateDoubleScalar(pEMVector->GetCameraPixelSizeX()));
    
    mxSetField(struct_array_ptr, 0, "CameraPixelSizeY", mxCreateDoubleScalar(pEMVector->GetCameraPixelSizeY()));
    
    mxSetField(struct_array_ptr, 0, "CameraOffsetX", mxCreateDoubleScalar(pEMVector->GetCameraOffsetX()));
    
    mxSetField(struct_array_ptr, 0, "CameraOffsetY", mxCreateDoubleScalar(pEMVector->GetCameraOffsetY()));
    
    mxSetField(struct_array_ptr, 0, "CameraBinningX", mxCreateDoubleScalar(pEMVector->GetCameraBinningX()));
    
    mxSetField(struct_array_ptr, 0, "CameraBinningY", mxCreateDoubleScalar(pEMVector->GetCameraBinningY()));
    
    mxSetField(struct_array_ptr, 0, "CameraExposureTime", mxCreateDoubleScalar(pEMVector->GetCameraExposureTime()));
    
    mxSetField(struct_array_ptr, 0, "CameraGain", mxCreateDoubleScalar(pEMVector->GetCameraGainFactor()));
    
    mxSetField(struct_array_ptr, 0, "CameraReadoutRate", mxCreateDoubleScalar(pEMVector->GetCameraReadoutRate()));
    
    mxSetField(struct_array_ptr, 0, "CameraSensitivity", mxCreateDoubleScalar(pEMVector->GetCameraSense()));
    
    mxSetField(struct_array_ptr, 0, "CameraDose", mxCreateDoubleScalar(pEMVector->GetCameraElectronDose()));
    
    mxSetField(struct_array_ptr, 0, "CamGroup", mxCreateDoubleScalar(pEMVector->GetCameraGroup()));
    
    mxSetField(struct_array_ptr, 0, "CamGainIndex", mxCreateDoubleScalar(pEMVector->GetGainIndex()));
    
    mxSetField(struct_array_ptr, 0, "CamSpeedIndex", mxCreateDoubleScalar(pEMVector->GetSpeedIndex()));
    
    mxSetField(struct_array_ptr, 0, "CamCorrectionNrImages", mxCreateDoubleScalar(pEMVector->GetCorrectionNrImages()));
    
    mxSetField(struct_array_ptr, 0, "CamDynamic", mxCreateDoubleScalar(pEMVector->GetCameraDynamic()));
    
    mxSetField(struct_array_ptr, 0, "CamBorderSize", mxCreateDoubleScalar(pEMVector->GetBorder()));
    
    field_value = mxCreateString(pEMVector->GetValue(CAM_SCX_AMPLIFIER).c_str());
    mxSetField(struct_array_ptr, 0, "CamSCXAmplifier", field_value); 
    
    mxSetField(struct_array_ptr, 0, "CamNumPorts", mxCreateDoubleScalar(pEMVector->GetNumberOfPorts()));
    
    field_value = mxCreateString(pEMVector->GetValue(CAM_READOUT_MODE).c_str());
    mxSetField(struct_array_ptr, 0, "CamReadoutMode", field_value); 
    
    field_value = mxCreateString(pEMVector->GetValue(CAM_GEOMETRY).c_str());
    mxSetField(struct_array_ptr, 0, "CamReadoutGeometry", field_value); 
    
    mxSetField(struct_array_ptr, 0, "CamFlatfieldNumber", mxCreateDoubleScalar(pEMVector->GetFlatfieldNumber()));
  
    mxSetField(struct_array_ptr, 0, "CamLSDiffractionCenterX", mxCreateDoubleScalar(pEMVector->GetLSDiffractionCenterX()));
    
    mxSetField(struct_array_ptr, 0, "CamLSDiffractionCenterY", mxCreateDoubleScalar(pEMVector->GetLSDiffractionCenterY()));
    
    field_value = mxCreateString(pEMVector->GetTecnaiMicInfo().c_str());
    mxSetField(struct_array_ptr, 0, "MicroscopeInfo", field_value); 
    
    field_value = mxCreateString(pEMVector->GetTecnaiSpecInfo().c_str());
    mxSetField(struct_array_ptr, 0, "SpecimenInfo", field_value); 

	
	plhs[1] = struct_array_ptr;
    
    //******************************** End ********************************************
    
    //************************* Interpret Image Data **********************************
    switch (lpDst->lType) {
        case DT_UCHAR:          // 1 -> 8Bit (just unsigned)
            BYTE *imgData;
            imgData = (BYTE*)lpDst->vpData;
            
            int i, j;
            for (i = 0; i < lpDst->lDim[0]; i++)
                for (j = 0; j < lpDst->lDim[1]; j++) 
                    outImage[i + lpDst->lDim[0] * j] = (double)imgData[j + lpDst->lDim[1] * i];
 
            break;
        case DT_SHORT:          // 3 -> 16Bit (signed and unsigned)
            if ( pEMVector->GetValue(IMAGE_DATA_TYPE) == "unsigned short ( 16 bit )" ) {
                unsigned short *imgData;
                imgData = (unsigned short*)lpDst->vpData;
                
                int i, j;
                for (i = 0; i < lpDst->lDim[0]; i++)
                    for (j = 0; j < lpDst->lDim[1]; j++) 
                        outImage[i + lpDst->lDim[0] * j] = (double)imgData[j + lpDst->lDim[1] * i];
                        
                break;
            } else if ( pEMVector->GetValue(IMAGE_DATA_TYPE) == "signed short ( 16 bit )" ) {
                signed short *imgData;
                imgData = (signed short*)lpDst->vpData;
                
                int i, j;
                for (i = 0; i < lpDst->lDim[0]; i++)
                    for (j = 0; j < lpDst->lDim[1]; j++) 
                        outImage[i + lpDst->lDim[0] * j] = (double)imgData[j + lpDst->lDim[1] * i];
                        
                break;  
            }
    }
    
/*
	cout << "lDim0:" << lpDst->lDim[0] << "\n";
	cout << "Image Name: " << pEMVector->GetValue(IMAGE_NAME) << "\n";
	cout << "Image Folder: " << pEMVector->GetValue(IMAGE_FOLDER) << "\n";
	cout << "Image Size X: " << pEMVector->GetValue(IMAGE_SIZE_X) << "\n";
	cout << "Image Size Y: " << pEMVector->GetValue(IMAGE_SIZE_Y) << "\n";
	cout << "Image Size Z: " << pEMVector->GetValue(IMAGE_SIZE_Z) << "\n";
	cout << "Image Size E: " << pEMVector->GetValue(IMAGE_SIZE_E) << "\n";
	cout << "Image Data Type: " << pEMVector->GetValue(IMAGE_DATA_TYPE) << "\n";
	cout << "Creation Date: " << pEMVector->GetValue(IMAGE_CREATION_DATE) << "\n";
	cout << "Creation Time: " << pEMVector->GetValue(IMAGE_CREATION_TIME) << "\n";
	cout << "Image Comment: " << pEMVector->GetValue(IMAGE_COMMENT) << "\n";
	cout << "Image History: " << pEMVector->GetImageHistory() << "\n";
	cout << "Image Scaling: " << pEMVector->GetValue(IMAGE_SCALING_LOWER) 
		<< " - " << pEMVector->GetValue(IMAGE_SCALING_UPPER) << "\n"; // Array!!!!!!!!!!!!!!!!!
	cout << "Image Type: " << pEMVector->GetValue(IMAGE_TYPE_PARAM) << "\n";
	cout << "Display Type: " << pEMVector->GetDisplayType() << "\n";
	cout << "Image Distance X: " << pEMVector->GetImageDistX() << "\n";
	cout << "Image Distance Y: " << pEMVector->GetImageDistY() << "\n";
	cout << "Image Distance Z: " << pEMVector->GetImageDistZ() << "\n";
	cout << "Image Distance E: " << pEMVector->GetImageDistE() << "\n";
	cout << "TEM Type: " << pEMVector->GetValue(TEM_TYPE) << "\n";
	cout << "High Tension: " << pEMVector->GetValue(TEM_HT) << "\n";
	cout << "Abberations Cs: " << pEMVector->GetValue(TEM_ABBERATIONS_CS) << "\n";
	cout << "Abberations Cc: " << pEMVector->GetValue(TEM_ABBERATIONS_CC) << "\n";
	cout << "Energy Parameters: " << pEMVector->GetValue(TEM_ENERGY) << "\n";
	cout << "Tem Mode: " << pEMVector->GetValue(TEM_MODE) << "\n";
	cout << "current magnification relative to plate: " << pEMVector->GetValue(TEM_MAGNIFICATION) << "\n";
	cout << "Corrected Magnification: " << pEMVector->GetValue(TEM_MAG_COR) << "\n";
	cout << "Post Magnification: " << pEMVector->GetValue(TEM_POSTMAG) << "\n";
	cout << "Stage Type: " << pEMVector->GetStageType() << "\n";
  */  
    return;
}

