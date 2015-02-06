/***********************************************************************//**
 * \file tom_mex_rotate.cpp
 * \brief MEX-Function
 * \author  Thomas Haller
 * \version 0.2
 * \date    12.12.2007
 **************************************************************************/

#include <sstream>
#include <iostream>
#include <typeinfo>

#include <mex.h>


#include <tom/tools/mex_helpfcn.h>

#include <tom/volume.hpp>
#include <tom/transf/transform.hpp>



#define PI 3.141592653589793238512808959406186204433

/***********************************************************************//**
 * \brief
 **************************************************************************/
std::string getStringFromMxArray(const mxArray *m) {

    std::string res;

    if (mxIsChar(m)) {
        mwSize numel = mxGetNumberOfElements(m);
        std::vector<char> c(numel+3);
        mxGetString(m, &c.at(0), numel+2);
        res = &c.at(0);
    }
    return res;
}



/***********************************************************************//**
 * \brief
 **************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    std::stringstream synopsis;
    std::stringstream ss;
    synopsis << "[rotvol] = " << mexFunctionName() << "(vol, rotmatrix, [type, center, shift])";

    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis.str().c_str());
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ "\n");
        return;
    }

    if (nlhs >= 2 || nrhs>5 || nrhs<2) {
        ss.str(""); ss << "SYNOPSIS: " << synopsis.str();
        mexErrMsgTxt(ss.str().c_str());
    }

    const mxArray *arg;

    std::auto_ptr<tom::Volume<float > > v_f;
    std::auto_ptr<tom::Volume<double> > v_d;
    std::auto_ptr<tom::Volume<float > > rv_f;
    std::auto_ptr<tom::Volume<double> > rv_d;

    arg = prhs[0];
    mwSize dims[3];
    if (!mxIsNumeric(arg) || mxGetNumberOfDimensions(arg)>3 || mxGetNumberOfElements(arg)<1 || mxIsComplex(arg)) {
        ss.str(""); ss << "The first parameter must be a real volume.";
        mexErrMsgTxt(ss.str().c_str());
    }
    {
        const mwSize *dims2 = mxGetDimensions(arg);
        dims[0] = dims2[0];
        dims[1] = dims2[1];
        if (mxGetNumberOfDimensions(arg) == 3) {
            dims[2] = dims2[2];
        } else {
            dims[2] = 1;
        }
    }

    mxArray *plhs0;

    bool isdouble = false;
    if (mxIsDouble(arg)) {
        isdouble = true;
        plhs0 = mxCreateNumericArray(3, dims,   mxDOUBLE_CLASS, mxREAL);
        v_d.reset(new tom::Volume<double>(mxGetPr(arg), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
        rv_d.reset(new tom::Volume<double>(mxGetPr(plhs0), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
    } else if (mxIsSingle(arg)) {
        plhs0 = mxCreateNumericArray(3, dims,   mxSINGLE_CLASS, mxREAL);
        v_f.reset(new tom::Volume<float>(reinterpret_cast<float *>(mxGetData(arg)), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
        rv_f.reset(new tom::Volume<float>(reinterpret_cast<float *>(mxGetData(plhs0)), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
    } else {
        ss.str(""); ss << "The first parameter must be a floating point, real volume.";
        mexErrMsgTxt(ss.str().c_str());
    }

    double P[16] = { 0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 1. };
    arg = prhs[1];
    if (!mxIsNumeric(arg) || mxGetNumberOfDimensions(arg)!=2 || mxIsComplex(arg)) {
        ss.str(""); ss << mexFunctionName() << ": The second parameter must be the rotation angles (3-vector) or a rotation matrix.";
        mexErrMsgTxt(ss.str().c_str());
    } else {
        const mwSize *dims = mxGetDimensions(arg);
        mxArray *mxp;
        if ((dims[0]==1 && dims[1]==3 ||
             dims[0]==3 && dims[1]==1 ||
             dims[0]==3 && dims[1]==3) &&
             (mxp=getDoubleArray(arg))) {
            const double *PP = mxGetPr(mxp);
            if (dims[0]*dims[1]==3) {
                const int axes[3] = { 2, 0, 2 };
                const double angles[3] = { -PP[0]*(PI/180.), -PP[2]*(PI/180.), -PP[1]*(PI/180.) };
                tom::transf::sum_rotation(P, true, false, 0, 0, 3, angles, axes, 0);
            } else {
                P[0] = PP[0];
                P[1] = PP[1];
                P[2] = PP[2];
                P[4] = PP[3];
                P[5] = PP[4];
                P[6] = PP[5];
                P[8] = PP[6];
                P[9] = PP[7];
                P[10] = PP[8];
                tom::transf::invert_4x4(P, P, 4);
            }
        } else {
            ss.str(""); ss << mexFunctionName() << ": The second parameter must be the rotation angles (3-vector) or a rotation matrix.";
            mexErrMsgTxt(ss.str().c_str());
        }
    }

    #define INTPOL_TRILINEAR                1
    #define INTPOL_NEARESTNEIGHBOUR         2
    #define INTPOL_INT3D_NEAREST            3
    #define INTPOL_INT3D_LINEAR             4
    #define INTPOL_INT3D_CUBIC              5
    #define INTPOL_INT3D_CSPLINE            6

    int intp_type = INTPOL_TRILINEAR;
    short interpolate3d_method;

    if (nrhs>=3) {
        arg = prhs[2];
        std::vector<char> cbuff(30);
        if (mxIsChar(arg) && mxGetNumberOfElements(arg)<cbuff.size()-1) {
            mxGetString(arg, &cbuff[0], cbuff.size()-1);
            std::string type = &cbuff[0];
            if (type == "linear") {
                intp_type = INTPOL_TRILINEAR;
            } else if (type == "nearest") {
                intp_type = INTPOL_NEARESTNEIGHBOUR;
            } else {
                ss.str(""); ss << mexFunctionName() << ": The third parameter is the type of interpolation: currently 'linear', 'nearest' implemented.";
                mexErrMsgTxt(ss.str().c_str());
            }
        } else if ((mxIsChar(arg)||mxIsNumeric(arg)) && mxGetNumberOfElements(arg)==0) {

        } else {
            ss.str(""); ss << mexFunctionName() << ": The third parameter is the type of interpolation: currently 'linear' and 'nearest' implemented.";
            mexErrMsgTxt(ss.str().c_str());
        }
    }


    double center[3] = { dims[0]/2, dims[1]/2, dims[2]/2 };
    if (nrhs>=4) {
        arg = prhs[3];
        if (!mxIsNumeric(arg) || mxGetNumberOfDimensions(arg)!=2 || mxIsComplex(arg)) {
            ss.str(""); ss << mexFunctionName() << ": The fourth parameter must be the center of the rotation (3-vector).";
            mexErrMsgTxt(ss.str().c_str());
        } else {
            const mwSize *dims = mxGetDimensions(arg);
            mxArray *mxp;
            if ((dims[0]==1 && dims[1]==3 ||
                 dims[0]==3 && dims[1]==1) &&
                (mxp=getDoubleArray(arg))) {
                center[0] = mxGetPr(arg)[0];
                center[1] = mxGetPr(arg)[1];
                center[2] = mxGetPr(arg)[2];
            } else {
                ss.str(""); ss << mexFunctionName() << ": The fourth parameter must be the center of the rotation (3-vector).";
                mexErrMsgTxt(ss.str().c_str());
            }
        }
    }


    double shift[3] = { 0, 0, 0 };
    if (nrhs>=5) {
        arg = prhs[4];
        if (!mxIsNumeric(arg) || mxGetNumberOfDimensions(arg)!=2 || mxIsComplex(arg)) {
            ss.str(""); ss << mexFunctionName() << ": The fifth parameter must be the shift after rotating (3-vector).";
            mexErrMsgTxt(ss.str().c_str());
        } else {
            const mwSize *dims = mxGetDimensions(arg);
            mxArray *mxp;
            if ((dims[0]==1 && dims[1]==3 ||
                 dims[0]==3 && dims[1]==1) &&
                (mxp=getDoubleArray(arg))) {
                shift[0] = mxGetPr(arg)[0];
                shift[1] = mxGetPr(arg)[1];
                shift[2] = mxGetPr(arg)[2];
            } else {
                ss.str(""); ss << mexFunctionName() << ": The fifth parameter must be the shift after rotating (3-vector).";
                mexErrMsgTxt(ss.str().c_str());
            }
        }
    }


    P[ 3] = shift[0] + center[0] + P[ 0]*(-center[0]) + P[ 1]*(-center[1]) + P[ 2]*(-center[2]);
    P[ 7] = shift[1] + center[1] + P[ 4]*(-center[0]) + P[ 5]*(-center[1]) + P[ 6]*(-center[2]);
    P[11] = shift[2] + center[2] + P[ 8]*(-center[0]) + P[ 9]*(-center[1]) + P[10]*(-center[2]);

    #if 0
    std::cout << "rot_mat = [" << std::endl;
    std::cout << "  [" << P[0] << ", " << P[1] << ", " << P[ 2] << ", " << P[ 3] << "];" << std::endl;
    std::cout << "  [" << P[4] << ", " << P[5] << ", " << P[ 6] << ", " << P[ 7] << "];" << std::endl;
    std::cout << "  [" << P[8] << ", " << P[9] << ", " << P[10] << ", " << P[11] << "];" << std::endl;
    std::cout << " ];" << std::endl;
    #endif


    if (isdouble) {
        if (intp_type == INTPOL_TRILINEAR) {
            tom::transf::transform<double, tom::transf::InterpolTriLinear<double> >(*v_d, *rv_d, P, true, 0, tom::transf::InterpolTriLinear<double>(0));
        } else if (intp_type == INTPOL_NEARESTNEIGHBOUR) {
            tom::transf::transform<double, tom::transf::InterpolNearestNeighbour<double> >(*v_d, *rv_d, P, true, 0, tom::transf::InterpolNearestNeighbour<double>(0));
        }
    } else {
        if (intp_type == INTPOL_TRILINEAR) {
            tom::transf::transform<float , tom::transf::InterpolTriLinear<float > >(*v_f, *rv_f, P, true, 0, tom::transf::InterpolTriLinear<float >(0));
        } else if (intp_type == INTPOL_NEARESTNEIGHBOUR) {
            tom::transf::transform<float , tom::transf::InterpolNearestNeighbour<float > >(*v_f, *rv_f, P, true, 0, tom::transf::InterpolNearestNeighbour<float >(0));
        }
    }

    plhs[0] = plhs0;

}




