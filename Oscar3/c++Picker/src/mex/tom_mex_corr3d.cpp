/***********************************************************************//**
 * \file tom_mex_corr3d.cpp
 * \brief MEX-Function for align2 in tom_mex_corr3d.hpp
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.12.2007
 **************************************************************************/

#include <stdio.h>
#ifdef _WIN32
    #define snprintf _snprintf
#endif


#include <sstream>
#include <typeinfo>
#include <fftw3.h>

#include "mex.h"


#include "tom_mex_helpfcn.h"

#include "helper/auto_vector.hpp"

#include "tom/corr/corr.hpp"
#include "tom/core/volume_fcn.hpp"



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
template<typename T>
void typed_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    std::stringstream synopsis;
    std::stringstream ss;
    synopsis << "[out] = " << mexFunctionName() << "([datatype], size, templates, particles, angles, mask, wedge_ref, wedge_part, type, [...])";


    if (tom::is_double<T>()) {
        tom::setFcnMem(fftw_malloc, fftw_free);
    } else if (tom::is_float<T>()) {
        tom::setFcnMem(fftwf_malloc, fftwf_free);
    } else {
        throw std::runtime_error("Template function is only for single or double floating point");
    }

    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis.str().c_str());
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ "\n");
        return;
    }


    if (nrhs < 7) {
        ss.str(""); ss << mexFunctionName() << " requires arguments: " << synopsis.str();
        mexErrMsgTxt(ss.str().c_str());
    }


    const mxArray *arg;


    std::vector<tom::VolumeContainer<T> *> vintemplates;
    std::vector<tom::VolumeContainer<T> *> vinparticles;
    auto_vector<tom::VolumeContainer<T> > intemplates;
    auto_vector<tom::VolumeContainer<T> > inparticles;
    std::vector<tom::Wedge<T> *> *pwedge_particles;
    std::vector<tom::Wedge<T> *> wedge_particles;
    auto_vector<tom::Wedge<T> > wedge_particles_release;
    std::vector<tom::Wedge<T> *> *pwedge_templates;
    std::vector<tom::Wedge<T> *> wedge_templates;
    auto_vector<tom::Wedge<T> > wedge_templates_release;
    std::auto_ptr<tom::CorrelationHandler<T> > corr_hdl;
    std::size_t dims[3];
    std::size_t ntemplates, nparticles, nangles;

    mxArray *retval_cc = NULL;
    auto_vector<tom::Volume<T> > release_vols_cc;


    std::auto_ptr<tom::Volume<T> > mask;
    mxArray *pangles;
    std::auto_ptr<tom::Volume<double> > angles;

    #define min(a,b) ((a) < (b) ? (a) : (b))
    #define PI 3.141592653589793238512808959406186204433


    ///////////////////////////////////////////////////////////////////
    // Process input parameters.
    ///////////////////////////////////////////////////////////////////
    {
        std::vector<const mxArray *> elements;
        int cnt_ok;
        std::size_t size;
        const mxArray *el;
        const mwSize *dims_;

        ///////////////////////////////////////////////////////////////////
        // The size of the volumes.
        ///////////////////////////////////////////////////////////////////
        arg = prhs[0];
        if (!mxIsNumeric(arg) || mxGetNumberOfElements(arg)!=3 || mxIsComplex(arg)) {
            ss.str(""); ss << mexFunctionName() << ": The first parameter must be the size of the volumes (i.e. a 3-vector).";
            mexErrMsgTxt(ss.str().c_str());
        } else {
            mxArray *el = getDoubleArray(arg);
            const double *d = mxGetPr(el);
            for (int i=0; i<3; i++) {
                if (d[i] < 1 || d[i]!=round(d[i]) || d[i]>0xFFFFFFFF) {
                    ss.str(""); ss << mexFunctionName() << ": The parameter SIZE must be the dimension of the volumes (i.e. a 3-vector of positive integers).";
                    mexErrMsgTxt(ss.str().c_str());
                }
                dims[i] = static_cast<std::size_t>(d[i]);
            }
            mxDestroyArray(el);
        }

        ///////////////////////////////////////////////////////////////////
        // The template(s)
        ///////////////////////////////////////////////////////////////////
        arg = prhs[1];
        if (mxIsChar(arg) || mxIsNumeric(arg)) {
            elements.push_back(arg);
        } else if (mxIsCell(arg)) {
            mwSize numel = mxGetNumberOfElements(arg);
            for (int i=0; i<numel; i++) {
                elements.push_back(mxGetCell(arg, i));
            }
        }
        cnt_ok = 0;
        size = elements.size();
        intemplates.resize(size);
        for (std::size_t i=0; i<size; i++) {
            el = elements.at(i);
            if (mxIsChar(el)) {
                intemplates.assign_direct(i, new tom::VolumeContainerEM<T>(getStringFromMxArray(el)));
                cnt_ok++;
                continue;
            } else if (mxIsNumeric(el) && !mxIsComplex(el) && mxGetNumberOfDimensions(el)==3) {
                dims_ = mxGetDimensions(el);
                ss.str(""); ss << "template " << i;
                if (dims[0]==dims_[0] && dims[1]==dims_[1] && dims[2]==dims_[2]) {
                    if (mxIsSingle(el)) {
                        tom::Volume<float > v(reinterpret_cast<float *>(mxGetData(el)), dims[0], dims[1], dims[2], 0,0,0, false, NULL);
                        intemplates.assign_direct(i, new tom::VolumeContainerMem<T>(ss.str(), v, false));
                        cnt_ok++;
                        continue;
                    } else if (mxIsDouble(el)) {
                        tom::Volume<double> v(mxGetPr(el), dims[0], dims[1], dims[2], 0,0,0, false, NULL);
                        intemplates.assign_direct(i, new tom::VolumeContainerMem<T>(ss.str(), v, false));
                        cnt_ok++;
                        continue;
                    }
                }
            }
            break;
        }
        if (cnt_ok != size || size<1) {
            ss.str(""); ss << mexFunctionName() << ": The content of the parameter TEMPLATES must be either em-filename(s) or real floating point volume(s) of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << ".";
            mexErrMsgTxt(ss.str().c_str());
        }
        ntemplates = intemplates.size();



        ///////////////////////////////////////////////////////////////////
        // The particle(s)
        ///////////////////////////////////////////////////////////////////
        arg = prhs[2];
        elements.clear();
        if (mxIsChar(arg) || mxIsNumeric(arg)) {
            elements.push_back(arg);
        } else if (mxIsCell(arg)) {
            mwSize numel = mxGetNumberOfElements(arg);
            for (int i=0; i<numel; i++) {
                elements.push_back(mxGetCell(arg, i));
            }
        }
        size = elements.size();
        inparticles.resize(size);
        cnt_ok = 0;
        for (std::size_t i=0; i<size; i++) {
            el = elements.at(i);
            if (mxIsChar(el)) {
                std::string filename = getStringFromMxArray(el);
                inparticles.assign_direct(i, new tom::VolumeContainerEM<T>(filename));
                cnt_ok++;
                continue;
            } else if (mxIsNumeric(el) && !mxIsComplex(el) && mxGetNumberOfDimensions(el)==3) {
                dims_ = mxGetDimensions(el);
                ss.str(""); ss << "particle " << i;
                if (dims[0]==dims_[0] && dims[1]==dims_[1] && dims[2]==dims_[2]) {
                    if (mxIsSingle(el)) {
                        tom::Volume<float > v(reinterpret_cast<float *>(mxGetData(el)), dims[0], dims[1], dims[2], 0,0,0, false, NULL);
                        inparticles.assign_direct(i, new tom::VolumeContainerMem<T>(ss.str(), v, false));
                        cnt_ok++;
                        continue;
                    } else if (mxIsDouble(el)) {
                        tom::Volume<double> v(mxGetPr(el), dims[0], dims[1], dims[2], 0,0,0, false, NULL);
                        inparticles.assign_direct(i, new tom::VolumeContainerMem<T>(ss.str(), v, false));
                        cnt_ok++;
                        continue;
                    }
                }
            }
            break;
        }
        if (cnt_ok != size || size<1) {
            ss.str(""); ss << mexFunctionName() << ": The content of the parameter PARTICLES must be either em-filename(s) or real floating point volume(s) of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << ".";
            mexErrMsgTxt(ss.str().c_str());
        }
        nparticles = inparticles.size();

        ///////////////////////////////////////////////////////////////////
        // The angles
        ///////////////////////////////////////////////////////////////////
        arg = prhs[3];
        dims_ = mxGetDimensions(arg);
        if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg)>2 || mxGetNumberOfElements(arg)<1 || dims_[0]!=3) {
            ss.str(""); ss << mexFunctionName() << ": The parameter ANGLES must be a Nx3 matrix of euler-angles [phi, psi, theta] in degrees.";
            mexErrMsgTxt(ss.str().c_str());
        }
        pangles = getDoubleArray(arg);
        angles.reset(new tom::Volume<double>(mxGetPr(pangles), 3, mxGetN(pangles), 1, 0,0,0, false, NULL));
        #define PI 3.141592653589793238512808959406186204433
        angles->template shift_scale<double>(0, PI/180.);
        nangles = angles->getSizeY();


        ///////////////////////////////////////////////////////////////////
        // The mask
        ///////////////////////////////////////////////////////////////////
        arg = prhs[4];
        dims_ = mxGetDimensions(arg);
        if (!mxIsNumeric(arg) || mxIsComplex(arg) ||
            !(mxGetNumberOfElements(arg)==0 || mxGetNumberOfElements(arg)==1 ||
            (mxGetNumberOfDimensions(arg)==3 && dims_[0]==dims[0] && dims_[1]==dims[1] && dims_[2]==dims[2]))) {
            ss.str(""); ss << mexFunctionName() << ": The parameter MASK must be a numeric mask of size " << dims[0] << "x" << dims[1] << "x" << dims[2];
            mexErrMsgTxt(ss.str().c_str());
        }
        if (mxGetNumberOfElements(arg) == 0) {
        } else if (mxGetNumberOfElements(arg) == 1) {
            mask.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], NULL, NULL));
            tom::init_spheremask(*mask, mxGetScalar(arg), 0, 0);
        } else {
            if (mxIsDouble(arg)) {
                if (tom::is_double<T>()) {
                    mask.reset(new tom::Volume<T>(reinterpret_cast<T *>(mxGetData(arg)), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
                } else if (tom::is_float<T>()) {
                    mask.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], NULL, NULL));
                    mask->setValues(tom::Volume<double>(mxGetPr(arg), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
                } else {
                    assert(0);
                }
            } else if (mxIsSingle(arg)) {
                if (tom::is_float<T>()) {
                    mask.reset(new tom::Volume<T>(reinterpret_cast<T *>(mxGetData(arg)), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
                } else if (tom::is_double<T>()) {
                    mask.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], NULL, NULL));
                    mask->setValues(tom::Volume<float>(reinterpret_cast<float *>(mxGetData(arg)), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
                } else {
                    assert(0);
                }
            } else {
                ss.str(""); ss << mexFunctionName() << ": The parameter MASK must be a numeric mask of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << " (single or double).";
                mexErrMsgTxt(ss.str().c_str());
            }
        }

        ///////////////////////////////////////////////////////////////////
        // Wedge of the reference.
        ///////////////////////////////////////////////////////////////////
        arg = prhs[5];
        size = intemplates.size();
        if (mxIsNumeric(arg) && mxGetNumberOfElements(arg)==0) {
            pwedge_templates = NULL;
        } else {
            wedge_templates.resize(size);
            wedge_templates_release.resize(size);
            elements.clear();
            if (mxIsChar(arg) || mxIsNumeric(arg)) {
                elements.push_back(arg);
            } else if (mxIsCell(arg)) {
                mwSize numel = mxGetNumberOfElements(arg);
                for (int i=0; i<numel; i++) {
                    elements.push_back(mxGetCell(arg, i));
                }
            } else {
                ss.str(""); ss << mexFunctionName() << ": WEDGE_REF must be either [] (no wedge), a scalar (angle of the wedge), a volume (the hole weighting in fourier space, fftshifted) or a cell array of these.";
                mexErrMsgTxt(ss.str().c_str());
            }
            if (elements.size()!=1 && elements.size()!=size) {
                ss.str(""); ss << mexFunctionName() << ": The number of WEDGE_REF must be either 0 (none), 1 (the same wedge for all templates) or " << size << " (individual wedge)";
                mexErrMsgTxt(ss.str().c_str());
            }
            std::size_t nwedges = elements.size();
            cnt_ok = 0;
            for (std::size_t i=0; i<nwedges; i++) {
                el = elements.at(i);
                if (mxIsChar(el)) {
                    wedge_templates_release.assign_direct(i, new tom::VolumeEMWedge<T>(getStringFromMxArray(el), 1));
                    cnt_ok++;
                    continue;
                } else if (mxIsNumeric(el) && !mxIsComplex(el)) {
                    if (mxGetNumberOfElements(el) == 0) {
                        wedge_templates_release.assign_direct(i, new tom::NoWedge<T>());
                        cnt_ok++;
                        continue;
                    } else if (mxGetNumberOfElements(el) == 1) {
                        double angle = mxGetScalar(el);
                        wedge_templates_release.assign_direct(i, new tom::SimpleWedge<T>(angle*PI/180., (min(dims[0], min(dims[1],dims[2]))-1) / 2 ));
                        cnt_ok++;
                        continue;
                    } else {
                        dims_ = mxGetDimensions(el);
                        if ((!mxIsSingle(el) && !mxIsDouble(el)) || mxGetNumberOfDimensions(el)!=3 || dims_[0]!=dims[0] || dims_[1]!=dims[1] || dims_[2]!=dims[2]) {
                            ss.str(""); ss << mexFunctionName() << ": If WEDGE_REF is specified as numerical volume, it must be floating point of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << " (single or double).";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        if (mxIsSingle(el)) {
                            tom::Volume<float > v((float *)mxGetData(el), dims[0], dims[1], dims[2], 0,0,0, false, NULL);
                            wedge_templates_release.assign_direct(i, new tom::VolumeWedge<T>(v ,false));
                        } else if (mxIsDouble(el)) {
                            tom::Volume<double> v(mxGetPr(el), dims[0], dims[1], dims[2], 0,0,0, false, NULL);
                            wedge_templates_release.assign_direct(i, new tom::VolumeWedge<T>(v, false));
                        }
                        cnt_ok++;
                        continue;
                    }
                }
                break;
            }
            if (cnt_ok!=nwedges) {
                ss.str(""); ss << mexFunctionName() << ": WEDGE_REF must be either [], angle(s) in degree, em-filename(s) or real floating point volume(s) of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << ".";
                mexErrMsgTxt(ss.str().c_str());
            }
            pwedge_templates = &wedge_templates;
            if (cnt_ok==1 && size>1) {
                for (std::size_t i=0; i<size; i++) { wedge_templates.at(i) = wedge_templates_release.at(0); }
            } else {
                for (std::size_t i=0; i<size; i++) { wedge_templates.at(i) = wedge_templates_release.at(i); }
            }
        }


        ///////////////////////////////////////////////////////////////////
        // Wedge of the particle.
        ///////////////////////////////////////////////////////////////////
        arg = prhs[6];
        size = inparticles.size();
        if (mxIsNumeric(arg) && mxGetNumberOfElements(arg)==0) {
            pwedge_particles = NULL;
        } else {
            wedge_particles.resize(size);
            wedge_particles_release.resize(size);
            elements.clear();
            if (mxIsChar(arg) || mxIsNumeric(arg)) {
                elements.push_back(arg);
            } else if (mxIsCell(arg)) {
                mwSize numel = mxGetNumberOfElements(arg);
                for (int i=0; i<numel; i++) {
                    elements.push_back(mxGetCell(arg, i));
                }
            } else {
                ss.str(""); ss << mexFunctionName() << ": WEDGE_PART must be either [] (no wedge), a scalar (angle of the wedge), a volume (the hole weighting in fourier space, fftshifted) or a cell array of these.";
                mexErrMsgTxt(ss.str().c_str());
            }
            if (elements.size()!=1 && elements.size()!=size) {
                ss.str(""); ss << mexFunctionName() << ": The number of WEDGE_PART must be either 0 (none), 1 (the same wedge for all particles) or " << size << " (individual wedge)";
                mexErrMsgTxt(ss.str().c_str());
            }
            cnt_ok = 0;
            std::size_t nwedges = elements.size();
            for (std::size_t i=0; i<nwedges; i++) {
                el = elements.at(i);
                if (mxIsChar(el)) {
                    wedge_particles_release.assign_direct(i, new tom::VolumeEMWedge<T>(getStringFromMxArray(el), 1));
                    cnt_ok++;
                    continue;
                } else if (mxIsNumeric(el) && !mxIsComplex(el)) {
                    if (mxGetNumberOfElements(el) == 0) {
                        wedge_particles_release.assign_direct(i, new tom::NoWedge<T>());
                        cnt_ok++;
                        continue;
                    } else if (mxGetNumberOfElements(el) == 1) {
                        double angle = mxGetScalar(el);
                        wedge_particles_release.assign_direct(i, new tom::SimpleWedge<T>(angle*PI/180., (min(dims[0], min(dims[1],dims[2]))-1) / 2 ));
                        cnt_ok++;
                        continue;
                    } else {
                        dims_ = mxGetDimensions(el);
                        if ((!mxIsSingle(el) && !mxIsDouble(el)) || mxGetNumberOfDimensions(el)!=3 || dims_[0]!=dims[0] || dims_[1]!=dims[1] || dims_[2]!=dims[2]) {
                            ss.str(""); ss << mexFunctionName() << ": If WEDGE_PART is specified as numerical volume, it must be floating point of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << " (single or double).";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        if (mxIsSingle(el)) {
                            wedge_particles_release.assign_direct(i, new tom::VolumeWedge<T>(tom::Volume<float >((float *)mxGetData(el), dims[0], dims[1], dims[2], 0,0,0, false, NULL), false ));
                        } else if (mxIsDouble(el)) {
                            wedge_particles_release.assign_direct(i, new tom::VolumeWedge<T>(tom::Volume<double>(         mxGetPr  (el), dims[0], dims[1], dims[2], 0,0,0, false, NULL), false ));
                        }
                        cnt_ok++;
                        continue;
                    }
                }
                break;
            }
            if (cnt_ok!=nwedges) {
                ss.str(""); ss << mexFunctionName() << ": WEDGE_PART must be either [], angle(s) in degree, em-filename(s) or real floating point volume(s) of size " << dims[0] << "x" << dims[1] << "x" << dims[2] << ".";
                mexErrMsgTxt(ss.str().c_str());
            }
            pwedge_particles = &wedge_particles;
            if (cnt_ok==1 && size>1) {
                for (std::size_t i=0; i<size; i++) { wedge_particles.at(i) = wedge_particles_release.at(0); }
            } else {
                for (std::size_t i=0; i<size; i++) { wedge_particles.at(i) = wedge_particles_release.at(i); }
            }
        }


        ///////////////////////////////////////////////////////////////////
        // The type...
        ///////////////////////////////////////////////////////////////////
        #define _____TYPE_NONE  ((int)0)
        #define _____TYPE_PEAK  ((int)1)
        #define _____TYPE_CC    ((int)2)
        #define _____TYPE_CC2EM ((int)4)
        {
            int type = _____TYPE_NONE;
            if (nrhs < 8) {
                type = _____TYPE_PEAK;
            } else {
                arg = prhs[7];
                if (mxIsChar(arg)) {
                    std::string t = getStringFromMxArray(arg);
                    if (t == "peak") {
                        type = _____TYPE_PEAK;
                    } else if (t == "cc") {
                        type = _____TYPE_CC;
                    } else if (t == "cc2em") {
                        type = _____TYPE_CC2EM;
                    }
                }
                if (type == _____TYPE_NONE) {
                    ss.str(""); ss << mexFunctionName() << ": TYPE must be a string specifying what to do (currently implemented 'peak' (default), 'cc', 'cc2em').";
                    mexErrMsgTxt(ss.str().c_str());
                }
            }
            switch (type) {
                case _____TYPE_CC:
                    {
                        bool do_fftshift = true;
                        if (nlhs > 1) {
                            ss.str(""); ss << mexFunctionName() << ": Too many output arguments.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        if (nrhs == 9) {
                            arg = prhs[8];
                            if (mxIsNumeric(arg) && mxGetNumberOfElements(arg)==0) {
                            } else if (mxIsLogicalScalar(arg)) {
                                do_fftshift = mxIsLogicalScalarTrue(arg);
                            } else {
                                ss.str(""); ss << mexFunctionName() << ": Operation 'cc' needs as 9th parameter a logical flag whether to get the fft-shifted correlation volume. Ommitting it defaults to true (i.e. do fft-shift).";
                                mexErrMsgTxt(ss.str().c_str());
                            }
                        } else if (nrhs > 9) {
                            ss.str(""); ss << mexFunctionName() << ": Operation 'cc' can only have one additional parameter: a boolean flag whether to do an fftshift before returning the volumes.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        mwSize dims_[3] = { ntemplates, nparticles, nangles };
                        if (!(retval_cc = mxCreateCellArray(3, dims_))) {
                            ss.str(""); ss << mexFunctionName() << ": Error allocating memory for the result.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        std::size_t i,j,k, ii;
                        dims_[0] = dims[0];
                        dims_[1] = dims[1];
                        dims_[2] = dims[2];
                        mxClassID classid = tom::is_double<T>() ? mxDOUBLE_CLASS : mxSINGLE_CLASS;
                        ii = 0;
                        tom::CorrelationHandlerCC<T> *corr_hdl_;
                        corr_hdl.reset(corr_hdl_ = new tom::CorrelationHandlerCC<T>(ntemplates, nparticles, nangles, do_fftshift, dims[0], dims[1], dims[2]));
                        mxArray *el;
                        release_vols_cc.resize(ntemplates*nparticles*nangles);
                        std::auto_ptr<tom::Volume<T> > vv;
                        for (k=0; k<nangles; k++) {
                            for (j=0; j<nparticles; j++) {
                                for (i=0; i<ntemplates; i++) {
                                    if (!(el = mxCreateNumericArray(3, dims_, classid, mxREAL))) {
                                        ss.str(""); ss << mexFunctionName() << ": Error allocating memory for the result.";
                                        mexErrMsgTxt(ss.str().c_str());
                                    }
                                    mxSetCell(retval_cc, ii, el);
                                    vv.reset(new tom::Volume<T>(reinterpret_cast<T *>(mxGetData(el)), dims[0], dims[1], dims[2], 0,0,0, false, NULL));
                                    corr_hdl_->reset(*vv, i, j, k);

                                    release_vols_cc.push_back(vv);
                                    ii++;
                                }
                            }
                        }
                    }
                    break;
                case _____TYPE_CC2EM:
                    {
                        bool do_fftshift = true;
                        if (nlhs > 1) {
                            ss.str(""); ss << mexFunctionName() << ": Too many output arguments.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        if (nrhs < 9) {
                            ss.str(""); ss << mexFunctionName() << ": Operation 'cc2em' needs as parameter a " << ntemplates << "x" << nparticles << " cell array of directory names where to save the cc-volumes.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        arg = prhs[8];
                        if (!mxIsCell(arg) || mxGetNumberOfDimensions(arg)!=2 || mxGetM(arg)!=ntemplates || mxGetN(arg)!=nparticles) {
                            ss.str(""); ss << mexFunctionName() << ": Operation 'cc2em' needs as parameter a " << ntemplates << "x" << nparticles << " cell array of directory names where to save the cc-volumes.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        std::vector<std::string> dirnames(ntemplates*nparticles);
                        {
                            int i,j,k;
                            for (i=0, k=0; i<ntemplates; i++) {
                                for (j=0; j<nparticles; j++) {
                                    el = mxGetCell(arg, k);
                                    if (!mxIsChar(el)) {
                                        ss.str(""); ss << mexFunctionName() << ": The parameter DIRNAME of 'cc2em' must contain all strings (directory names).";
                                        mexErrMsgTxt(ss.str().c_str());
                                    }
                                    std::string n = getStringFromMxArray(el);
                                    dirnames.at(k) = n;
                                    k++;
                                }
                            }
                        }
                        if (nrhs == 10) {
                            arg = prhs[9];
                            if (mxIsNumeric(arg) && mxGetNumberOfElements(arg)==0) {
                            } else if (mxIsLogicalScalar(arg)) {
                                do_fftshift = mxIsLogicalScalarTrue(arg);
                            } else {
                                ss.str(""); ss << mexFunctionName() << ": Operation 'cc2em' needs as 10th parameter a logical flag whether to save the fft-shifted correlation volume. Ommitting it defaults to true (i.e. do fft-shift before saving).";
                                mexErrMsgTxt(ss.str().c_str());
                            }
                        } else if (nrhs > 10) {
                            ss.str(""); ss << mexFunctionName() << ": Operation 'cc2em' can only have two additional parameters: a " << ntemplates << "x" << nparticles << " cell array of directory names and a boolean flag whether to save fft-shifted volumes.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        corr_hdl.reset(new tom::CorrelationHandlerCC2EM<T>(dirnames, ntemplates, nparticles, *angles, do_fftshift));
                    }
                    break;
                default:
                    {
                        if (nlhs > 1) {
                            ss.str(""); ss << mexFunctionName() << ": Too many output arguments.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        if (nrhs < 9) {
                            corr_hdl.reset(new tom::CorrelationHandlerPeaks<T>(intemplates.size(), inparticles.size()));
                        } else {
                            arg = prhs[8];
                            std::auto_ptr<tom::Volume<T> > mask_cc;
                            bool empty = false;
                            if (mxIsNumeric(arg) && !mxIsComplex(arg)) {
                                if (mxGetNumberOfElements(arg) == 0) {
                                    empty = true;
                                } else if (mxGetNumberOfElements(arg) == 1) {
                                    double radi = mxGetScalar(arg);
                                    mask_cc.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], NULL, NULL));
                                    tom::init_spheremask(*mask_cc, radi, 0, 0);
                                } else {
                                    dims_ = mxGetDimensions(el);
                                    if (mxGetNumberOfDimensions(arg)==3 && dims_[0]==dims[0] && dims_[1]==dims[1] && dims_[2]==dims[2]) {
                                        if (mxIsDouble(arg)) {
                                            if (tom::is_double<T>()) {
                                                mask_cc.reset(new tom::Volume<T>(reinterpret_cast<T *>(mxGetPr(arg)), dims[0],dims[1],dims[2], 0,0,0, false, NULL));
                                            } else {
                                                mask_cc.reset(new tom::Volume<T>(tom::Volume<double>(mxGetPr(arg), dims[0],dims[1],dims[2], 0,0,0, false, NULL)));
                                            }
                                        } else if (mxIsSingle(arg)) {
                                            if (tom::is_float<T>()) {
                                                mask_cc.reset(new tom::Volume<T>(reinterpret_cast<T *>(mxGetData(arg)), dims[0],dims[1],dims[2], 0,0,0, false, NULL));
                                            } else {
                                                mask_cc.reset(new tom::Volume<T>(tom::Volume<float>(reinterpret_cast<float *>(mxGetData(arg)), dims[0],dims[1],dims[2], 0,0,0, false, NULL)));
                                            }
                                        }
                                    }
                                }
                            }
                            if (empty) {
                                corr_hdl.reset(new tom::CorrelationHandlerPeaks<T>(intemplates.size(), inparticles.size()));
                            } else {
                                if (!mask_cc.get()) {
                                    ss.str(""); ss << mexFunctionName() << ": TYPE 'peak' needs a correlation mask as parameter: either [], the radius of a sphere in voxels, or the full floating volume.";
                                    mexErrMsgTxt(ss.str().c_str());
                                }
                                corr_hdl.reset(new tom::CorrelationHandlerPeaks<T>(intemplates.size(), inparticles.size(), *mask_cc, true));
                            }
                        }
                    }
                    break;
            }
        }
    }


    for (std::size_t i=0; i<inparticles.size(); i++) {
        vinparticles.push_back(inparticles.at(i));
    }
    for (std::size_t i=0; i<intemplates.size(); i++) {
        vintemplates.push_back(intemplates.at(i));
    }
    unsigned fftw_flags = FFTW_MEASURE;


    //tom::corr3d<T, T>(dims[0], dims[1], dims[2], vinparticles, 0, vintemplates, 0, *angles, mask.get(), pwedge_particles, pwedge_templates, corr_hdl.get(), fftw_flags);
    tom::corr3d_batch<T, T>(dims[0], dims[1], dims[2], vinparticles, 2, vintemplates, 4, *angles, mask.get(), pwedge_particles, pwedge_templates, corr_hdl.get(), fftw_flags);




    ///////////////////////////////////////////////////////////////////
    // Pack the result for matlab...
    ///////////////////////////////////////////////////////////////////
    if (typeid(*corr_hdl.get()) == typeid(tom::CorrelationHandlerPeaks<T>)) {
        // The peak...
        tom::CorrelationHandlerPeaks<T> *corr_hdl_ = dynamic_cast<tom::CorrelationHandlerPeaks<T> *>(corr_hdl.get());

        std::vector<std::vector<tom::cc_peak<T> > > list_peak = corr_hdl_->getPeakList();

        const char *fieldnames[6] = {"peak", "ncc_val", "angle", "angle_idx", "template", "particle"};

        const mwSize dims_res[2] = { ntemplates, nparticles };
        if (!(plhs[0] = mxCreateStructArray(2, dims_res, 6, fieldnames))) {
            ss.str(""); ss << mexFunctionName() << ": Error allocating memory for the result.";
            mexErrMsgTxt(ss.str().c_str());
        }
        std::size_t i, j, k;
        mxArray *mxpeak, *mxncc_val, *mxangle, *mxangle_idx, *mxtemplate, *mxtemplate_idx, *mxparticle, *mxparticle_idx;
        k = 0;
        for (j=0; j<nparticles; j++) {
            for (i=0; i<ntemplates; i++) {
                if (!(mxpeak = mxCreateNumericMatrix(1, 3, mxDOUBLE_CLASS, mxREAL)) ||
                    !(mxncc_val = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL)) ||
                    !(mxangle = mxCreateNumericMatrix(1, 3, mxDOUBLE_CLASS, mxREAL)) ||
                    !(mxangle_idx = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL)) ||
                    !(mxtemplate = mxCreateString(intemplates.at(i)->getName().c_str())) ||
                    !(mxparticle = mxCreateString(inparticles.at(j)->getName().c_str()))) {
                    ss.str(""); ss << mexFunctionName() << ": Error allocating memory for the result.";
                    mexErrMsgTxt(ss.str().c_str());
                }
                if (list_peak.at(i).at(j).angle_idx >= 0) {
                    reinterpret_cast<double *>(mxGetData(mxpeak))[0] = list_peak.at(i).at(j).x;
                    reinterpret_cast<double *>(mxGetData(mxpeak))[1] = list_peak.at(i).at(j).y;
                    reinterpret_cast<double *>(mxGetData(mxpeak))[2] = list_peak.at(i).at(j).z;
                    mxGetPr(mxncc_val)[0] = list_peak.at(i).at(j).val;
                    std::size_t idx = list_peak.at(i).at(j).angle_idx;
                    reinterpret_cast<uint32_t *>(mxGetData(mxangle_idx))[0] = idx + 1;
                    mxGetPr(mxangle)[0] = angles->get(0,idx,0) * 180. / PI;
                    mxGetPr(mxangle)[1] = angles->get(1,idx,0) * 180. / PI;
                    mxGetPr(mxangle)[2] = angles->get(2,idx,0) * 180. / PI;
                } else {
                    reinterpret_cast<double *>(mxGetData(mxpeak))[0] = 0;
                    reinterpret_cast<double *>(mxGetData(mxpeak))[1] = 0;
                    reinterpret_cast<double *>(mxGetData(mxpeak))[2] = 0;
                    mxGetPr(mxncc_val)[0] = mxGetNaN();
                    reinterpret_cast<uint32_t *>(mxGetData(mxangle_idx))[0] = 0;
                    mxGetPr(mxangle)[0] = mxGetNaN();
                    mxGetPr(mxangle)[1] = mxGetNaN();
                    mxGetPr(mxangle)[2] = mxGetNaN();
                }

                mxSetField(plhs[0], k, "peak", mxpeak);
                mxSetField(plhs[0], k, "ncc_val", mxncc_val);
                mxSetField(plhs[0], k, "angle_idx", mxangle_idx);
                mxSetField(plhs[0], k, "angle", mxangle);
                mxSetField(plhs[0], k, "template", mxtemplate);
                mxSetField(plhs[0], k, "particle", mxparticle);
                k++;
            }
        }
    } else if (typeid(*corr_hdl.get()) == typeid(tom::CorrelationHandlerCC2EM<T>)) {
        const char *fieldnames[6] = { "filename", "status", "template", "particle", "angle", "angle_idx" };

        const mwSize dims_res[3] = { ntemplates, nparticles, nangles };
        if (!(plhs[0] = mxCreateStructArray(3, dims_res, 6, fieldnames))) {
            ss.str(""); ss << mexFunctionName() << ": Error allocating memory for the result.";
            mexErrMsgTxt(ss.str().c_str());
        }
        {
            std::size_t i,j,k,ii;
            std::pair<std::string, bool> fs;
            tom::CorrelationHandlerCC2EM<T> *corr_hdl_ = dynamic_cast<tom::CorrelationHandlerCC2EM<T> *>(corr_hdl.get());
            mxArray *mxfilename, *mxstatus, *mxtemplate, *mxparticle, *mxangle, *mxangle_idx;

            ii = 0;
            for (k=0; k<nangles; k++) {
                for (j=0; j<nparticles; j++) {
                    for (i=0; i<ntemplates; i++) {
                        fs = corr_hdl_->getFileStatus(i, j, k);
                        if (!(mxfilename = mxCreateString(fs.first.c_str())) ||
                            !(mxstatus = mxCreateLogicalScalar(fs.second)) ||
                            !(mxtemplate = mxCreateString(intemplates.at(i)->getName().c_str())) ||
                            !(mxparticle = mxCreateString(inparticles.at(j)->getName().c_str())) ||
                            !(mxangle = mxCreateNumericMatrix(1, 3, mxDOUBLE_CLASS, mxREAL)) ||
                            !(mxangle_idx = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL))) {
                            ss.str(""); ss << mexFunctionName() << ": Error allocating memory for the result.";
                            mexErrMsgTxt(ss.str().c_str());
                        }
                        reinterpret_cast<uint32_t *>(mxGetData(mxangle_idx))[0] = k+1;
                        mxGetPr(mxangle)[0] = angles->get(0,k,0) * 180. / PI;
                        mxGetPr(mxangle)[1] = angles->get(1,k,0) * 180. / PI;
                        mxGetPr(mxangle)[2] = angles->get(2,k,0) * 180. / PI;

                        mxSetField(plhs[0], ii, "filename", mxfilename);
                        mxSetField(plhs[0], ii, "status", mxstatus);
                        mxSetField(plhs[0], ii, "angle_idx", mxangle_idx);
                        mxSetField(plhs[0], ii, "angle", mxangle);
                        mxSetField(plhs[0], ii, "template", mxtemplate);
                        mxSetField(plhs[0], ii, "particle", mxparticle);
                        ii++;
                    }
                }
            }
        }
    } else if (typeid(*corr_hdl.get()) == typeid(tom::CorrelationHandlerCC<T>)) {

        tom::CorrelationHandlerCC<T> *corr_hdl_ = static_cast<tom::CorrelationHandlerCC<T> *>(corr_hdl.get());
        mxArray *el;
        release_vols_cc.resize(ntemplates*nparticles*nangles);
        std::auto_ptr<tom::Volume<T> > vv;
        std::size_t i,j,k,ii;
        ii = 0;
        for (k=0; k<nangles; k++) {
            for (j=0; j<nparticles; j++) {
                for (i=0; i<ntemplates; i++) {
                    if (!corr_hdl_->get_success(i,j,k)) {
                        mxSetCell(retval_cc, ii, NULL);
                    }
                    ii++;
                }
            }
        }

        plhs[0] = retval_cc;
    }

}




/***********************************************************************//**
 * \brief
 **************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    typedef float Tdefault;
    if (nrhs > 0 && mxIsChar(prhs[0])) {
        std::string param1 = getStringFromMxArray(prhs[0]);
        if (param1 == "double") {
            typed_mexFunction<double>(nlhs, plhs, nrhs-1, &prhs[1]);
            return;
        } else if (param1 == "single") {
            typed_mexFunction<float >(nlhs, plhs, nrhs-1, &prhs[1]);
            return;
        } else {
            std::stringstream ss;
            ss.str(""); ss << mexFunctionName() << ": The first parameter can be (optionally) a string specifying the precision of the computation (either 'single' or 'double').";
            if (typeid(Tdefault) == typeid(double)) {
                ss << " Ommitting defaults to double.";
            } else if (typeid(Tdefault) == typeid(float)) {
                ss << " Ommitting defaults to single.";
            }
            mexErrMsgTxt(ss.str().c_str());
        }
    }

    typed_mexFunction<Tdefault>(nlhs, plhs, nrhs, prhs);

}



