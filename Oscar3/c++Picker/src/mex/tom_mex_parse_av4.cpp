/***********************************************************************//**
 * \file tom_mex_parse_av4.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    04.03.2007
 **************************************************************************/


#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>

#include <mex.h>


#include <tom/corr/config_files.hpp>
#include "tom_mex_helpfcn.h"

/***********************************************************************//**
 *
 **************************************************************************/
typedef void (*t_mexFunction)(const std::string command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


/***********************************************************************//**
 *
 **************************************************************************/
void mexFunction_parse_average_log(const std::string command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    std::string filename;
    std::vector<tom::avg::st_average_result> av;
    std::stringstream ss;

    if (nrhs != 2 || !prhs[1] || !mxIsChar(prhs[1])) {
        ss << "Option '" << command << "' needs one filename as parameter.";
        mexErrMsgTxt(ss.str().c_str());
    }
    if (nlhs > 1) {
        ss << "Option '" << command << "' has too many output arguments.";
        mexErrMsgTxt(ss.str().c_str());
    }

    {
        std::vector<char> c_buff(mxGetNumberOfElements(prhs[1])+2);
        mxGetString(prhs[1], &c_buff[0], c_buff.size()-1);
        filename = &c_buff[0];
    }

    try {
        tom::avg::parse_average(filename, av);
    } catch (std::exception &e) {
        mexErrMsgTxt(e.what());
    }

    #define NFIELDNAMES 6
    const char *fieldnames[NFIELDNAMES] = { "ccval_s", "ccval_1", "ccval_sqrtn", "ccval_n", "use_idx", "fsc" };
    plhs[0] = mxCreateStructMatrix(1, av.size(), NFIELDNAMES, fieldnames);
    #undef NFIELDNAMES
    if (!plhs[0]) {
        ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
    }

    std::size_t i, j;
    mxArray *v;
    for (i=0; i<av.size(); i++) {
        if (av[i].ccval_s.get()) {
            if (!(v = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL))) {
                ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
            }
            *mxGetPr(v) = *av[i].ccval_s;
            mxSetField(plhs[0], i, "ccval_s", v);
        }
        if (av[i].ccval_1.get()) {
            if (!(v = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL))) {
                ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
            }
            *mxGetPr(v) = *av[i].ccval_1;
            mxSetField(plhs[0], i, "ccval_1", v);
        }
        if (av[i].ccval_sqrtn.get()) {
            if (!(v = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL))) {
                ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
            }
            *mxGetPr(v) = *av[i].ccval_sqrtn;
            mxSetField(plhs[0], i, "ccval_sqrtn", v);
        }
        if (av[i].ccval_n.get()) {
            if (!(v = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL))) {
                ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
            }
            *mxGetPr(v) = *av[i].ccval_n;
            mxSetField(plhs[0], i, "ccval_n", v);
        }
        {
            const std::size_t n = av[i].use_idx.size();
            if (!(v = mxCreateNumericMatrix(1, n, mxUINT32_CLASS, mxREAL))) {
                ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
            }
            if (n) {
                uint32_t *puint32 = reinterpret_cast<uint32_t *>(mxGetData(v));
                const std::size_t *psrc = &av[i].use_idx[0];
                for (j=0; j<n; j++) {
                    puint32[j] = psrc[j];
                }
            }
            mxSetField(plhs[0], i, "use_idx", v);
        }
        if (!av[i].fsc.empty()) {
            const std::size_t n = av[i].fsc.size();
            if (!(v = mxCreateNumericMatrix(n, 2, mxDOUBLE_CLASS, mxREAL))) {
                ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
            }
            const std::pair<double, double> *psrc = &av[i].fsc[0];
            double *pdst = mxGetPr(v);
            for (j=0; j<n; j++) {
                pdst[0*n + j] = 1./psrc[j].first;
                pdst[1*n + j] = psrc[j].second;
            }
            mxSetField(plhs[0], i, "fsc", v);
        }
    }
}





/***********************************************************************//**
 *
 **************************************************************************/
void mexFunction_parse_peakfile(const std::string command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    std::string filename;
    std::stringstream ss;

    std::vector<tom::cc_peak<double> > peak_list;
    std::vector<helper::triple<double, double, double> > anglesv;

    std::size_t ntemplates, nparticles;


    if (nrhs != 4) {
        ss << "Option '" << command << "' needs 3 parameters: filename, number_of_templates, number_of_particles.";
        mexErrMsgTxt(ss.str().c_str());
    }
    if (nlhs > 1) {
        ss << "Option '" << command << "' has too many output arguments.";
        mexErrMsgTxt(ss.str().c_str());
    }

    {
        std::vector<char> c_buff(mxGetNumberOfElements(prhs[1])+2);
        mxGetString(prhs[1], &c_buff[0], c_buff.size()-1);
        filename = &c_buff[0];
    }
    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1 ||
        !mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1) {
        ss << "Needs the number of particles and templates as parameter.";
        mexErrMsgTxt(ss.str().c_str());
    }

    try {
        ntemplates = getScalar<std::size_t>(prhs[2]);
        nparticles = getScalar<std::size_t>(prhs[3]);
    } catch (std::bad_cast) {
        mexErrMsgTxt("Error converting the number of tempaltes or particles to an integer.");
    }


    try {
        tom::corr::parse_peaklist(filename, peak_list, anglesv, ntemplates, nparticles);
    } catch (std::exception &e) {
        mexErrMsgTxt(e.what());
    }

    if (peak_list.size() != ntemplates*nparticles || anglesv.size()!=ntemplates*nparticles) {
        mexErrMsgTxt("Unexpected error: peaklist size does not match.");
    }

    #define NFIELDNAMES 4
    const char *fieldnames[NFIELDNAMES] = { "shift", "angle_index", "angles", "ccval" };
    plhs[0] = mxCreateStructMatrix(ntemplates, nparticles, NFIELDNAMES, fieldnames);
    #undef NFIELDNAMES
    if (!plhs[0]) {
        ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
    }


    std::size_t i, j, k;
    mxArray *v;
    std::vector<tom::cc_peak<double> >::const_iterator peak_list_it = peak_list.begin();
    std::vector<helper::triple<double, double, double> >::const_iterator anglesv_it = anglesv.begin();

    k = 0;
    double *pd;
    for (i=0; i<ntemplates; i++) {
        for (j=0; j<nparticles; j++) {
            const tom::cc_peak<double> &p = *peak_list_it;
            const helper::triple<double, double, double> &av = *anglesv_it;
            k = j*ntemplates + i;
            if (p.angle_idx >= 0) {
                if (!(v = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL))) {
                    ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
                }
                *mxGetPr(v) = p.angle_idx+1;
                mxSetField(plhs[0], k, "angle_index", v);

                if (!(v = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL))) {
                    ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
                }
                *mxGetPr(v) = p.val;
                mxSetField(plhs[0], k, "ccval", v);

                if (!(v = mxCreateNumericMatrix(1, 3, mxDOUBLE_CLASS, mxREAL))) {
                    ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
                }
                pd = mxGetPr(v);
                pd[0] = p.x;
                pd[1] = p.y;
                pd[2] = p.z;
                mxSetField(plhs[0], k, "shift", v);

                if (!(v = mxCreateNumericMatrix(1, 3, mxDOUBLE_CLASS, mxREAL))) {
                    ss << "Could not allocate memory for result (" << __LINE__ << ")."; mexErrMsgTxt(ss.str().c_str());
                }
                pd = mxGetPr(v);
                pd[0] = av.x * 57.29577951308232087665461840231273527024;
                pd[1] = av.y * 57.29577951308232087665461840231273527024;
                pd[2] = av.z * 57.29577951308232087665461840231273527024;
                mxSetField(plhs[0], k, "angles", v);
            }
            peak_list_it++;
            anglesv_it++;
        }
    }
}







/***********************************************************************//**
 *
 **************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    std::stringstream ss;
    std::string synopsis;
    std::string command;

    ss << "[OUT =] " << mexFunctionName() << "(command, [IN])";
    synopsis = ss.str();
    ss.clear();

    std::map<std::string, t_mexFunction> m;
    std::map<std::string, t_mexFunction>::const_iterator mit;
    m["parse_average_log"] = &mexFunction_parse_average_log;
    m["parse_peakfile"] = &mexFunction_parse_peakfile;
    int i;

    if (nrhs==0 && nlhs==0) {
        /* Print help */
        mexPrintf("SYNOPSIS: %s\n", synopsis.c_str());
        mexPrintf("mex-file compiled from file \"" __FILE__ "\" at " __DATE__ ", " __TIME__ "\n");
        mexPrintf("valid values for \"command\" are:\n");

        for (i=1, mit=m.begin(); mit!=m.end(); mit++, i++) {
            mexPrintf("   %d: '%s'\n", i, mit->first.c_str());
        }
        return;
    }


    t_mexFunction f_mexFunction = NULL;
    bool invalid_command_text = false;

    if (nrhs >= 1 && mxIsChar(prhs[0])) {
        std::vector<char> c_buff(mxGetNumberOfElements(prhs[0]) + 2);
        if (mxGetString(prhs[0], &c_buff[0], c_buff.size()-1)) {
            mexErrMsgTxt("Unexpected error in mxGetString.\n");
        }
        command = &c_buff[0];

        mit = m.find(command);
        if (mit == m.end()) {
            invalid_command_text = true;
        } else {
            f_mexFunction = mit->second;
        }
    }

    if (!f_mexFunction) {
        ss <<   "SYNOPSIS: " << synopsis.c_str() << "\n"
                "valid values for \"command\" are:\n";
        for (i=1, mit=m.begin(); mit!=m.end(); mit++, i++) {
            ss << "   " << i << ": '" << mit->first << "'.\n";
        }
        mexErrMsgTxt(ss.str().c_str());
    }


    // call the function for the parameter...
    f_mexFunction(command, nlhs, plhs, nrhs, prhs);

}




