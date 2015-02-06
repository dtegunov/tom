/****************************************************************************//**
 * \file fftw_plan.cpp
 * \brief Contains implementations of the functions related to the fourier-transformation.
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.11.2007
 *
 *******************************************************************************/
#include "tom/core/fftw_plan.hpp"



#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <exception>
#include <iostream>
#include <cctype>
#include <algorithm>

#include <boost/lexical_cast.hpp>



#define PARAM2STR(x) #x
#define DEF2STR(x) PARAM2STR(x)
#define DEFCONCAT(prefix, name) prefix ## name


#include "tom/core/volume_fcn.hpp"
#include <boost/shared_ptr.hpp>


/****************************************************************************//**
 * \brief Variable for saving the wisdom-filename.
 *******************************************************************************/
char *wisdom_filename = NULL;



namespace {
template<typename T>
class plan_destroyer {
public:
	void operator()(void *p);
};
template<>
class plan_destroyer<float > {
public:
	void operator()(void *p) {
		if (p) {
			fftwf_destroy_plan(reinterpret_cast<fftwf_plan>(p));
		}
	}
};
template<>
class plan_destroyer<double> {
public:
	void operator()(void *p) {
		if (p) {
			fftw_destroy_plan (reinterpret_cast<fftw_plan >(p));
		}
	}
};

}






/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::fftw::flag2str(unsigned flag) {
    flag = flag & (FFTW_ESTIMATE | FFTW_MEASURE | FFTW_PATIENT | FFTW_EXHAUSTIVE);
    if (flag == FFTW_ESTIMATE) {
        return "FFTW_ESTIMATE";
    } else if (flag == FFTW_MEASURE) {
        return "FFTW_MEASURE";
    } else if (flag == FFTW_PATIENT) {
        return "FFTW_PATIENT";
    } else if (flag == FFTW_EXHAUSTIVE) {
        return "FFTW_EXHAUSTIVE";
    } else {
        throw std::runtime_error("unrecognised fftw-flag \"" + boost::lexical_cast<std::string>(flag) + "\"");
    }
}


/****************************************************************************//**
 *
 *******************************************************************************/
unsigned tom::fftw::str2flag(const std::string &s) {
    std::string S(s);
    std::transform(S.begin(), S.end(), S.begin(), &toupper);
    if (S == "FFTW_PATIENT") {
        return FFTW_PATIENT;
    } else if (S == "FFTW_MEASURE") {
        return FFTW_MEASURE;
    } else if (S == "FFTW_ESTIMATE") {
        return FFTW_ESTIMATE;
    } else if (S == "FFTW_EXHAUSTIVE") {
        return FFTW_EXHAUSTIVE;
    } else {
        throw std::runtime_error("unrecognised string for fftw-flag \"" + s + "\"");
    }
}


/****************************************************************************//**
 * \brief Returns the currently set wisdom filename.
 *
 * Do not modify the resulting pointer.\n
 * These wisdom-functions are not thread save!
 *******************************************************************************/
const char *tom::fftw::get_wisdom_name() {
    return wisdom_filename;
}


/****************************************************************************//**
 * \brief Load the wisdom from file and memorize its name in
 *
 * \param[in] fname The name of the wisdom file. Although this name is
 *   saved, the handled argument can be savely freed/modified afterwards.
 *   It is save to pass get_fftw_wisdom_name() as input (self-assignment).
 * \param[in] load If true, the wisdom is reloaded from the file. Otherwise
 *   its name is just remembered.
 * \return 1 in case of success, 0 otherwise. In case of failure,
 *   probably load was true, but the file contains invalid data.
 *   Or maybe (less likely), allocating memory for the copy of fname failed.
 *
 * For example you can set the name, and load it later be calling
 * \code
 * set_fftw_wisdom_name(name, 0);
 * ...
 * set_fftw_wisdom_name(get_fftw_wisdom_name(), 1);
 * \endcode
 *
 * These wisdom-functions are not thread save!
 *******************************************************************************/
int tom::fftw::set_wisdom_name(const char *fname, int load) {

    int res = 0;
    char *wisdom_filename_local = NULL;


    if (fname) {
        wisdom_filename_local = (char *)malloc(sizeof(char)*strlen(fname)+1);
        if (!wisdom_filename_local) {
            return 0;
        }
        strcpy(wisdom_filename_local, fname);
    }
    tom::fftw::clear_wisdom_name();
    if ((wisdom_filename=wisdom_filename_local) && load) { /* This assignment is intended :) */
        FILE *f = fopen(fname, "rb");
        if (f) {
            res = fftw_import_wisdom_from_file(f);
            fclose(f);
        }
    }
    return res;
}


/****************************************************************************//**
 * \brief Saves the wisdom to the file with name get_fftw_wisdom_name()
 *
 * \return 1 in case of success, 0 otherwise. Failure can indicate, that
 *   the filename was not set with set_fftw_wisdom_name or the file could not
 *   be written.
 *
 * These wisdom-functions are not thread save!
 *******************************************************************************/
int tom::fftw::save_wisdom() {
    int res = 0;
    if (wisdom_filename) {
        FILE *f = fopen(wisdom_filename, "wb");
        if (f) {
            fftw_export_wisdom_to_file(f);
            res = !ferror(f);
            fclose(f);
        }
    }
    return res;
}



/****************************************************************************//**
 * \brief Saves the wisdom to the file with name get_fftw_wisdom_name()
 *
 * \return 1 in case of success, 0 otherwise. Failure can indicate, that
 *   the filename was not set with set_fftw_wisdom_name or the file could not
 *   be written.
 *
 * These wisdom-functions are not thread save!
 *******************************************************************************/
void tom::fftw::clear_wisdom_name() {
    if (wisdom_filename) {
        free(wisdom_filename);
        wisdom_filename = NULL;
    }
}








/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
int tom::fftw::Plan<T>::is_valid_dft_r2c(const tom::Volume<T> &vsrc, const tom::Volume<std::complex<T> > &vdst) {
    if (vsrc.getStrideX()%sizeof(T) == 0 &&
        vsrc.getStrideY()%sizeof(T) == 0 &&
        vsrc.getStrideZ()%sizeof(T) == 0 &&
        vdst.getStrideX()%sizeof(std::complex<T>) == 0 &&
        vdst.getStrideY()%sizeof(std::complex<T>) == 0 &&
        vdst.getStrideZ()%sizeof(std::complex<T>) == 0) {
        if (vsrc.getSizeX()      == vdst.getSizeX() &&
             vsrc.getSizeY()      == vdst.getSizeY() &&
            (vsrc.getSizeZ()/2+1) == vdst.getSizeZ()) {
            return TOM_FFTW_PLAN_3D;
        } else if (vsrc.getSizeX()      == vdst.getSizeX() &&
                   (vsrc.getSizeY()/2+1) == vdst.getSizeY() &&
                   vsrc.getSizeZ()==1 && vdst.getSizeZ()==1) {
            return TOM_FFTW_PLAN_2D;
        }
    }
    return TOM_FFTW_PLAN_INVALID;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
int tom::fftw::Plan<T>::is_valid_dft_c2r(const tom::Volume<std::complex<T> > &vsrc, const tom::Volume<T> &vdst) {
    if (vsrc.getStrideX()%sizeof(std::complex<T>) == 0 &&
        vsrc.getStrideY()%sizeof(std::complex<T>) == 0 &&
        vsrc.getStrideZ()%sizeof(std::complex<T>) == 0 &&
        vdst.getStrideX()%sizeof(T) == 0 &&
        vdst.getStrideY()%sizeof(T) == 0 &&
        vdst.getStrideZ()%sizeof(T) == 0) {
        if (vsrc.getSizeX() == vdst.getSizeX() &&
            vsrc.getSizeY() == vdst.getSizeY() &&
            vsrc.getSizeZ() == (vdst.getSizeZ()/2+1)) {
            return TOM_FFTW_PLAN_3D;
        } else if (vsrc.getSizeX() == vdst.getSizeX() &&
                   vsrc.getSizeY() == (vdst.getSizeY()/2+1) &&
                   vsrc.getSizeZ()==1 && vdst.getSizeZ()==1) {
            return TOM_FFTW_PLAN_2D;
        }
    }
    return TOM_FFTW_PLAN_INVALID;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::Plan(tom::fftw::Plan<T> &originalPlan){


	if(originalPlan.plan.get() == NULL)
		throw std::invalid_argument("tom::fftw::Plan(tom::fftw::Plan<T> &originalPlan) - The original plan is empty!");

	std::auto_ptr<tom::Volume<T               > > p_time0(new tom::Volume<T               >(*originalPlan.p_time0,NULL,originalPlan.p_time0->getSizeX(),originalPlan.p_time0->getSizeY(),originalPlan.p_time0->getSizeZ(),originalPlan.p_time0->getStrideX(),originalPlan.p_time0->getStrideY(),originalPlan.p_time0->getStrideZ()));

    std::auto_ptr<tom::Volume<std::complex<T> > > p_freq0(new tom::Volume<std::complex<T> >(*originalPlan.p_freq0,NULL,originalPlan.p_freq0->getSizeX(),originalPlan.p_freq0->getSizeY(),originalPlan.p_freq0->getSizeZ(),originalPlan.p_freq0->getStrideX(), originalPlan.p_freq0->getStrideY(), originalPlan.p_freq0->getStrideZ()));

	this->p_time0 = p_time0.release();
	this->p_freq0 = p_freq0.release();
	this->p_time1 = NULL;
	this->p_freq1 = NULL;

	this->type = originalPlan.type;
	this->plan = originalPlan.plan;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::Plan(tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst, unsigned flags){

    const int type = this->is_valid_dft_r2c(vsrc, vdst);
    if (type == TOM_FFTW_PLAN_INVALID) {
        throw std::invalid_argument("tom::fftw::Plan() - The input volumes have not the right size or alignment for dft_r2c");
    }


    int rank = (type==TOM_FFTW_PLAN_3D ? 3 : 2);
    int howmany_rank = 0;
    fftw_iodim iodim[3];
    iodim[0].n  = vsrc.getSizeX();
    iodim[0].is = vsrc.getStrideX()/sizeof(T);
    iodim[0].os = vdst.getStrideX()/sizeof(std::complex<T>);
    iodim[1].n  = vsrc.getSizeY();
    iodim[1].is = vsrc.getStrideY()/sizeof(T);
    iodim[1].os = vdst.getStrideY()/sizeof(std::complex<T>);
    iodim[2].n  = vsrc.getSizeZ();
    iodim[2].is = vsrc.getStrideZ()/sizeof(T);
    iodim[2].os = vdst.getStrideZ()/sizeof(std::complex<T>);
    fftw_iodim *howmany_dims = NULL;


	boost::shared_ptr<void> p(static_cast<void *>(0), ::plan_destroyer<T>());
	if (tom::is_float<T>()) {
		p = boost::shared_ptr<void>(static_cast<void *>(fftwf_plan_guru_dft_r2c(rank, iodim, howmany_rank, howmany_dims, (float  *)&vsrc.get(), (fftwf_complex *)&vdst.get(), flags)), ::plan_destroyer<T>());
	} else {
		p = boost::shared_ptr<void>(static_cast<void *>(fftw_plan_guru_dft_r2c (rank, iodim, howmany_rank, howmany_dims, (double *)&vsrc.get(), (fftw_complex  *)&vdst.get(), flags)), ::plan_destroyer<T>());
	}
	if (!p.get()) {
		throw std::runtime_error("tom::fftw::Plan() - could not create fftw plan dft_r2c");
	}

    std::auto_ptr<tom::Volume<std::complex<T> > > p_freq0(new tom::Volume<std::complex<T> >(vdst, NULL, vdst.getSizeX(), vdst.getSizeY(), vdst.getSizeZ(), vdst.getStrideX(), vdst.getStrideY(), vdst.getStrideZ()));
    std::auto_ptr<tom::Volume<T               > > p_time0(new tom::Volume<T               >(vsrc, NULL, vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ(), vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ()));

	this->plan.swap(p);
    this->type = tom::fftw::Plan<T>::FFTW_DFT_R2C;
    this->p_freq0 = p_freq0.release();
    this->p_time0 = p_time0.release();
    this->p_freq1 = NULL;
    this->p_time1 = NULL;
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::~Plan() {
    delete this->p_freq0;
    delete this->p_time0;
    delete this->p_freq1;
    delete this->p_time1;

}

/****************************************************************************//**
 * \brief Constructor of a C2R plan.
 *
 * \param[in,out] vsrc The source volume.
 * \param[in,out] vdst The destination volume.
 * \param[in] flags Flags passed directly to fftw_plan_guru_dft_c2r.
 *
 * Calles fftw_plan_guru_dft_c2r to create an fftw_plan. Depending on the flags
 * the input volumes are destroyed. See the documentation of fftw_plan_guru_dft_c2r.\n
 * Unfortunately there seems to be a bug in fftw3.1 concerning fftw_plan_guru_dft_c2r with
 * FFTW_PATIENT and FFTW_EXHAUSTIVE. The function writes outside the boundaries of
 * vdst. One sollution could be to initialize vdst with last dimension vsrc.getSizeZ()
 * large, instead of vsrc.getSizeZ()/2+1. If there is already wisdom about this type
 * of transformation, the error also don't happens.\n
 * It is not possible to use fftw_plan_dft_c2r_3d instead of the guru version because
 * fftw_plan_dft_c2r_3d expects the first dimension (vsrc.getSizeX()) to be reduced
 * instead of the last. Also the dimensions x,y,z are swapped. I think thats
 * an error in the FFTW-api, because its not consistent with the guru interface nor
 * the documentation.
 * With fftw3.2 this bug seams to be fixed.
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::Plan(tom::Volume<std::complex<T> > &vsrc, tom::Volume<T> &vdst, unsigned flags) {

    const int type = this->is_valid_dft_c2r(vsrc, vdst);
    if (type == TOM_FFTW_PLAN_INVALID) {
        throw std::invalid_argument("tom::fftw::Plan() - The input volumes have not the right size or alignment for dft_c2r");
    }

    const int rank = type==TOM_FFTW_PLAN_3D ? 3 : 2;
    int howmany_rank = 0;
    fftw_iodim iodim[3];
    iodim[0].n  = vdst.getSizeX();
    iodim[0].is = vsrc.getStrideX()/sizeof(std::complex<T>);
    iodim[0].os = vdst.getStrideX()/sizeof(T);
    iodim[1].n  = vdst.getSizeY();
    iodim[1].is = vsrc.getStrideY()/sizeof(std::complex<T>);
    iodim[1].os = vdst.getStrideY()/sizeof(T);
    iodim[2].n  = vdst.getSizeZ();
    iodim[2].is = vsrc.getStrideZ()/sizeof(std::complex<T>);
    iodim[2].os = vdst.getStrideZ()/sizeof(T);
    fftw_iodim *howmany_dims = NULL;
    std::auto_ptr<tom::Volume<std::complex<T> > > auto_vol_fourier;
    std::auto_ptr<tom::Volume<T               > > auto_vol_real;


	boost::shared_ptr<void> p(static_cast<void *>(0), ::plan_destroyer<T>());
	if (tom::is_float<T>()) {
		p = boost::shared_ptr<void>(static_cast<void *>(fftwf_plan_guru_dft_c2r(rank, iodim, howmany_rank, howmany_dims, (fftwf_complex *)&vsrc.get(), (float  *)&vdst.get(), flags)), ::plan_destroyer<T>());
	} else {
		p = boost::shared_ptr<void>(static_cast<void *>(fftw_plan_guru_dft_c2r (rank, iodim, howmany_rank, howmany_dims, (fftw_complex  *)&vsrc.get(), (double *)&vdst.get(), flags)), ::plan_destroyer<T>());
	}
	if (!p.get()) {
		throw std::runtime_error("tom::fftw::Plan() - could not create fftw plan dft_c2r");
	}

	auto_vol_fourier.reset(new tom::Volume<std::complex<T> >(vsrc, NULL, vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ(), vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ()));
	auto_vol_real   .reset(new tom::Volume<T               >(vdst, NULL, vdst.getSizeX(), vdst.getSizeY(), vdst.getSizeZ(), vdst.getStrideX(), vdst.getStrideY(), vdst.getStrideZ()));

    this->plan.swap(p);
    this->type = tom::fftw::Plan<T>::FFTW_DFT_C2R;
    this->p_freq0 = auto_vol_fourier.release();
    this->p_time0 = auto_vol_real   .release();
	this->p_freq1 = NULL;
    this->p_time1 = NULL;
}





/****************************************************************************//**
 * \brief Execute the plan
 *
 * \param[in,out] vsrc The source volume.
 * \param[in,out] vdst The destination volume.
 *
 * It calles simply fftw_execute with the plan created upon instantiation.
 * The parameters are not acctually needed, its more for making clear what
 * is going to happen. i.e. you can not pass different volumes than at
 * construction time!
 *******************************************************************************/
template<typename T>
void tom::fftw::Plan<T>::execute(tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst) const {
	assert(this->type!=tom::fftw::Plan<T>::FFTW_DFT_R2C || (this->p_time0&&this->p_freq0));
    if (this->type != tom::fftw::Plan<T>::FFTW_DFT_R2C ||
        !this->p_freq0->equal_memory(vdst) ||
        !this->p_time0->equal_memory(vsrc)) {
        throw std::runtime_error("tom::Volume::execute (R2C) - Call execute with the volume upon construction of the plan.");
    }
    if (tom::is_float<T>()) {
        fftwf_execute((fftwf_plan)this->plan.get());
    } else {
        fftw_execute ((fftw_plan )this->plan.get());
    }
}

/****************************************************************************//**
 * \brief Execute the plan
 *
 * \param[in,out] vsrc The source volume.
 * \param[in,out] vdst The destination volume.
 *
 * It calles simply fftw_execute with the plan created upon instantiation.
 * The parameters are not acctually needed, its more for making clear what
 * is going to happen. i.e. you can not pass different volumes than at
 * construction time!
 *******************************************************************************/
template<typename T>
void tom::fftw::Plan<T>::execute(tom::Volume<std::complex<T> > &vsrc, tom::Volume<T> &vdst) const {
    if (this->type != tom::fftw::Plan<T>::FFTW_DFT_C2R ||
        !this->p_freq0->equal_memory(vsrc) ||
        !this->p_time0->equal_memory(vdst)) {
        throw std::runtime_error("tom::Volume::execute (C2R) - Call execute with the volume upon construction of the plan.");
    }
    if (tom::is_float<T>()) {
        fftwf_execute(reinterpret_cast<fftwf_plan>(this->plan.get()));
    } else {
        fftw_execute (reinterpret_cast<fftw_plan>(this->plan.get()));
    }
}






template class tom::fftw::Plan<float >;
template class tom::fftw::Plan<double>;













