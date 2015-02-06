/***********************************************************************//**
 * \file correlation_handler.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.12.2007
 **************************************************************************/
#include "tom/corr/correlation_handler.hpp"



#include "tom/core/volume_fcn.hpp"



/****************************************************************************//**
 * \brief Destructor of the class.
 *
 * Deletes the local mask.
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerPeaks<T>::~CorrelationHandlerPeaks() {
    delete this->mask;
}


/****************************************************************************//**
 * \brief Constructor of the class.
 *
 * \param[in] ntemplates Number of template to expect.
 * \param[in] number of particles to expect.
 *
 * The mask is not set. That means, every voxel of the correlation volume is
 * considered.
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerPeaks<T>::CorrelationHandlerPeaks(std::size_t ntemplates, std::size_t nparticles):
    list_peak(ntemplates),
    mask(NULL) {
    std::size_t itemplate;
    for (itemplate=0; itemplate<ntemplates; itemplate++) {
        list_peak.at(itemplate) = std::vector<tom::cc_peak<T> >(nparticles);
    }
}



/****************************************************************************//**
 * \brief Constructor of the class.
 *
 * \param[in] ntemplates Number of template to expect.
 * \param[in] number of particles to expect.
 * \param[in] mask Mask for the correlation volume to ignore certain voxels
 *   (for example those close to the end of the volume). A local copy of the
 *   input parameter is kept inside the object, thus the input parameter can be
 *   freed. \n As convention, the method \c process gets an non-fft-shifted
 *   correlation volume as input. But this input volume mask here, is expected to
 *   be fftshifted. Internally the constructor performs an inverse fftshift on its
 *   local copy of the mask. This is done this way, because it is more natural to
 *   pass for example a sphere mask to the constructor, but at the same time its
 *   non necessary to fftshift to volume before \c process. A voxel of mask, set to
 *   0 means to ignore peaks at the corresponding position in the correlation volume.
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerPeaks<T>::CorrelationHandlerPeaks(std::size_t ntemplates, std::size_t nparticles, const tom::Volume<T> &mask, bool is_fft_shifted):
    list_peak(ntemplates),
    mask(NULL) {
    std::size_t itemplate;
    for (itemplate=0; itemplate<ntemplates; itemplate++) {
        list_peak.at(itemplate) = std::vector<tom::cc_peak<T> >(nparticles);
    }

    std::auto_ptr<tom::Volume<T> > lmask(new tom::Volume<T>(mask.getSizeX(), mask.getSizeY(), mask.getSizeZ(), NULL, NULL));
    if (is_fft_shifted) {
        tom::fftshift(mask, *lmask, true);
    } else {
        lmask->setValues(mask);
    }
    this->mask = lmask.release();
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::vector<std::vector<tom::cc_peak<T> > > tom::CorrelationHandlerPeaks<T>::getPeakList() const {
    return this->list_peak;
}

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::CorrelationHandlerPeaks<T>::process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle) {


    std::vector<tom::st_idx> peaks;
    if (this->mask) {
       peaks  = tom::peak(v, *this->mask);
    } else {
       peaks  = tom::peak(v);
    }

    if (peaks.size() >= 1) {
        tom::cc_peak<T> peak = list_peak.at(itemplate).at(iparticle);
        const std::size_t &px = peaks[0].x;
        const std::size_t &py = peaks[0].y;
        const std::size_t &pz = peaks[0].z;
        peak.x = px;
        peak.y = py;
        peak.z = pz;
        T val = v.get(px, py, pz);
        if (val > peak.val || peak.angle_idx<0) {
            peak.val = val;
            peak.angle_idx = iangle;
            list_peak.at(itemplate).at(iparticle) = peak;
        }
    }
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerCC2EM<T>::CorrelationHandlerCC2EM(const std::vector<std::string> &dirnames, std::size_t ntemplates, std::size_t nparticles, const tom::Volume<double> &angles, bool do_fftshift):
    dirnames(dirnames),
    names(ntemplates*nparticles*angles.getSizeY(), std::string("")),
    success(ntemplates*nparticles*angles.getSizeY(), false),
    nparticles(nparticles),
    ntemplates(ntemplates),
    angles(angles),
    do_fftshift(do_fftshift) {
    if (ntemplates<1 || nparticles<1 || angles.getSizeX()!=3 || angles.getSizeZ()!=1) {
        throw std::invalid_argument("ntemplates and nparticles must be positive and angles must be a 3xNx1 volume.");
    }
    if (dirnames.size() != ntemplates*nparticles) {
        throw std::invalid_argument("The number of elements in the dirnames list must equal templates x particles.");
    }
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerCC2EM<T>::~CorrelationHandlerCC2EM() {
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::size_t tom::CorrelationHandlerCC2EM<T>::get_idx_all(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const {

    if (iangle>=this->angles.getSizeY()) {
        throw std::out_of_range("CorrelationHandlerCC2EM access out of range");
    }

    return this->get_idx_dirnames(itemplate, iparticle) * this->angles.getSizeY() + iangle;

}

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::size_t tom::CorrelationHandlerCC2EM<T>::get_idx_dirnames(std::size_t itemplate, std::size_t iparticle) const {

    if (itemplate>=this->ntemplates || iparticle>=this->nparticles) {
        throw std::out_of_range("CorrelationHandlerCC2EM access out of range");
    }

    return (itemplate * this->nparticles) + iparticle;

}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::CorrelationHandlerCC2EM<T>::getFileName(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const {

    const std::size_t idx = this->get_idx_dirnames(itemplate, iparticle);

    double a[3] = { this->angles.get(0, iangle, 0), this->angles.get(1, iangle, 0), this->angles.get(2, iangle, 0) };

    #define PI 3.141592653589793238512808959406186204433
    a[0] = round(a[0] * 180./PI * 100.) / 100.;
    a[1] = round(a[1] * 180./PI * 100.) / 100.;
    a[2] = round(a[2] * 180./PI * 100.) / 100.;

    std::string dirname = this->dirnames.at(idx);
    if (dirname == "") {
        dirname = ".";
    }

    std::stringstream ss;
    ss << dirname << (dirname[dirname.length()-1]=='/' ? "" : "/") << itemplate << "_" << iparticle << "_" << a[0] << 'x' << a[1] << 'x' << a[2] << ".em";

    return ss.str();

}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::pair<std::string, bool> tom::CorrelationHandlerCC2EM<T>::getFileStatus(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const {

    const std::size_t idx = this->get_idx_all(itemplate, iparticle, iangle);

    std::pair<std::string, bool> res;
    res.first = this->names.at(idx);
    res.second = this->success.at(idx);

    return res;

}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::CorrelationHandlerCC2EM<T>::process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle) {

    const std::size_t idx = this->get_idx_all(itemplate, iparticle, iangle);

    std::string n = this->getFileName(itemplate, iparticle, iangle);

    this->names.at(idx) = n;
    this->success.at(idx) = false;

    std::auto_ptr<tom::Volume<T> > v_shift;
    const tom::Volume<T> *pv = &v;
    if (this->do_fftshift) {
        v_shift.reset(new tom::Volume<T>(v.getSizeX(), v.getSizeY(), v.getSizeZ(), NULL, NULL));
        tom::fftshift(v, *v_shift, false);
        pv = v_shift.get();
    }

    try {
        pv->write_to_em(n, NULL);
    } catch (...) {
        return;
    }

    this->success.at(idx) = true;
}







/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerCC<T>::CorrelationHandlerCC(std::size_t ntemplates, std::size_t nparticles, std::size_t nangles, bool do_fftshift, std::size_t sizex, std::size_t sizey, std::size_t sizez):
    vols(ntemplates*nparticles*nangles, NULL),
    vols_local(ntemplates*nparticles*nangles, NULL),
    success(ntemplates*nparticles*nangles, false),
    nparticles(nparticles),
    ntemplates(ntemplates),
    nangles(nangles),
    sizex(sizex),
    sizey(sizey),
    sizez(sizez),
    do_fftshift(do_fftshift) {
    if (ntemplates<1 || nparticles<1 || nangles<1 || sizex<1 || sizey<1 || sizez<1) {
        throw std::invalid_argument("ntemplates, nparticles, nangles and the volume sizes must be positive. ");
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerCC<T>::~CorrelationHandlerCC() {
    const std::size_t l = this->vols_local.size();
    std::size_t i;
    for (i=0; i<l; i++) {
         delete this->vols_local[i];
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::size_t tom::CorrelationHandlerCC<T>::get_idx(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const {

    if (itemplate>=this->ntemplates || iparticle>=this->nparticles || iangle>=this->nangles) {
        throw std::out_of_range("CorrelationHandlerCC access out of range");
    }
    return ((itemplate * this->nparticles) + iparticle) * this->nangles + iangle;
}





/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::CorrelationHandlerCC<T>::reset(tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle) {

    if (v.getSizeX()!=this->sizex || v.getSizeY()!=this->sizey || v.getSizeZ()!=this->sizez) {
        throw std::invalid_argument("The input volume of CorrelationHandlerCC::reset has the wrong size.");
    }

    std::size_t idx = this->get_idx(itemplate, iparticle, iangle);

    delete this->vols_local[idx];
    this->vols_local[idx] = NULL;

    this->vols[idx] = &v;
    this->success[idx] = false;

}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
const tom::Volume<T> &tom::CorrelationHandlerCC<T>::get(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const {
    std::size_t idx = this->get_idx(itemplate, iparticle, iangle);
    tom::Volume<T> *v = this->vols[idx];
    if (!v) {
        throw std::runtime_error("The volume is still unset.");
    }
    return *v;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
bool tom::CorrelationHandlerCC<T>::get_success(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const {
    return this->success[this->get_idx(itemplate, iparticle, iangle)];
}





/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::CorrelationHandlerCC<T>::process(const tom::Volume<T> &volume, std::size_t itemplate, std::size_t iparticle, std::size_t iangle) {

    const std::size_t idx = this->get_idx(itemplate, iparticle, iangle);

    tom::Volume<T> *v = this->vols[idx];
    if (!v) {
        std::auto_ptr<tom::Volume<T> > autov(new tom::Volume<T>(this->sizex, this->sizey, this->sizez, NULL, NULL));
        this->vols[idx] = v = autov.get();
        this->vols_local[idx] = autov.release();
    }

    if (this->do_fftshift) {
        tom::fftshift(volume, *v, false);
    } else {
        v->setValues(volume);
    }

    this->success[idx] = true;
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::CorrelationHandlerPeakVolume<T>::process(const tom::Volume<T> &volume, std::size_t itemplate, std::size_t iparticle, std::size_t iangle) {

    if (iparticle >= this->nparticles || itemplate>=this->ntemplates) {
        throw std::out_of_range("CorrelationHandlerPeakVolume access out of range");
    }
    const std::size_t idx = itemplate*this->nparticles + iparticle;

    std::pair<boost::shared_ptr<tom::Volume<T> >, boost::shared_ptr<tom::Volume<int32_t> > > &v = this->vols[idx];
    if (!v.first.get() || !v.second.get()) {
        std::auto_ptr<tom::Volume<T      > > ccv(new tom::Volume<T      >(this->sizex, this->sizey, this->sizez, NULL,NULL));
        std::auto_ptr<tom::Volume<int32_t> > cci(new tom::Volume<int32_t>(this->sizex, this->sizey, this->sizez, NULL,NULL));
        ccv->setValues(volume);
        cci->setValues(iangle);
        v.first = ccv;
        v.second = cci;
    } else {
        tom::update_correlation_volume<T, int32_t>(volume, *v.first, *v.second, iangle);
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::CorrelationHandlerPeakVolume<T>::CorrelationHandlerPeakVolume(std::size_t ntemplates_, std::size_t nparticles_, std::size_t sizex_, std::size_t sizey_, std::size_t sizez_)
    :   CorrelationHandler<T>(),
        vols(ntemplates_*nparticles_),
        sizex(sizex_),
        sizey(sizey_),
        sizez(sizez_),
        ntemplates(ntemplates_),
        nparticles(nparticles_) {
}







// Template instantiations
template class tom::CorrelationHandler<float >;
template class tom::CorrelationHandler<double>;
template class tom::CorrelationHandlerPeaks<float >;
template class tom::CorrelationHandlerPeaks<double>;
template class tom::CorrelationHandlerPeakVolume<float >;
template class tom::CorrelationHandlerPeakVolume<double>;
template class tom::CorrelationHandlerCC<float >;
template class tom::CorrelationHandlerCC<double>;
template class tom::CorrelationHandlerCC2EM<float >;
template class tom::CorrelationHandlerCC2EM<double>;










