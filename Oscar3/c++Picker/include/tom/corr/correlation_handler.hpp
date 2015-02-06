/****************************************************************************//**
 * \file correlation_handler.hpp
 * \brief Contains the class CorrelationHandler and its derived classes.
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.12.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__CORRELATION_HANDLER_HPP__
#define ___INCLUDE_CORR__CORRELATION_HANDLER_HPP__


#include <vector>

#include <boost/shared_ptr.hpp>

#include "tom/core/volume.hpp"
#include "tom/core/cc_peak.hpp"



namespace tom {




/****************************************************************************//**
 * \brief Class to process the correlation volume in tom::corr3d.
 *
 * This abstract base class has the virtual method process which is called in
 * tom::corr3d with each correlation volume. Derived classes process the data
 * accordingly.
 *******************************************************************************/
template<typename T>
class CorrelationHandler {
public:
    /// Virtual empty destructor of the class.
    virtual ~CorrelationHandler() {};

    /// abstract method process which is called in tom::corr3d.
    virtual void process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle) = 0;
};


/****************************************************************************//**
 * \brief Class to process the correlation volume in tom::corr3d.
 *
 * This process method of this class finds the peak in each correlation volume
 * and saves its maximum. The peak can be requested at the end.
 *******************************************************************************/
template<typename T>
class CorrelationHandlerPeaks: public CorrelationHandler<T> {
public:
    CorrelationHandlerPeaks(std::size_t ntemplates, std::size_t nparticles);
    CorrelationHandlerPeaks(std::size_t ntemplates, std::size_t nparticles, const tom::Volume<T> &mask, bool is_fft_shifted);
    virtual ~CorrelationHandlerPeaks();
    virtual void process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle);

    std::vector<std::vector<tom::cc_peak<T> > > getPeakList() const;

private:
    std::vector<std::vector<tom::cc_peak<T> > > list_peak; ///< The list of peaks. Each call or @c process may update it.
    tom::Volume<T> *mask;                                  ///< Saves the boolean mask.

};




template<typename T>
class CorrelationHandlerPeakVolume: public CorrelationHandler<T> {
public:
    CorrelationHandlerPeakVolume(std::size_t ntemplates, std::size_t nparticles, std::size_t sizex, std::size_t sizey, std::size_t sizez);
    virtual ~CorrelationHandlerPeakVolume() { };
    virtual void process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle);

    const std::vector<std::pair<boost::shared_ptr<tom::Volume<T> >, boost::shared_ptr<tom::Volume<int32_t> > > > &getVols() const { return this->vols; };
    std::vector<std::pair<boost::shared_ptr<tom::Volume<T> >, boost::shared_ptr<tom::Volume<int32_t> > > > getVols() { return this->vols; };

private:
    std::vector<std::pair<boost::shared_ptr<tom::Volume<T> >, boost::shared_ptr<tom::Volume<int32_t> > > > vols;
    std::size_t sizex, sizey, sizez;
    std::size_t ntemplates;
    std::size_t nparticles;
};




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
class CorrelationHandlerCC: public CorrelationHandler<T> {
public:
    CorrelationHandlerCC(std::size_t ntemplates, std::size_t nparticles, std::size_t nangles, bool do_fftshift, std::size_t sizex, std::size_t sizey, std::size_t sizez);
    virtual ~CorrelationHandlerCC();
    virtual void process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle);

    void reset(tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle);
    const tom::Volume<T> &get(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const;
    bool get_success(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const;


private:
    std::size_t get_idx(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const;

    std::vector<tom::Volume<T> *> vols;
    std::vector<tom::Volume<T> *> vols_local;
    std::vector<bool> success;


    std::size_t nparticles, ntemplates, nangles;
    std::size_t sizex, sizey, sizez;
    bool do_fftshift;
};




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
class CorrelationHandlerCC2EM: public CorrelationHandler<T> {
public:
    CorrelationHandlerCC2EM(const std::vector<std::string> &dirnames, std::size_t ntemplates, std::size_t nparticles, const tom::Volume<double> &angles, bool do_fftshift);
    virtual ~CorrelationHandlerCC2EM();
    virtual void process(const tom::Volume<T> &v, std::size_t itemplate, std::size_t iparticle, std::size_t iangle);

    std::pair<std::string, bool> getFileStatus(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const;
    std::string getFileName(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const;

private:
    std::size_t get_idx_all(std::size_t itemplate, std::size_t iparticle, std::size_t iangle) const;
    std::size_t get_idx_dirnames(std::size_t itemplate, std::size_t iparticle) const;

    std::vector<std::string> dirnames;
    std::vector<std::string> names;
    std::vector<bool> success;
    std::size_t nparticles, ntemplates;
    tom::Volume<double> angles;
    bool do_fftshift;
};



}


#endif

