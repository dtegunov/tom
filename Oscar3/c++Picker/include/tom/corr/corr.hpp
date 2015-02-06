/****************************************************************************//**
 * \file corr.hpp
 * \brief Contains the correlation function corr3d.
 * \author  Thomas Haller
 * \version 0.1.2
 * \date    28.11.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__CORR_HPP__
#define ___INCLUDE_CORR__CORR_HPP__


#include <vector>
#include <cstddef>

#include "tom/core/wedge.hpp"
#include "tom/core/volume_container.hpp"
#include "tom/corr/correlation_handler.hpp"

#include "helper/triple.hpp"

namespace tom {



/****************************************************************************//**
 * Structure for a 3-tupel. It can be used for the 3-euler angles in rotation
 * or as 3D coordinate (vector).
 *******************************************************************************/
typedef helper::triple<double, double, double> double3;





template<typename T, typename T2>
void transf_particle(const tom::Volume<T> &src, tom::Volume<T> &dst, const tom::Volume<T> *mask, double phi, double psi, double theta, double shiftx, double shifty, double shiftz);


template<typename T, typename T2>
void average(   const std::vector<std::string> &fparticles,
                std::vector<double3> angles,
                std::vector<double3> shifts,
                const tom::Volume<T2> *mask,
                tom::Volume<T> &sum);


template <typename T, typename T2>
void align( const std::vector<std::string> &fparticles,
            const std::vector<std::string> &ftemplates,
            const tom::Volume<double> &vangles,
            const tom::Volume<T> *mask,
            const tom::Volume<T2> &mask_cc,
            const std::vector<tom::Wedge<T> *> wedge_particles,
            const std::vector<tom::Wedge<T> *> wedge_templates,
            std::vector<std::vector<tom::cc_peak<T> > > &list_peak, unsigned fftw_flags);
template <typename T, typename T2>
void corr3d(    std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                const std::vector<tom::VolumeContainer<T> *> &inparticles,
                std::size_t nparticles_amount,
                const std::vector<tom::VolumeContainer<T> *> &intemplates,
                std::size_t ntemplates_amount,
                const tom::Volume<double> &vangles,
                const tom::Volume<T2> *mask,
                const std::vector<tom::Wedge<T> *> *wedge_particles,
                const std::vector<tom::Wedge<T> *> *wedge_templates,
                tom::CorrelationHandler<T> *corrhdl,
                unsigned fftw_flags);


template<typename T, typename T2>
void corr3d_batch(      std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                        const std::vector<tom::VolumeContainer<T> *> &inparticles,
                        std::size_t nparticles_amount,
                        const std::vector<tom::VolumeContainer<T> *> &intemplates,
                        std::size_t ntemplates_amount,
                        const tom::Volume<double> &vangles,
                        const tom::Volume<T2> *mask,
                        const std::vector<tom::Wedge<T> *> *wedge_particles,
                        const std::vector<tom::Wedge<T> *> *wedge_templates,
                        tom::CorrelationHandler<T> *corrhdl,
                        unsigned fftw_flags);


template <typename T, typename T2>
void corr3d_single( const std::vector<const tom::Volume<T > *>               &templates,
                    const std::vector<const tom::Volume<std::complex<T> > *> &particles_fs,
                    const std::vector<double> &variance_particles,
                    const tom::Volume<double> &vangles,
                    std::size_t itemplates_baseidx,
                    std::size_t iparticles_baseidx,
                    const char *process_pair,
                    const tom::Volume<T2> *mask,
                    const std::vector<tom::Wedge<T> *> *wedge_templates,
                    const std::vector<tom::Wedge<T> *> *wedge_particles,
                    tom::CorrelationHandler<T> *corrhdl,
                    unsigned fftw_flags,
                    bool reapply_mask_after_wedge);


template<typename T>
bool average_wedges(const std::vector<tom::Wedge<T> *> &wedge_particles,
                    std::vector<double3> angles,
                    tom::Volume<T> &sum,
                    std::size_t sizez);


}


#endif



