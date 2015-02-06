/****************************************************************************//**
 * \file average.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    15.02.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__AVERAGE_HPP__
#define ___INCLUDE_CORR__AVERAGE_HPP__


#include <string>
#include <boost/shared_ptr.hpp>



#include <helper/triple.hpp>

#include <tom/core/volume.hpp>
#include <tom/core/wedge.hpp>
#include <tom/core/cc_peak.hpp>



namespace tom {
namespace avg {







template<typename T> void apply_avgwedge(const tom::Volume<T> &avg0, const tom::Volume<T> &wedge0_half, const tom::Volume<T> *mask_sphere, unsigned fftw_flag, tom::Volume<T> &avg);

template<typename TFLOAT>
void average(   std::size_t volsize,
                const std::vector<std::string> &fparticles,
                const std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > &wedge_particles,
                const std::vector<std::size_t> &idx,
                std::size_t itemplate,
                const tom::Volume<TFLOAT> *mask_sphere,
                const std::vector<tom::cc_peak<TFLOAT> > &peak_list,
                const std::vector<helper::triple<double, double, double> > &anglesv,
                bool force_files_exist,
                std::ostream *clog,
                tom::Volume<TFLOAT> &avg,
                tom::Volume<TFLOAT> &avgwedge);


template<typename TFLOAT> std::auto_ptr<tom::Volume<TFLOAT> > read_particle(const std::string &filename, std::size_t volsize, std::size_t iparticle, bool force_files_exist, std::ostream *clog);



} // namespace avg
} // namespace tom





#endif



