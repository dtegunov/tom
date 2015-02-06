/****************************************************************************//**
 * \file average.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    06.02.2008
 *******************************************************************************/
#include <tom/corr/average.hpp>



#include <set>
#include <iomanip>
#include <boost/lexical_cast.hpp>



#include <helper/filesystem.hpp>

#include <tom/corr/config_files.hpp>
#include <tom/core/volume_fcn.hpp>






#ifndef THREAD_SAFE
#  error define THREAD_SAFE.
#endif
#if THREAD_SAFE
#  define __MAKE_STATIC__
#else
#  define __MAKE_STATIC__ static
#endif





/****************************************************************************//**
 * \brief Reads a particle from EM-file.
 *
 * \param[in] filename The name of the EM-file.
 * \param[in] volsize The expected size of the volume. If the size differs,
 *    and error is written to the log file and NULL is returned (or an
 *    exception is thrown depending on \c force_files_exist).
 * \param[in] iparticle The index of the particle. Only necessary for the
 *    logfile output.
 * \param[in] force_files_exist Boolean value. If true, and an error occures,
 *    the function throws an exception. Otherwise it returns NULL.
 * \param[in] clog Pointer to output stream where to write the log messages.
 *    If NULL, nothing is logged.
 * \return An \a auto_ptr of the particle. Before returning, the mean of
 *    the particle is set to zero.
 *******************************************************************************/
template<typename TFLOAT>
std::auto_ptr<tom::Volume<TFLOAT> > tom::avg::read_particle(const std::string &filename, std::size_t volsize, std::size_t iparticle, bool force_files_exist, std::ostream *clog) {

    assert(volsize>0);

    tom::Volume<TFLOAT> *pvol;
    std::auto_ptr<tom::Volume<TFLOAT> > vol;
    std::stringstream ss;

    // Read the particle.
    try {
        tom::read_from_em(pvol, filename, NULL,NULL,NULL, NULL, NULL,NULL);
        vol.reset(pvol);
    } catch (std::exception &e) {
        if (force_files_exist) {
            if (clog) { *clog << "# ERROR: could not read particle #" << iparticle << " (" << filename << ")." << std::endl; }
            throw;
        }
        if (clog) { *clog << "# WARNING: could not read particle #" << iparticle << " (" << filename << "). SKIP" << std::endl; }
        return std::auto_ptr<tom::Volume<TFLOAT> >();
    }

    if (!vol->is_equal_size(volsize)) {
        ss << "particle #" << iparticle << " (" << filename << ") has the wrong size.";
        if (force_files_exist) {
            if (clog) { *clog << "# ERROR: " << ss.str() << std::endl; }
            throw std::invalid_argument(ss.str());
        }
        if (clog) { *clog << "# WARNING: " << ss.str() << std::endl; }
        return std::auto_ptr<tom::Volume<TFLOAT> >();
    }

    double mean, variance;

    // Compute the mean and the variance of the particle.
    vol->stat(mean, variance, false);
    if (!variance) {
        // The variance is used, to abort in case of variance zero.
        ss.clear(); ss << "particle #" << iparticle << " (" << filename << ") has variance zero.";
        if (force_files_exist) {
            if (clog) { *clog << "# ERROR: " << ss.str() << std::endl; }
            throw std::invalid_argument(ss.str());
        }
        if (clog) { *clog << "# WARNING: " << ss.str() << std::endl; }
        return std::auto_ptr<tom::Volume<TFLOAT> >();
    }
    // The particle is set to mean==0 so that the rotation (which inserts zeros)
    // does not make a sharp edge at the border of the volume.
    vol->template shift_scale<TFLOAT>(-mean, 1);

    return vol;
}






/****************************************************************************//**
 * \brief Averages the particles.
 *
 * \param[in] volsize The expected size of (all) particles. Particles with
 *    different sizes are erroneous.
 * \param[in] fparticles The filename of all particles.
 * \param[in] wedge_particles The \c WedgeDescriptor for the wedge of each
 *    particles. This vector must have the same length as \a fparticles.
 * \param[in] idx Index-vector of the particles to use. Not all particles
 *    in \a fparticles are averaged. Only those, which index is in this vector.
 * \param[in] itemplate The number of the template, for which the average is
 *    computed.
 * \param[in] mask_sphere A pointer to the sphere mask applied to the volume.
 *    NULL means, no mask at all.
 * \param[in] peak_list The list of found peaks to know the shift of the volume.
 *    It has the size of the number of particles (<t>== fparticles.size()</t>)
 *    times the number of templates. It has the format as parsed by
 *    \a parse_peaklist.
 * \param[in] anglesv A vector with the rotation angles. Since the \a peak_list
 *    contains only the index of the angles, the best found angle is given in
 *    this field. Each element is [phi,psi,theta] in radians.
 * \param[in] force_files_exist boolean value. If true, and an error occures the
 *    function throws an exception.
 * \param[in] clog Pointer to output stream where to write log messages. If NULL,
 *    nothing is logged.
 * \param[out] avg Volume which contains the average after return.
 *    \a avg is not normalised (mean 0 and std 1) but all particles are, before
 *    averaging.
 * \param[out] avgwedge Average of the (rotated) wedges from the particles.
 *******************************************************************************/
template<typename TFLOAT>
void tom::avg::average(   std::size_t volsize,
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
                tom::Volume<TFLOAT> &avgwedge) {

    const std::size_t nparticles = fparticles.size();

    assert(volsize>0 &&
           nparticles*(peak_list.size() / nparticles)==peak_list.size() &&
           nparticles*(peak_list.size() / nparticles)==anglesv.size() &&
           itemplate<(peak_list.size() / nparticles) &&
           nparticles==wedge_particles.size() &&
           (!mask_sphere || mask_sphere->is_equal_size(volsize)) &&
           avg.is_equal_size(volsize) &&
           avgwedge.is_equal_size(volsize, volsize, volsize/2+1));

    std::size_t iparticle, i;

    double phi, psi, theta, shiftx, shifty, shiftz;
    double mean, variance;

    std::auto_ptr<tom::Volume<TFLOAT> > vol;
    std::auto_ptr<tom::Wedge<TFLOAT> > wedge;

    tom::Volume<TFLOAT> vol_tmp(volsize, volsize, volsize, NULL,NULL);


    std::stringstream ss;
    const std::size_t volcenter = static_cast<std::size_t>(ceil((volsize-1)/2.));

    std::size_t avgwedges_addval = 0;

    avg.setValues(0);
    avgwedge.setValues(0);

    assert(std::set<std::size_t>(idx.begin(), idx.end()).size() == idx.size());
    for (std::vector<std::size_t>::const_iterator it=idx.begin(); it!=idx.end(); it++) {
        iparticle = *it;
        i = itemplate*nparticles + iparticle;
        assert(iparticle < nparticles);

        phi = - anglesv[i].y;
        psi   = - anglesv[i].x;
        theta = - anglesv[i].z;
        shiftx = - peak_list[i].x;
        shifty = - peak_list[i].y;
        shiftz = - peak_list[i].z;

        assert(peak_list[i].angle_idx>=0 && fabs(shiftx)<=volsize/2 && fabs(shifty)<=volsize/2 && fabs(shiftz)<=volsize/2);


        // Read the particle.
        vol = read_particle<TFLOAT>(fparticles[iparticle], volsize, iparticle, force_files_exist, clog);
        if (!vol.get()) { continue; }


        // Load the wedge.
        if (wedge_particles[iparticle].get()) {
            wedge = wedge_particles[iparticle]->createWedge();
        } else {
            wedge.reset(NULL);
        }


        tom::rotate(*vol, vol_tmp, phi, psi, theta, volcenter, volcenter, volcenter, shiftx, shifty, shiftz, 0,0,0 );

        // Applay mask to particle and/or normalise
        if (mask_sphere) {
            tom::norm_mask<TFLOAT, TFLOAT, double>(vol_tmp, *mask_sphere, tom::norm::NORM_STDDEV_1, &variance, false);
        } else {
            vol_tmp.stat(mean, variance, false);
            if (variance) {
                vol_tmp.template shift_scale<TFLOAT>(-mean, 1/sqrt(variance));
            }
        }

        if (!variance) {
            if (clog) { *clog << "# WARNING: particle #" << iparticle << " (for template #" << itemplate << ") has variance zero after shifting and rotating. SKIP." << std::endl; }
            continue;
        }


        // Average the wedge...
        if (wedge.get() && wedge->is_active()) {
            wedge->rotate(phi, psi, theta);
            const tom::Volume<TFLOAT> *pwedge = wedge->get_wedge(volsize, volsize, volsize/2+1, volsize);
            if (pwedge) {
                tom::element_wise_add<TFLOAT, TFLOAT>(avgwedge, *pwedge);
            } else {
                avgwedges_addval ++;
            }
        } else {
            avgwedges_addval ++;
        }

        tom::element_wise_add(avg, vol_tmp);
    }
    avgwedge.template shift_scale<double>(avgwedges_addval, 1);

}



/****************************************************************************//**
 * \brief Scales every frequency of the average by its wedge.
 *
 * \param[in] avg0 The average of the particles as returned by \c average.
 * \param[in] wedge0_half The summed wedge. It must be in reduced complex form
 *    (i.e. the Z-dimension is only half) and NOT fftshifted.
 * \param[in] mask_sphere The sphere mask as used during \c average. While
 *    it was taken as a floating point weighting during \c average, it is now a
 *    boolean mask, to only set the voxels outside the mask to zero. NULL, means
 *    no mask.
 * \param[in] fftw_flag Flags for creating the fftw-plan.
 * \param[out] avg The resulting average.
 *
 * The average \a avg0 as computed by \c average results from particles
 * with missing frequencies in fourier space. Thus simply summing over these
 * particles (as happens in \c average) weigths some frequencies stronger then
 * others. If all particles would have no missing wedges, the summe as computed
 * by \c average needs to be divided by the number of particles. But because the
 * frequency space is not filled up, each frequency of \a avg0 is weighted
 * different (by dividing through the corresponding weight in \a wedge0_half).
 * \a avg0 and \a avg can be the same volume.
 *******************************************************************************/
template<typename T>
void tom::avg::apply_avgwedge(const tom::Volume<T> &avg0, const tom::Volume<T> &wedge0_half, const tom::Volume<T> *mask_sphere, unsigned fftw_flag, tom::Volume<T> &avg) {

    assert(wedge0_half.is_equal_size(avg0.getSizeX(), avg0.getSizeY(), avg0.getSizeZ()/2+1));

    __MAKE_STATIC__ std::auto_ptr<tom::Volume<T> > a_avg;
    __MAKE_STATIC__ std::auto_ptr<tom::Volume<std::complex<T> > > a_avg_ft;
    __MAKE_STATIC__ std::auto_ptr<tom::fftw::Plan<T> > a_plan_fwd;
    __MAKE_STATIC__ std::auto_ptr<tom::fftw::Plan<T> > a_plan_bwd;
    {
        bool change = false;
        if (!a_avg.get() || !a_avg->is_equal_size(avg0)) {
            a_avg.reset(new tom::Volume<T>(avg0.getSizeX(), avg0.getSizeY(), avg0.getSizeZ(), tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>()));
            change = true;
        }
        if (!a_avg_ft.get() || !a_avg_ft->is_equal_size(avg0.getSizeX(), avg0.getSizeY(), avg0.getSizeZ()/2+1)) {
            a_avg_ft.reset(new tom::Volume<std::complex<T> >(avg0.getSizeX(), avg0.getSizeY(), avg0.getSizeZ()/2+1, tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>()));
            change = true;
        }
        if (!a_plan_fwd.get() || change) {
            a_plan_fwd.reset(new tom::fftw::Plan<T>(*a_avg, *a_avg_ft, fftw_flag | FFTW_DESTROY_INPUT));
        }
        if (!a_plan_bwd.get() || change) {
            a_plan_bwd.reset(new tom::fftw::Plan<T>(*a_avg_ft, *a_avg, fftw_flag | FFTW_DESTROY_INPUT));
        }
    }

    a_avg->setValues(avg0);
    a_plan_fwd->execute(*a_avg, *a_avg_ft);
    tom::element_wise_div<std::complex<T>, T>(*a_avg_ft, wedge0_half, std::complex<T>(0,0));
    a_avg_ft->get(0,0,0) = std::complex<T>(0, 0);
    a_plan_bwd->execute(*a_avg_ft, *a_avg);


    if (mask_sphere) {
        tom::norm_mask<T, T, double>(*a_avg, *mask_sphere, tom::norm::NORM_STDDEV_1, NULL, true);
    } else {
        double v = a_avg->variance_mean_free(false);
        if (v) {
            a_avg->template shift_scale<T>(0, 1./sqrt(v));
        }
    }

    avg.setValues(*a_avg);
}









// template instatiations.
template void tom::avg::apply_avgwedge<float >(const tom::Volume<float > &avg0, const tom::Volume<float > &wedge0_half, const tom::Volume<float > *mask_sphere, unsigned fftw_flag, tom::Volume<float > &avg);
template void tom::avg::apply_avgwedge<double>(const tom::Volume<double> &avg0, const tom::Volume<double> &wedge0_half, const tom::Volume<double> *mask_sphere, unsigned fftw_flag, tom::Volume<double> &avg);


template void tom::avg::average<float >(std::size_t volsize, const std::vector<std::string> &fparticles,
                const std::vector<boost::shared_ptr<tom::WedgeDescriptor<float > > > &wedge_particles,
                const std::vector<std::size_t> &idx,
                std::size_t itemplate,
                const tom::Volume<float > *mask_sphere,
                const std::vector<tom::cc_peak<float > > &peak_list,
                const std::vector<helper::triple<double, double, double> > &anglesv,
                bool force_files_exist,
                std::ostream *clog,
                tom::Volume<float > &avg,
                tom::Volume<float > &avgwedge);
template void tom::avg::average<double>(std::size_t volsize, const std::vector<std::string> &fparticles,
                const std::vector<boost::shared_ptr<tom::WedgeDescriptor<double> > > &wedge_particles,
                const std::vector<std::size_t> &idx,
                std::size_t itemplate,
                const tom::Volume<double> *mask_sphere,
                const std::vector<tom::cc_peak<double> > &peak_list,
                const std::vector<helper::triple<double, double, double> > &anglesv,
                bool force_files_exist,
                std::ostream *clog,
                tom::Volume<double> &avg,
                tom::Volume<double> &avgwedge);


template std::auto_ptr<tom::Volume<float > > tom::avg::read_particle<float >(const std::string &filename, std::size_t volsize, std::size_t iparticle, bool force_files_exist, std::ostream *clog);
template std::auto_ptr<tom::Volume<double> > tom::avg::read_particle<double>(const std::string &filename, std::size_t volsize, std::size_t iparticle, bool force_files_exist, std::ostream *clog);





