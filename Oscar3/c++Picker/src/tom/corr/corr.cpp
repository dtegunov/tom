/****************************************************************************//**
 * \file corr.cpp
 * \brief Contains implementations of the alignment functions.
 * \author  Thomas Haller
 * \version 0.1
 * \date    19.11.2007
 *******************************************************************************/
#include "tom/corr/corr.hpp"


#include <fftw3.h>
#include <iostream>
#include <string>
#include <sstream>
#include <libgen.h>
#include <typeinfo>

#define PI 3.141592653589793238512808959406186204433


#include "tom/core/io.h"
#include "tom/core/wedge.hpp"
#include "tom/core/volume_fcn.hpp"
#include "helper/auto_vector.hpp"
#include "helper/snippets.hpp"






/****************************************************************************//**
 * \brief Split a filename in basename, directory and extension.
 *******************************************************************************/
void split_filename(const std::string &filename, std::string &dirname_, std::string &basename_, std::string &extension) {

    std::vector<char> f2(filename.size() + 2);

    strcpy(&f2.at(0), filename.c_str());
    dirname_ = dirname(&f2.at(0));

    strcpy(&f2.at(0), filename.c_str());
    basename_ = basename(&f2.at(0));
    extension = "";

    std::size_t i = basename_.find_last_of('.');

    if (i != std::string::npos && i+1<basename_.size()) {
        extension = basename_.substr(i+1);
        basename_.resize(i);
    }
}


#if 0
#define DBGWRITE_BASE "/fs/home/haller/haller/DA/data/vols/run_"
#define DBGWRITE(x) if (true) { x } else (void)0
#else
#define DBGWRITE(x)
#endif


/****************************************************************************//**
 * \brief Average wedges.
 *******************************************************************************/
template<typename T>
bool tom::average_wedges(const std::vector<tom::Wedge<T> *> &wedge_particles,
                         std::vector<double3> angles,
                         tom::Volume<T> &sum,
                         std::size_t sizez) {

    const std::size_t nwedges = wedge_particles.size();
    std::size_t i;

    if (angles.size() != nwedges) {
        throw std::invalid_argument("The number of wedges and angles differs.");
    }

    bool reduced = false;
    if (sizez != sum.getSizeZ()) {
        if (sizez/2+1 != sum.getSizeZ()) {
            throw std::invalid_argument("either sizez or sizez/2+1 must be equal to the size of the wedge-volume");
        }
        reduced = true;
    }
    bool res = false;

    tom::Volume<T> sum_reduced(sum, NULL, sum.getSizeX(), sum.getSizeY(), sizez/2+1, sum.getStrideX(), sum.getStrideY(), sum.getStrideZ());

    const tom::Volume<T> *pwedge;
    std::size_t add_value = 0;
    for (i=0; i<nwedges; i++) {
        if (wedge_particles.at(i) && wedge_particles.at(i)->is_active()) {
            const double3 &a = angles[i];
            wedge_particles.at(i)->rotate(a.first, a.second, a.third);
            pwedge = wedge_particles.at(i)->get_wedge(sum_reduced.getSizeX(), sum_reduced.getSizeY(), sum_reduced.getSizeZ(), sizez);
        } else {
            pwedge = NULL;
        }

        if (pwedge) {
            if (!res) {
                sum_reduced.setValues(*pwedge);
                res = true;
            } else {
                tom::element_wise_add<T, T>(sum_reduced, *pwedge);
            }
        } else {
            add_value++;
        }
    }
    if (!res) {
        sum_reduced.setValues(add_value);
    } else {
        sum_reduced.template shift_scale<T>(add_value, 1.);
    }

    if (!reduced) {
        ::tom::hermitian_symmetry_to_full(sum);
    }

    return res;
}


/****************************************************************************//**
 * \brief Average paricles.
 *******************************************************************************/
template<typename T, typename T2>
void tom::average(  const std::vector<std::string> &fparticles,
                    std::vector<double3> angles,
                    std::vector<double3> shifts,
                    const tom::Volume<T2> *mask,
                    tom::Volume<T> &sum) {


    const std::size_t dims[3] = { sum.getSizeX(), sum.getSizeY(), sum.getSizeZ() };
    const std::size_t nparticles = fparticles.size();
    std::size_t iparticle;

    if ((mask && (dims[0]!=mask->getSizeX() || dims[1]!=mask->getSizeY() || dims[2]!=mask->getSizeZ()))) {
        throw std::invalid_argument("All volumes must have the same size.");
    }
    if (angles.size() != nparticles ||
        shifts.size() != nparticles) {
        throw std::invalid_argument("The number of particles and its angles and its wedges must correspond.");
    }

    sum.setValues(0.);

    std::auto_ptr<tom::Volume<T> > particle;
    tom::Volume<T> particle_rot(dims[0], dims[1], dims[2], NULL, NULL);

    double phi, psi, theta, shiftx, shifty, shiftz;

    for (iparticle=0; iparticle<nparticles; iparticle++) {
        tom::Volume<T> *pvol;
        tom::read_from_em<T>(pvol, fparticles.at(iparticle), NULL, NULL, NULL, NULL, NULL, NULL);
        particle.reset(pvol);

        if (dims[0]!=particle->getSizeX() || dims[1]!=particle->getSizeY() || dims[2]!=particle->getSizeZ()) {
            throw std::invalid_argument("All volumes must have the same dimension. ");
        }

        //particle->write_to_em("/fs/home/haller/haller/DA/data/vols/run_particle.em", NULL);

        phi = angles.at(iparticle).first;
        psi = angles.at(iparticle).second;
        theta = angles.at(iparticle).third;
        shiftx = shifts.at(iparticle).x;
        shifty = shifts.at(iparticle).y;
        shiftz = shifts.at(iparticle).z;
        tom::transf_particle<T, T2>(*particle, particle_rot, mask, phi, psi, theta, shiftx, shifty ,shiftz);

        //particle_rot.write_to_em("/fs/home/haller/haller/DA/data/vols/run_particle_rot.em", NULL);

        tom::element_wise_add<T,T>(sum, particle_rot);
    }



}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T, typename T2>
void tom::transf_particle(const tom::Volume<T> &src, tom::Volume<T> &dst, const tom::Volume<T> *mask, double phi, double psi, double theta, double shiftx, double shifty, double shiftz) {

    const std::size_t dims[3] = { src.getSizeX(), src.getSizeY(), src.getSizeZ() };
    if (dims[0]!=dst.getSizeX() || dims[1]!=dst.getSizeY() || dims[2]!=dst.getSizeZ() ||
        (mask && (dims[0]!=mask->getSizeX() || dims[1]!=mask->getSizeY() || dims[2]!=mask->getSizeZ()))) {
        throw std::invalid_argument("The input, output and mask volume must all have the same size.");
    }

    tom::Volume<T> v(src);

    double mean, variance;
    mean = v.mean();
    v.shift_scale(-mean, 1.);

    tom::rotate(v, dst, phi, psi, theta, ceil((dims[0]-1)/2.), ceil((dims[1]-1)/2.), ceil((dims[2]-1)/2.), shiftx, shifty, shiftz, 0,0,0 );

    if (mask) {
        tom::norm_mask<T, T, double>(dst, *mask, tom::norm::NORM_STDDEV_1, NULL, false);
    } else {
        v.stat(mean, variance, false);
        if (!variance) {
            variance = 1;
        } else {
            variance = 1./sqrt(variance);
        }
        dst.template shift_scale<double>(-mean, variance);
    }
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T, typename T2>
void tom::align(  const std::vector<std::string> &fparticles,
            const std::vector<std::string> &ftemplates,
            const tom::Volume<double> &vangles,
            const tom::Volume<T> *mask,
            const tom::Volume<T2> &mask_cc,
            const std::vector<tom::Wedge<T> *> wedge_particles,
            const std::vector<tom::Wedge<T> *> wedge_templates,
            std::vector<std::vector<tom::cc_peak<T> > > &list_peak, unsigned fftw_flags) {

    const size_t ntemplates = ftemplates.size();
    const size_t nparticles = fparticles.size();
    const size_t nangles = vangles.getSizeY();
    int itemplate, iparticle, iangle;

    std::size_t dims[3] = { mask_cc.getSizeX(), mask_cc.getSizeY(), mask_cc.getSizeZ() };
    std::size_t numel = dims[0]*dims[1]*dims[2];


    tom_io_em_header header;


    {
        // Initialise output vector.
        tom::st_idx idx_peak;
        idx_peak.x = 0; idx_peak.y = 0; idx_peak.z = 0;
        list_peak.resize(ntemplates);
        for (itemplate=0; (size_t)itemplate<ntemplates; itemplate++) {
            list_peak.at(itemplate) = std::vector<tom::cc_peak<T> >(nparticles);
        }
    }

    DBGWRITE( vangles.write_to_em(DBGWRITE_BASE "angles.em", NULL); );

    double variance_particle, variance_particle_orig, variance_template, mean;


    // Set memory allocation functions.
    if (tom::is_double<T>()) {
        tom::setFcnMem(fftw_malloc, fftw_free);
    } else if (tom::is_float<T>()) {
        tom::setFcnMem(fftwf_malloc, fftwf_free);
    } else {
        throw std::runtime_error("Template function is only for single or double floating point");
    }

    // check input parameters.
    if (ntemplates != wedge_templates.size()) { throw std::invalid_argument("The number of templates and its wedges must correspond."); }
    if (nparticles != wedge_particles.size()) { throw std::invalid_argument("The number of particles and its wedges must correspond."); }
    if (vangles.getSizeX() != 3) { throw std::invalid_argument("Angles must be a Nx3 matrix containing all angles."); }
    if (mask && (mask->getSizeX()!=dims[0] || mask->getSizeY()!=dims[1] || mask->getSizeZ()!=dims[2])) {
        throw std::invalid_argument("Angles must be a Nx3 matrix containing all angles.");
    }

    std::cout << "Try for " << nangles << " angles" << std::endl;
    vangles.printInfo("ANGLES");


    // Read all templates...
    auto_vector<tom::Volume<T> > templates(ntemplates);
    for (itemplate=0; (size_t)itemplate<ntemplates; itemplate++) {
        std::cout << (itemplate+1) << ": \"" << ftemplates.at(itemplate) << "\" - ";
        try {
            tom::Volume<T> *pvol;
            tom::read_from_em<T>(pvol, ftemplates.at(itemplate), NULL, NULL, NULL, &header, NULL, NULL);
            templates.assign_direct(itemplate, pvol);
        } catch (int &i) {
            std::cout << "Error reading file (" << i << "). Skip." << std::endl;
            continue;
        } catch (std::exception &e) {
            std::cout << "Error reading file (" << e.what() << "). Skip." << std::endl;
            continue;
        }
        std::cout << "File read: " << (int)header.machine << "," << (int)header.type << " - [" <<
            header.dims[0] << "," << header.dims[1] << "," << header.dims[2] << "] - [" <<
            templates.at(itemplate)->getSizeX() << "," << templates.at(itemplate)->getSizeY() << "," << templates.at(itemplate)->getSizeZ() << "]" << std::endl;
        templates.at(itemplate)->printInfo("template");

        if (dims[0]!=templates.at(itemplate)->getSizeX() || dims[1]!=templates.at(itemplate)->getSizeY() || dims[2]!=templates.at(itemplate)->getSizeZ()) {
            throw std::invalid_argument("All volumes must have the same dimension. ");
        }

        {
            templates.at(itemplate)->stat(mean, variance_template, false);
            if (!variance_template) {
                std::cout << "  HAS VARIANCE == 0... skip." << std::endl;
                templates.assign_direct(itemplate, NULL);
                continue;
            }
            templates.at(itemplate)->shift_scale(-mean, 1./ sqrt(variance_template));
        }

        DBGWRITE( templates.at(itemplate)->write_to_em(DBGWRITE_BASE "template_"+stringify(itemplate)+".em", NULL); );

    }

    std::cout << "------------------------------------------" << std::endl << std::endl;

    if (mask) {
        std::cout << "Used mask: " << std::endl;
        mask->printInfo("mask");

        DBGWRITE( mask->write_to_em(DBGWRITE_BASE "mask.em", NULL); );

    } else {
        std::cout << "Used mask: NONE" << std::endl;
    }

    tom::Volume<T2> mask_cc_ifftshift(dims[0],dims[1],dims[2],NULL,NULL);
    tom::fftshift(mask_cc, mask_cc_ifftshift, true);


    tom::Volume<T> ccf_max(dims[0], dims[1], dims[2], NULL, NULL);
    tom::Volume<int> ccf_max_angles(dims[0], dims[1], dims[2], NULL, NULL);

    tom::Volume<T> particle(dims[0], dims[1], dims[2], NULL, NULL);
    tom::Volume<T> particle2(dims[0], dims[1], dims[2], NULL, NULL);
    tom::Volume<T> template_rot(dims[0], dims[1], dims[2], NULL, NULL);
    tom::Volume<std::complex<T> > particle_fs(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > template_rot_fs(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > template_rot_fs2(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > particle_fs_wedge(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > particle_fs_wedge2(dims[0], dims[1], dims[2]/2+1, NULL, NULL);


    tom::fftw::Plan<T> plan_template_r2c(template_rot, template_rot_fs, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_template_c2r(template_rot_fs, template_rot, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_template_c2r2(template_rot_fs2, template_rot, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_particle_r2c(particle, particle_fs, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_particle_c2r2(particle_fs_wedge2, particle2, FFTW_DESTROY_INPUT | fftw_flags);





    for (iparticle=0; (size_t)iparticle<nparticles; iparticle++) {

        std::cout << "------------------------------------------" << std::endl << std::endl;
        std::cout << "Processing particle #" << (iparticle) << ": \"" << fparticles.at(iparticle) << "\" - ";

        {
            // Read particle...
            std::auto_ptr<tom::Volume<T> > auto_pvol;
            try {
                tom::Volume<T> *pvol;
                tom::read_from_em<T>(pvol, fparticles.at(iparticle), NULL, NULL, NULL, &header, NULL, NULL);
                auto_pvol.reset(pvol);
            } catch (int &i) {
                std::cout << "Error reading file (" << i << "). Skip." << std::endl;
                continue;
            } catch (std::exception &e) {
                std::cout << "Error reading file (" << e.what() << "). Skip." << std::endl;
                continue;
            }
            std::cout << "File read: " << (int)header.machine << "," << (int)header.type << " - ["
                        << header.dims[0] << "," << header.dims[1] << "," << header.dims[2] << "] - ["
                        << auto_pvol->getSizeX() << "," << auto_pvol->getSizeY() << "," << auto_pvol->getSizeZ() << "]" << std::endl;
            auto_pvol->printInfo("particle");

            if (dims[0] && dims[0]!=auto_pvol->getSizeX() || dims[1]!=auto_pvol->getSizeY() || dims[2]!=auto_pvol->getSizeZ()) {
                throw std::invalid_argument("All volumes must have the same dimension. ");
            }
            particle.setValues(*auto_pvol);
        }

        DBGWRITE( particle.write_to_em(DBGWRITE_BASE "particle_"+stringify(iparticle)+".em", NULL); );


        // Applay mask to particle.
        if (mask) {
            tom::norm_mask<T, T, double>(particle, *mask, tom::norm::NORM_NO_NORM, &variance_particle_orig, false);
        } else {
            particle.stat(mean, variance_particle_orig, false);
            particle.template shift_scale<double>(-mean, 1);
        }
        if (!variance_particle_orig) {
            std::cout << "  HAS VARIANCE == 0... skip." << std::endl;
            continue;
        }


        DBGWRITE( particle.write_to_em(DBGWRITE_BASE "particle_masked_"+stringify(iparticle)+".em", NULL); );


        // FFT particle.
        plan_particle_r2c.execute(particle, particle_fs);
        particle_fs.get(0,0,0) = 0;

        DBGWRITE( {
            tom::Volume<T>(particle_fs,true ).write_to_em(DBGWRITE_BASE "particle_fs_"+stringify(iparticle)+"_re.em", NULL);
            tom::Volume<T>(particle_fs,false).write_to_em(DBGWRITE_BASE "particle_fs_"+stringify(iparticle)+"_im.em", NULL);
        });

        for (itemplate=0; (size_t)itemplate<ntemplates; itemplate++) {
            std::cout << "with template #" << (itemplate+1) << ": \"" << ftemplates.at(itemplate) << "\"";
            if (!templates.at(itemplate)) {
                printf(" - not loaded. Skip\n");
                continue;
            }
            std::cout << std::endl;


            bool first_pass = true;
            for (iangle=0; iangle<(int)nangles; iangle++) {
                const double angles[3] = { vangles.get(0,iangle,0), vangles.get(1,iangle,0), vangles.get(2,iangle,0) };
                //const double angles[3] = { 1,1.5,0.75 };

                //std::cout << "    " << iangle << ": [" << (angles[0]*(180./PI)) << "," << (angles[1]*(180./PI)) << "," << (angles[2]*(180./PI)) << "]" << std::endl;
                //std::cout << "."; // << std::flush;
                //std::cout << "." << std::flush;

                // Applay wedge of the template to the particle
                variance_particle = variance_particle_orig;
                if (wedge_templates.at(itemplate) && wedge_templates.at(itemplate)->is_active()) {
                    particle_fs_wedge.setValues(particle_fs); // Initialise the local copy of the particle in fourier space.

                    wedge_templates.at(itemplate)->rotate(angles[0], angles[1], angles[2]); // Rotate the wedge of the template.


                    const bool changed = wedge_templates.at(itemplate)->apply(particle_fs_wedge, dims[2]);

                    DBGWRITE( throw std::runtime_error("write DBGWRITE for wedge of template"); );

                    particle_fs_wedge.get(0,0,0) = 0; // Set the mean to zero in fourier space.
                    if (changed) {
                        // the wedge changed the volume and now its no longer normalised.
                        particle_fs_wedge2.setValues(particle_fs_wedge);
                        plan_particle_c2r2.execute(particle_fs_wedge2, particle2);
                        variance_particle = particle2.variance_mean_free(false) / ((double)numel*numel);
                    }
                } else if (first_pass) {
                    particle_fs_wedge.setValues(particle_fs);
                }

                DBGWRITE({
                    tom::Volume<T>(particle_fs_wedge,true ).write_to_em(DBGWRITE_BASE "particle_fs_wedge_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_re.em", NULL);
                    tom::Volume<T>(particle_fs_wedge,false).write_to_em(DBGWRITE_BASE "particle_fs_wedge_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_im.em", NULL);
                });

                // Rotate the template.
                tom::rotate(*templates.at(itemplate), template_rot, angles[0], angles[1], angles[2]);

                DBGWRITE( template_rot.write_to_em(DBGWRITE_BASE "template_rot_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+".em", NULL); );


                // Applay mask to template.
                if (mask) {
                    tom::norm_mask<T, T, double>(template_rot, *mask, tom::norm::NORM_NO_NORM, &variance_template, false);
                } else {
                    template_rot.stat(mean, variance_template, false);
                    template_rot.template shift_scale<T>(-mean, 1.);
                }
                DBGWRITE( template_rot.write_to_em(DBGWRITE_BASE "template_rot_masked_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+".em", NULL); );


                // FFT of the template.
                plan_template_r2c.execute(template_rot, template_rot_fs);

                DBGWRITE( {
                    tom::Volume<T>(template_rot_fs,true ).write_to_em(DBGWRITE_BASE "template_fs_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_re.em", NULL);
                    tom::Volume<T>(template_rot_fs,false).write_to_em(DBGWRITE_BASE "template_fs_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_im.em", NULL);
                });


                if (wedge_particles.at(iparticle) && wedge_particles.at(iparticle)->is_active()) {

                    const bool changed = wedge_particles.at(iparticle)->apply(template_rot_fs, dims[2]);

                    DBGWRITE({
                        tom::Volume<T>(template_rot_fs,true ).write_to_em(DBGWRITE_BASE "template_fs_wedge0_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_re.em", NULL);
                        tom::Volume<T>(template_rot_fs,false).write_to_em(DBGWRITE_BASE "template_fs_wedge0_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_im.em", NULL);
                        tom::Volume<std::complex<T> > v(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
                        v.setValues(std::complex<T>(1,1));
                        wedge_particles.at(itemplate)->apply(v, dims[2]);
                        tom::Volume<T>(v,true ).write_to_em(DBGWRITE_BASE "template_fs_wedge1_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_re.em", NULL);
                        tom::Volume<T>(v,false).write_to_em(DBGWRITE_BASE "template_fs_wedge1_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_im.em", NULL);
                    });

                    template_rot_fs.get(0,0,0) = 0; // Set the mean to zero in fourier space.
                    if (changed) {
                        // the wedge changed the volume and now its no longer normalised.
                        template_rot_fs2.setValues(template_rot_fs);
                        plan_template_c2r2.execute(template_rot_fs2, template_rot);
                        variance_template = template_rot.variance_mean_free(false) / ((double)numel*numel);

                        DBGWRITE( template_rot.write_to_em(DBGWRITE_BASE "template_rot_wedge2_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+".em", NULL); );
                    }
                }

                DBGWRITE({
                    tom::Volume<T>(template_rot_fs,true ).write_to_em(DBGWRITE_BASE "template_rot_wedge3_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_re.em", NULL);
                    tom::Volume<T>(template_rot_fs,false).write_to_em(DBGWRITE_BASE "template_rot_wedge3_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_im.em", NULL);
                });


                tom::element_wise_conj_multiply<std::complex<T>, std::complex<T> >(template_rot_fs, particle_fs_wedge);


                DBGWRITE({
                    tom::Volume<T>(template_rot_fs,true ).write_to_em(DBGWRITE_BASE "template_v5_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_re.em", NULL);
                    tom::Volume<T>(template_rot_fs,false).write_to_em(DBGWRITE_BASE "template_v5_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+"_im.em", NULL);
                });

                // ifft
                plan_template_c2r.execute(template_rot_fs, template_rot);

                template_rot.template shift_scale<double>(0, 1./sqrt(variance_particle*variance_template));

                DBGWRITE( template_rot.write_to_em(DBGWRITE_BASE "template_v6_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+".em", NULL); );



                // Update the maximum corrlation matrix.
                if (first_pass) {
                    ccf_max_angles.setValues(iangle);
                    ccf_max.setValues(template_rot);
                    first_pass = false;
                } else {
                    tom::update_correlation_volume(template_rot, ccf_max, ccf_max_angles, iangle);
                }


                DBGWRITE({
                    if (iangle>50) {
                        iangle = nangles;
                    } else {
                        iangle += 5;
                    }
                });

            }


            ccf_max.template shift_scale<double>(0., 1. / ((double)numel*(double)numel));

            ccf_max.printInfo("ccf_max");

            ccf_max.write_to_em("/fs/home/haller/haller/DA/data/vols/run_1_ccf_max.em", NULL);
            template_rot.setValues(ccf_max_angles);
            template_rot.write_to_em("/fs/home/haller/haller/DA/data/vols/run_1_ccf_max_angles.em", NULL);


            {
                std::vector<tom::st_idx> peaks = tom::peak(ccf_max, mask_cc_ifftshift);
                tom::cc_peak<T> peak;
                if (peaks.size() >= 1) {
                    std::cout << "  peaks at" << std::string(peaks.size()>1?" (take the first)":"") << ": " << std::endl;
                    for (std::size_t i=0; i<peaks.size(); i++) {
                        int j = ccf_max_angles.get(peaks[i].x,peaks[i].y,peaks[i].z);
                        peak.x = peaks[i].x;
                        peak.y = peaks[i].y;
                        peak.z = peaks[i].z;
                        if (dims[0] - peak.x < (std::size_t)peak.x) { peak.x -= dims[0]; }
                        if (dims[1] - peak.y < (std::size_t)peak.y) { peak.y -= dims[1]; }
                        if (dims[2] - peak.z < (std::size_t)peak.z) { peak.z -= dims[2]; }
                        peak.x = -peak.x;
                        peak.y = -peak.y;
                        peak.z = -peak.z;
                        std::cout << "    " << i << ": [" << peaks[i].x << "," << peaks[i].y << "," << peaks[i].z << "] = [" << peak.x << "," << peak.y << "," << peak.z << "]: " << ccf_max.get(peaks[i].x,peaks[i].y,peaks[i].z) << " (" << j << " = [" << (vangles.get(0,j,0)*180./PI) << "," << (vangles.get(1,j,0)*180./PI) << "," << (vangles.get(2,j,0)*180./PI) << "]deg)" << std::endl;
                    }
                    const std::size_t &px = peaks[0].x;
                    const std::size_t &py = peaks[0].y;
                    const std::size_t &pz = peaks[0].z;
                    peak.x = px;
                    peak.y = py;
                    peak.z = pz;
                    peak.angle_idx = ccf_max_angles.get(px, py, pz);
                    peak.val = ccf_max.get(px,py,pz);
                    if (dims[0] - px < px) { peak.x -= dims[0]; }
                    if (dims[1] - py < py) { peak.y -= dims[1]; }
                    if (dims[2] - pz < pz) { peak.z -= dims[2]; }
                    peak.x = -peak.x;
                    peak.y = -peak.y;
                    peak.z = -peak.z;
                }
                list_peak.at(itemplate).at(iparticle) = peak;
            }

            #if 0
            {
                // save the ccf_max volumes...
                std::string t_dirname, t_basename, t_extension, p_dirname, p_basename, p_extension;
                ::tom::split_filename(ftemplates.at(itemplate), t_dirname, t_basename, t_extension);
                ::tom::split_filename(fparticles.at(iparticle), p_dirname, p_basename, p_extension);
                std::stringstream s; s << "/fs/home/haller/haller/DA/data/vols/run__" << p_basename << "_with_" << t_basename;
                std::string fname = s.str().append(".em");
                std::string fname_idx = s.str().append("_idx.em");
                std::cout << "  save ccf volume to \"" << fname << "\" and \"" << fname_idx << "\"" << std::endl;

                ccf_max.write_to_em(fname, NULL);
                tom::Volume<uint32_t>(ccf_max_angles).write_to_em(fname_idx, NULL);
            }
            #endif

            DBGWRITE( itemplate = ntemplates; );
        }

        DBGWRITE( iparticle = nparticles; );
    }

}



#ifdef NDEBUG
#   define CORR3D_VERBOSE_DEBUG(x) static_cast<void>(0)
#else
#   define CORR3D_VERBOSE_DEBUG(x) x
#endif

#undef CORR3D_VERBOSE_DEBUG
#define CORR3D_VERBOSE_DEBUG(x)
/****************************************************************************//**
 * \brief
 *******************************************************************************/
template <typename T, typename T2>
void tom::corr3d(   std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                    const std::vector<tom::VolumeContainer<T> *> &inparticles,
                    std::size_t nparticles_amount,
                    const std::vector<tom::VolumeContainer<T> *> &intemplates,
                    std::size_t ntemplates_amount,
                    const tom::Volume<double> &vangles,
                    const tom::Volume<T2> *mask,
                    const std::vector<tom::Wedge<T> *> *wedge_particles,
                    const std::vector<tom::Wedge<T> *> *wedge_templates,
                    tom::CorrelationHandler<T> *corrhdl,
                    unsigned fftw_flags) {

    const std::size_t ntemplates = intemplates.size();
    const std::size_t nparticles = inparticles.size();
    const std::size_t nangles = vangles.getSizeY();
    std::size_t itemplate, iparticle, iangle, itemplates_amount, iparticles_amount;
    std::size_t itemplate_all, iparticle_all;

    if (nparticles_amount < 1 || nparticles_amount > nparticles) { nparticles_amount = nparticles; }
    if (ntemplates_amount < 1 || ntemplates_amount > ntemplates) { ntemplates_amount = ntemplates; }

    const std::size_t dims[3] = { sizex__, sizey__, sizez__ };
    const std::size_t numel = dims[0]*dims[1]*dims[2];


    double variance_particle, variance_template, variance_template_orig, mean;


    // Set memory allocation functions.
    if (tom::is_double<T>()) {
        tom::setFcnMem(fftw_malloc, fftw_free);
    } else if (tom::is_float<T>()) {
        tom::setFcnMem(fftwf_malloc, fftwf_free);
    } else {
        throw std::runtime_error("Template function is only for single or double floating point");
    }

    // check input parameters.
    if (!corrhdl) { throw std::invalid_argument("No handle for the correlation volume given. Thus nothing is returned/processed."); }
    if (wedge_templates && ntemplates != wedge_templates->size()) { throw std::invalid_argument("The number of templates and its wedges must correspond."); }
    if (wedge_particles && nparticles!=wedge_particles->size()) { throw std::invalid_argument("The number of particles and its wedges must correspond."); }
    if (vangles.getSizeX() != 3) { throw std::invalid_argument("Angles must be a Nx3 matrix containing all angles."); }
    if (mask && (mask->getSizeX()!=dims[0] || mask->getSizeY()!=dims[1] || mask->getSizeZ()!=dims[2])) {
        throw std::invalid_argument("The size of the mask does not match");
    }

    CORR3D_VERBOSE_DEBUG( std::cout << "Try for " << nangles << " angles" << std::endl; );



    tom::Volume<T> voltmp(dims[0], dims[1], dims[2], NULL, NULL);
    tom::Volume<std::complex<T> > voltmp_fs(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > particle_fs_wedge(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > voltmp_fs2(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > *pparticle_fs_wedge;
    tom::Volume<std::complex<T> > template_fs_wedge(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<T> template_rot(dims[0], dims[1], dims[2], NULL, NULL);
    auto_vector<tom::Volume<std::complex<T> > > particles_fs(nparticles_amount);
    for (iparticle=0; iparticle<nparticles_amount; iparticle++) {
        particles_fs.assign_direct(iparticle, new tom::Volume<std::complex<T> >(dims[0], dims[1], dims[2]/2+1, NULL, NULL));
    }
    std::vector<bool> particles_loaded(nparticles_amount);

    auto_vector<tom::Volume<T> > templates(ntemplates_amount);
    for (itemplate=0; itemplate<ntemplates_amount; itemplate++) {
        templates.assign_direct(itemplate, new tom::Volume<T>(dims[0], dims[1], dims[2], NULL, NULL));
    }
    std::vector<bool> templates_loaded(ntemplates_amount);
    CORR3D_VERBOSE_DEBUG( std::vector<std::string> templates_loaded_errormsg(ntemplates_amount); );
    std::vector<double> variance_particles(nparticles_amount);

    tom::fftw::Plan<T> plan_r2c1(voltmp, voltmp_fs, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_c2r1(voltmp_fs2, voltmp, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_c2r2(template_fs_wedge, voltmp, FFTW_DESTROY_INPUT | fftw_flags);



    for (itemplates_amount=0; itemplates_amount<(ntemplates+ntemplates_amount-1)/ntemplates_amount; itemplates_amount++) {
        const std::size_t itemplates_amount_length= std::min<int>(ntemplates_amount, (int)ntemplates - (int)itemplates_amount*(int)ntemplates_amount);
        CORR3D_VERBOSE_DEBUG( std::cout << "process " << (itemplates_amount+1) << "th part of the " << ntemplates << " templates (" << (itemplates_amount*ntemplates_amount) << ".." << (itemplates_amount*ntemplates_amount+itemplates_amount_length-1) << ")" << std::endl; );

        ////////////////////////////////////////////////////////////////////
        // Get the current template-set
        ////////////////////////////////////////////////////////////////////
        for (itemplate=0, itemplate_all=itemplates_amount*ntemplates_amount; itemplate<itemplates_amount_length; itemplate++, itemplate_all++) {
            templates_loaded.at(itemplate) = false;
            CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = ""; );
            if (!intemplates.at(itemplate_all)) {
                CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = "not loaded"; );
                continue;
            }
            try {
                templates.at(itemplate)->setValues(intemplates.at(itemplate_all)->getVolume());
                intemplates.at(itemplate_all)->clearCache();
            } catch (std::exception &e) {
                CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = e.what(); );
                intemplates.at(itemplate_all)->clearCache();
                continue;
            }

            // Normalise to mean 0 and std 1.
            templates.at(itemplate)->stat(mean, variance_template_orig, false);
            if (!variance_template_orig) {
                CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = "has variance zero"; );
                continue;
            }
            templates.at(itemplate)->shift_scale(-mean, 1./ sqrt(variance_template_orig));
            templates_loaded.at(itemplate) = true;
        }

        for (iparticles_amount=0; iparticles_amount<(nparticles+nparticles_amount-1)/nparticles_amount; iparticles_amount++) {
            const std::size_t iparticles_amount_length = std::min<int>(nparticles_amount, (int)nparticles - (int)iparticles_amount*(int)nparticles_amount);

            CORR3D_VERBOSE_DEBUG( std::cout << "  process " << (iparticles_amount+1) << "th part of the " << nparticles << " particles (" << (iparticles_amount*nparticles_amount) << ".." << (iparticles_amount*nparticles_amount+iparticles_amount_length-1) << ")" << std::endl; );

            ////////////////////////////////////////////////////////////////////
            // Get the current particle set.
            ////////////////////////////////////////////////////////////////////
            for (iparticle=0, iparticle_all=iparticles_amount*nparticles_amount; iparticle<iparticles_amount_length; iparticle++, iparticle_all++) {
                CORR3D_VERBOSE_DEBUG( std::cout << "    particle " << iparticle_all; );
                particles_loaded.at(iparticle) = false;
                if (!inparticles.at(iparticle_all)) {
                    CORR3D_VERBOSE_DEBUG( std::cout << ": NULL" << std::endl; );
                    continue;
                }
                CORR3D_VERBOSE_DEBUG( std::cout << " (\"" << inparticles.at(iparticle_all)->getName() << "\")"; );
                tom::fftw::Plan<T> plan(voltmp, *particles_fs.at(iparticle), FFTW_DESTROY_INPUT | fftw_flags);
                try {
                    voltmp.setValues(inparticles.at(iparticle_all)->getVolume());
                    inparticles.at(iparticle_all)->clearCache();
                } catch (std::exception &e) {
                    CORR3D_VERBOSE_DEBUG( std::cout << ": ERROR: " << e.what() << std::endl; );
                    inparticles.at(iparticle_all)->clearCache();
                    continue;
                }

                // Applay mask to particle and/or normalise
                if (mask) {
                    tom::norm_mask<T, T2, double>(voltmp, *mask, tom::norm::NORM_NO_NORM, &variance_particles.at(iparticle), false);
                } else {
                    voltmp.stat(mean, variance_particles.at(iparticle), false);
                    voltmp.template shift_scale<double>(-mean, 1);
                }
                if (!variance_particles.at(iparticle)) {
                    CORR3D_VERBOSE_DEBUG( std::cout << ": ERROR: has zero variance" << std::endl; );
                    continue;
                }
                plan.execute(voltmp, *particles_fs.at(iparticle));
                particles_fs.at(iparticle)->get(0,0,0) = std::complex<T>(0,0);
                particles_loaded.at(iparticle) = true;
                CORR3D_VERBOSE_DEBUG( std::cout << std::endl; );
            }


            for (itemplate=0, itemplate_all=itemplates_amount*ntemplates_amount; itemplate<itemplates_amount_length; itemplate++, itemplate_all++) {

                CORR3D_VERBOSE_DEBUG( std::cout << "    template " << itemplate_all; );
                if (!templates_loaded.at(itemplate)) {
                    CORR3D_VERBOSE_DEBUG( std::cout << " Not loaded: " << templates_loaded_errormsg.at(itemplate) << ": SKIP" << std::endl; );
                    continue;
                }
                CORR3D_VERBOSE_DEBUG( std::cout << ": " << intemplates.at(itemplate_all)->getName() << std::endl; );
            }


            ////////////////////////////////////////////////////////////////////
            // Loop over all angles...
            ////////////////////////////////////////////////////////////////////
            for (iangle=0; iangle<nangles; iangle++) {
                const double angles[3] = { vangles.get(0,iangle,0), vangles.get(1,iangle,0), vangles.get(2,iangle,0) };
                //CORR3D_VERBOSE_DEBUG( std::cout << "."; );
                //CORR3D_VERBOSE_DEBUG( if (!(iangle%20)) { std::cout << std::flush; } );
                //CORR3D_VERBOSE_DEBUG( if (!(iangle%30)) { std::cout << (round(10000.*iangle/nangles)/100.) << "%" << std::endl; } );


                for (itemplate=0, itemplate_all=itemplates_amount*ntemplates_amount; itemplate<itemplates_amount_length; itemplate++, itemplate_all++) {
                    if (!templates_loaded.at(itemplate)) {
                        continue;
                    }

                    // Rotate the template.
                    tom::rotate(*templates.at(itemplate), voltmp, angles[0], angles[1], angles[2]);
                    // Applay mask to template.
                    if (mask) {
                        tom::norm_mask<T, T2, double>(voltmp, *mask, tom::norm::NORM_NO_NORM, &variance_template_orig, false);
                    } else {
                        voltmp.stat(mean, variance_template_orig, false);
                        voltmp.template shift_scale<T>(-mean, 1.);
                    }
                    // FFT of the template.
                    plan_r2c1.execute(voltmp, voltmp_fs);
                    voltmp_fs.get(0,0,0) = std::complex<T>(0,0);


                    for (iparticle=0, iparticle_all=iparticles_amount*nparticles_amount; iparticle<iparticles_amount_length; iparticle++, iparticle_all++) {

                        if (!particles_loaded.at(iparticle)) {
                            continue;
                        }

                        // Applay wedge of the template to the particle
                        variance_particle = variance_particles.at(iparticle);
                        if (wedge_templates && wedge_templates->at(itemplate_all) && wedge_templates->at(itemplate_all)->is_active()) {

                            wedge_templates->at(itemplate_all)->rotate(angles[0], angles[1], angles[2]); // Rotate the wedge of the template.

                            particle_fs_wedge.setValues(*particles_fs.at(iparticle)); // Initialise the local copy of the particle in fourier space.

                            const bool changed = wedge_templates->at(itemplate_all)->apply(particle_fs_wedge, dims[2]);

                            particle_fs_wedge.get(0,0,0) = std::complex<T>(0,0);
                            if (changed) {
                                // the wedge changed the volume and now its no longer normalised.
                                voltmp_fs2.setValues(particle_fs_wedge);
                                plan_c2r1.execute(voltmp_fs2, voltmp);
                                variance_particle = voltmp.variance_mean_free(false) / ((double)numel*numel);
                            }
                            pparticle_fs_wedge = &particle_fs_wedge;
                        } else {
                            pparticle_fs_wedge = particles_fs.at(iparticle);
                        }

                        // Applay wedge of the particle to the template
                        template_fs_wedge.setValues(voltmp_fs);
                        variance_template = variance_template_orig;
                        if (wedge_particles && wedge_particles->at(iparticle_all) && wedge_particles->at(iparticle_all)->is_active()) {

                            const bool changed = wedge_particles->at(iparticle_all)->apply(template_fs_wedge, dims[2]);

                            template_fs_wedge.get(0,0,0) = 0; // Set the mean to zero in fourier space.
                            if (changed) {
                                // the wedge changed the volume and now its no longer normalised.
                                voltmp_fs2.setValues(template_fs_wedge);
                                plan_c2r1.execute(voltmp_fs2, voltmp);
                                variance_template = voltmp.variance_mean_free(false) / ((double)numel*numel);

                                DBGWRITE( template_rot.write_to_em(DBGWRITE_BASE "template_rot_wedge2_"+stringify(iparticle)+"_"+stringify(itemplate)+"_"+stringify(iangle)+".em", NULL); );
                            }
                        }

                        #if 0
                        if (iangle==0) {
                            tom::Volume<std::complex<T> > vf(dims[0],dims[1],dims[2]/2+1,NULL,NULL);
                            tom::Volume<T> vt(dims[0],dims[1],dims[2],NULL,NULL);
                            tom::fftw::Plan<T> p(vf, vt, FFTW_DESTROY_INPUT | fftw_flags);

                            vf.setValues(template_fs_wedge);
                            p.execute(vf, vt);
                            vt.template shift_scale<T>(0, 1./numel/sqrt(variance_template));
                            vt.printInfo("V_TEMPLATE: " + stringify(itemplate));

                            vf.setValues(*pparticle_fs_wedge);
                            p.execute(vf, vt);
                            vt.template shift_scale<T>(0, 1./numel/sqrt(variance_particle));
                            vt.printInfo("V_PARTICLE: " + stringify(iparticle));
                        }
                        #endif


                        tom::element_wise_conj_multiply<std::complex<T>, std::complex<T> >(template_fs_wedge, *pparticle_fs_wedge);


                        // ifft
                        plan_c2r2.execute(template_fs_wedge, voltmp);

                        voltmp.template shift_scale<double>(0, 1./(sqrt(variance_particle*variance_template)*numel*numel));

                        corrhdl->process(voltmp, itemplate_all, iparticle_all, iangle);
                    }
                }
            }
            CORR3D_VERBOSE_DEBUG( std::cout << std::endl; );
        }
    }
}



/****************************************************************************//**
 * \brief Correlates a set of templates with a set of particles.
 *
 * \param[in] templates A vector with pointers to input templates (references).
 *   Some templates can be undefined by leaving its corresponding vector
 *   element NULL. The vector contains pointers to tom::Volumes with the
 *   templates in real space. The all must have the same size and data-type.\n
 *   The templates don't have to be preprocessed in any way, but it may be
 *   convenient to make it mean-free. At least in the case where the @a mask
 *   is not defined, because during rotation of the templates, zeros are rotated
 *   into the template. Thus (if it is not mean-free) there appears a edge in
 *   the volume.\n
 *   The volumes are not changed or overwritten in any way.
 * \param[in] particle_fs A vector of pointers to particles. Vector elements
 *   set to NULL are skipped. The particles are already in fourier space and must
 *   be preprocessed: If the mask is present, it must be applied by calling \c
 *   tom::norm_mask (no scaling for std-deviation needed). In any case
 *   it must have a mean of zero. Then the particle must be fouier transformed
 *   where the complex volume contains only the upper half of the Z-dimension.
 * \param[in] variance_particles A vector containing the variances of each particle
 *   as it had before fourier transformation. The values are not checked to be
 *   positive.
 * \param[in] vangles A 3xN volume with the rotation angles [phi psi theta] in
 *   radians.
 * \param[in] itemplates_baseidx A base index for the template, added before calling
 *   corrhdl->process.
 * \param[in] iparticles_baseidx A base index for the particle, added before calling
 *   corrhdl->process.
 * \param[in] mask An optional mask which is applied to each template after
 *   rotating it. It must be the same mask as already applied to the particles
 *   before fourier transforming them (as given in @a particle_fs). Can be NULL
 *   for no mask.
 * \param[in] wedge_templates A pointer to a vector of wedges of the templates.
 *   If it is NULL, no template has a wedge. Otherwise it must have the same length
 *   as @a templates. If a single wedge in the vector is NULL, that specific
 *   template has not wedge. Otherwise the wedge of the template is multiplied
 *   on each particle before correlation. If all templates have the same wedge,
 *   it may be convenient to let each element of @a wedge_templates point to the
 *   same object. It depends on the specific implementation of the wedge-class
 *   if this is possible. Probably it is faster because rotations can be
 *   computed only once in the wedge object.
 * \param[in] wedge_particles Analog to @a wedge_templates.
 * \param[out] corrhdl Pointer to the correlation handler. For each computed
 *   correlation of the template, particles and angles the method \c process is
 *   called (with adding @a itemplates_baseidx and @a iparticles_baseidx to the
 *   index parameters).
 * \param[in] fftw_flags Flags of the fftw-plan. All in all 3 fftw-plans are
 *   created.
 *
 *******************************************************************************/
template <typename T, typename T2>
void tom::corr3d_single( const std::vector<const tom::Volume<T > *>          &templates,
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
                    bool reapply_mask_after_wedge) {

    const std::size_t ntemplates = templates.size();
    const std::size_t nparticles = particles_fs.size();
    const std::size_t nangles = vangles.getSizeY();

    #define CORR_DO_RENORMALISATION_MASK 1

    std::size_t itemplate, iparticle, iangle;
    double variance_particle, variance_template, variance_template_orig, mean;

    std::size_t numel = 0;
    std::size_t dims[3] = {0,0,0};


    CORR3D_VERBOSE_DEBUG( std::cerr << __FILE__ << ":" << __LINE__ << ": " << (void *)mask << std::endl; )
    CORR3D_VERBOSE_DEBUG( if (mask) { mask->write_to_em("/fs/home/haller/haller/DA/data/outputdir/mask"+helper::stringify(0)+"_s.em", NULL); } )

    // check input parameters.
    if (mask) {
        numel = (dims[0] = mask->getSizeX()) * (dims[1] = mask->getSizeY()) * (dims[2] = mask->getSizeZ());
    } else {
        for (itemplate=0; itemplate<ntemplates; itemplate++) {
            if (templates[itemplate]) {
                const tom::Volume<T> &v = *templates[itemplate];
                numel = (dims[0] = v.getSizeX()) * (dims[1] = v.getSizeY()) * (dims[2] = v.getSizeZ());
                break;
            }
        }
    }
    if (!corrhdl) { throw std::invalid_argument("No handle for the correlation volume given. Thus nothing is returned/processed."); }
    if (wedge_templates && ntemplates!=wedge_templates->size()) { throw std::invalid_argument("The number of templates and its wedges must correspond."); }
    if (wedge_particles && nparticles!=wedge_particles->size()) { throw std::invalid_argument("The number of particles and its wedges must correspond."); }
    if (vangles.getSizeX() != 3 || vangles.getSizeZ() != 1) { throw std::invalid_argument("Angles must be a Nx3 matrix containing all angles."); }
    if (numel < 1) { return; }

    tom::Volume<T> voltmp(dims[0], dims[1], dims[2], NULL, NULL);
    tom::Volume<std::complex<T> > template_rot_fs(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > particle_fs_wedge(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > voltmp_fs2(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<std::complex<T> > template_fs_wedge(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
    tom::Volume<T> template_rot(dims[0], dims[1], dims[2], NULL, NULL);
    const tom::Volume<std::complex<T> > *pparticle_fs_wedge;

    tom::fftw::Plan<T> plan_r2c1(voltmp, template_rot_fs, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_r2c2(voltmp, particle_fs_wedge, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_r2c3(voltmp, template_fs_wedge, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_c2r1(voltmp_fs2, voltmp, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_c2r2(template_fs_wedge, voltmp, FFTW_DESTROY_INPUT | fftw_flags);
    tom::fftw::Plan<T> plan_c2r3(particle_fs_wedge, voltmp, FFTW_DESTROY_INPUT | fftw_flags);




    std::vector<char> process_pair_v;
    if (!process_pair) {
        process_pair_v.resize(ntemplates*nparticles, 1);
    }
    const char * const pprocess_pair = (process_pair ? process_pair : &process_pair_v[0]);

    CORR3D_VERBOSE_DEBUG( tom_io_em_header header; );

    reapply_mask_after_wedge = reapply_mask_after_wedge && mask;

    ////////////////////////////////////////////////////////////////////
    // Loop over all angles...
    ////////////////////////////////////////////////////////////////////
    for (iangle=0; iangle<nangles; iangle++) {
        const double angles[3] = { vangles.get(0,iangle,0), vangles.get(1,iangle,0), vangles.get(2,iangle,0) };
        CORR3D_VERBOSE_DEBUG( std::cout << "."; if (!(iangle%20)) { std::cout << std::flush; } );
        //CORR3D_VERBOSE_DEBUG( if (!(iangle%30)) { std::cout << (round(10000.*iangle/nangles)/100.) << "%" << std::endl; } );

        for (itemplate=0 ; itemplate<ntemplates; itemplate++) {
            if (!templates.at(itemplate)) { continue; }

            // Rotate the template.
            tom::rotate(*templates.at(itemplate), voltmp, angles[0], angles[1], angles[2]);

            CORR3D_VERBOSE_DEBUG(
                header.machine = 6;
                header.emdata[0] = static_cast<int32_t>(round(angles[0]*100000));
                header.emdata[1] = static_cast<int32_t>(round(angles[1]*100000));
                header.emdata[2] = static_cast<int32_t>(round(angles[2]*100000));
                voltmp.write_to_em("/fs/home/haller/haller/DA/data/outputdir/templ_rot_"+helper::stringify(itemplate)+"_"+helper::stringify(iangle)+".em", &header);
            );



            // Applay mask to template.
            if (mask) {
                tom::norm_mask<T, T2, double>(voltmp, *mask, tom::norm::NORM_NO_NORM, &variance_template_orig, false);
            } else {
                voltmp.stat(mean, variance_template_orig, false);
                voltmp.template shift_scale<T>(-mean, 1.);
            }

            CORR3D_VERBOSE_DEBUG( voltmp.write_to_em("/fs/home/haller/haller/DA/data/outputdir/templ_rot_norm_"+helper::stringify(itemplate)+"_"+helper::stringify(iangle)+".em", &header); );

            // FFT of the template.
            plan_r2c1.execute(voltmp, template_rot_fs);
            template_rot_fs.get(0,0,0) = std::complex<T>(0,0);

            for (iparticle=0; iparticle<nparticles; iparticle++) {
                if (!particles_fs.at(iparticle) || !pprocess_pair[itemplate*nparticles + iparticle]) { continue; }

                // Applay wedge of the template to the particle
                variance_particle = variance_particles.at(iparticle);
                if (wedge_templates && wedge_templates->at(itemplate) && wedge_templates->at(itemplate)->is_active()) {
                    wedge_templates->at(itemplate)->rotate(angles[0], angles[1], angles[2]); // Rotate the wedge of the template.

                    particle_fs_wedge.setValues(*particles_fs.at(iparticle)); // Initialise the local copy of the particle in fourier space.
                    const bool changed = wedge_templates->at(itemplate)->apply(particle_fs_wedge, dims[2]);

                    particle_fs_wedge.get(0,0,0) = std::complex<T>(0,0);
                    if (changed) {
                        // the wedge changed the volume and now its no longer normalised.
                        if (reapply_mask_after_wedge) {
                            plan_c2r3.execute(particle_fs_wedge, voltmp);
                            tom::norm_mask<T, T2, double>(voltmp, *mask, tom::norm::NORM_NO_NORM, &variance_particle, true);
                            plan_r2c2.execute(voltmp, particle_fs_wedge);
                        } else {
                            // Only recompute the variance. No renormalisation.
                            voltmp_fs2.setValues(particle_fs_wedge);
                            plan_c2r1.execute(voltmp_fs2, voltmp);
                            variance_particle = voltmp.variance_mean_free(false) / ((double)numel*numel);
                        }

                    }
                    pparticle_fs_wedge = &particle_fs_wedge;
                } else {
                    pparticle_fs_wedge = particles_fs.at(iparticle);
                }


                #if 0
                {
                    tom::Volume<std::complex<T> > v_f(dims[0], dims[1], dims[2]/2+1, NULL,NULL);
                    tom::Volume<T> v_t(dims[0], dims[1], dims[2], NULL,NULL);
                    tom::fftw::Plan<T> plan(v_f, v_t, FFTW_DESTROY_INPUT | fftw_flags);
                    v_f.setValues(*pparticle_fs_wedge);
                    plan.execute(v_f, v_t);
                    v_t.template shift_scale<double>(0, 1./numel);
                    v_t.write_to_em("/fs/home/haller/haller/DA/data/outputdir/particle_"+helper::stringify(itemplate)+"_"+helper::stringify(iparticle)+"_"+helper::stringify(iangle)+".em", NULL);
                }
                #endif

                #if 0
                {
                    tom::Volume<std::complex<T> > v_f(dims[0], dims[1], dims[2]/2+1, NULL,NULL);
                    tom::Volume<T> v_t(dims[0], dims[1], dims[2], NULL,NULL);
                    tom::fftw::Plan<T> plan(v_f, v_t, FFTW_DESTROY_INPUT | fftw_flags);
                    v_f.setValues(template_rot_fs);
                    plan.execute(v_f, v_t);
                    v_t.template shift_scale<double>(0, 1./numel);
                    v_t.write_to_em("/fs/home/haller/haller/DA/data/outputdir/template_nonwedged"+helper::stringify(itemplate)+"_"+helper::stringify(iparticle)+"_"+helper::stringify(iangle)+".em", NULL);
                }
                #endif

                //mask->write_to_em("/fs/home/haller/haller/DA/data/outputdir/mask.em", NULL);


                CORR3D_VERBOSE_DEBUG( tom::Volume<T>(template_rot_fs, true ).write_to_em("/fs/home/haller/haller/DA/data/outputdir/templ_fs_wedge_r_"+helper::stringify(itemplate)+"_"+helper::stringify(iangle)+".em", &header); );
                CORR3D_VERBOSE_DEBUG( tom::Volume<T>(template_rot_fs, false).write_to_em("/fs/home/haller/haller/DA/data/outputdir/templ_fs_wedge_c_"+helper::stringify(itemplate)+"_"+helper::stringify(iangle)+".em", &header); );


                // Applay wedge of the particle to the template
                variance_template = variance_template_orig;
                bool inplace_multiplication = false;
                if (wedge_particles && wedge_particles->at(iparticle) && wedge_particles->at(iparticle)->is_active()) {
                    CORR3D_VERBOSE_DEBUG( std::cerr << __FILE__ ":" << __LINE__ << ": wedge_particles" << iparticle << " " << itemplate << " " << iangle << ": " << typeid(*wedge_templates->at(itemplate)).name() << std::endl; );

                    template_fs_wedge.setValues(template_rot_fs);
                    const bool changed = wedge_particles->at(iparticle)->apply(template_fs_wedge, dims[2]);

                    template_fs_wedge.get(0,0,0) = 0; // Set the mean to zero in fourier space.
                    if (changed) {
                        // The application of the wedge changed the normalisation of the
                        // wedge.
                        if (reapply_mask_after_wedge) {
                            plan_c2r2.execute(template_fs_wedge, voltmp);
                            tom::norm_mask<T, T2, double>(voltmp, *mask, tom::norm::NORM_NO_NORM, &variance_template, true);
                            plan_r2c3.execute(voltmp, template_fs_wedge);
                        } else {
                            // Only recompute the variance. No renormalisation.
                            voltmp_fs2.setValues(template_fs_wedge);
                            plan_c2r1.execute(voltmp_fs2, voltmp);
                            variance_template = voltmp.variance_mean_free(false) / ((double)numel*numel);
                        }
                        inplace_multiplication = true;
                    }
                    #if 0
                    {
                        tom::Volume<std::complex<T> > v_f(dims[0], dims[1], dims[2]/2+1, NULL,NULL);
                        tom::Volume<T> v_t(dims[0], dims[1], dims[2], NULL,NULL);
                        tom::fftw::Plan<T> plan(v_f, v_t, FFTW_DESTROY_INPUT | fftw_flags);
                        if (inplace_multiplication) {
                            v_f.setValues(template_fs_wedge);
                        } else {
                            v_f.setValues(template_rot_fs);
                        }
                        plan.execute(v_f, v_t);
                        v_t.template shift_scale<double>(0, 1./numel);
                        v_t.write_to_em("/fs/home/haller/haller/DA/data/outputdir/template_wedged"+helper::stringify(itemplate)+"_"+helper::stringify(iparticle)+"_"+helper::stringify(iangle)+".em", NULL);
                    }
                    #endif
                }

                if (inplace_multiplication) {
                    tom::element_wise_conj_multiply<std::complex<T>, std::complex<T> >(template_fs_wedge, *pparticle_fs_wedge);
                } else {
                    tom::element_wise_conj_multiply<std::complex<T>, std::complex<T>, std::complex<T> >(template_fs_wedge, template_rot_fs, *pparticle_fs_wedge);
                }

                template_fs_wedge.get(0,0,0) = 0;


                // ifft
                plan_c2r2.execute(template_fs_wedge, voltmp);

                voltmp.template shift_scale<double>(0, 1./(sqrt(variance_particle*variance_template)*numel*numel));

                CORR3D_VERBOSE_DEBUG( std::cerr << __FILE__ ":" << __LINE__ << ": " <<  typeid(*corrhdl).name() << std::endl; );

                corrhdl->process(voltmp, itemplate+itemplates_baseidx, iparticle+iparticles_baseidx, iangle);
            }
        }
    }
    CORR3D_VERBOSE_DEBUG( std::cout << std::endl; );
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template <typename T, typename T2>
void tom::corr3d_batch( std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                        const std::vector<tom::VolumeContainer<T> *> &inparticles,
                        std::size_t nparticles_amount,
                        const std::vector<tom::VolumeContainer<T> *> &intemplates,
                        std::size_t ntemplates_amount,
                        const tom::Volume<double> &vangles,
                        const tom::Volume<T2> *mask,
                        const std::vector<tom::Wedge<T> *> *wedge_particles,
                        const std::vector<tom::Wedge<T> *> *wedge_templates,
                        tom::CorrelationHandler<T> *corrhdl,
                        unsigned fftw_flags) {

    const std::size_t ntemplates = intemplates.size();
    const std::size_t nparticles = inparticles.size();

    std::size_t itemplate, iparticle, itemplates_amount, iparticles_amount;
    std::size_t itemplate_all, iparticle_all;

    if (nparticles_amount < 1 || nparticles_amount > nparticles) { nparticles_amount = nparticles; }
    if (ntemplates_amount < 1 || ntemplates_amount > ntemplates) { ntemplates_amount = ntemplates; }

    const std::size_t dims[3] = { sizex__, sizey__, sizez__ };

    // Set memory allocation functions.
    tom::fftw::setFcnMemFFTW<T>();

    // check input parameters.
    if (!corrhdl) { throw std::invalid_argument("No handle for the correlation volume given. Thus nothing is returned/processed."); }
    if (wedge_templates && ntemplates != wedge_templates->size()) { throw std::invalid_argument("The number of templates and its wedges must correspond."); }
    if (wedge_particles && nparticles!=wedge_particles->size()) { throw std::invalid_argument("The number of particles and its wedges must correspond."); }
    if (vangles.getSizeX() != 3) { throw std::invalid_argument("Angles must be a Nx3 matrix containing all angles."); }
    if (mask && (mask->getSizeX()!=dims[0] || mask->getSizeY()!=dims[1] || mask->getSizeZ()!=dims[2])) {
        throw std::invalid_argument("The size of the mask does not match");
    }

    CORR3D_VERBOSE_DEBUG( std::cout << "Try for " << vangles.getSizeY() << " angles" << std::endl; );

    tom::Volume<T> voltmp(dims[0], dims[1], dims[2], NULL, NULL);

    auto_vector<tom::fftw::Plan<T> > plan(nparticles_amount);

    auto_vector<tom::Volume<std::complex<T> > > particles_fs(nparticles_amount);
    for (iparticle=0; iparticle<nparticles_amount; iparticle++) {
        particles_fs.assign_direct(iparticle, new tom::Volume<std::complex<T> >(dims[0], dims[1], dims[2]/2+1, NULL, NULL));
        plan.assign_direct(iparticle, new tom::fftw::Plan<T>(voltmp, *particles_fs.at(iparticle), FFTW_DESTROY_INPUT | fftw_flags));
    }
    auto_vector<tom::Volume<T> > templates(ntemplates_amount);
    for (itemplate=0; itemplate<ntemplates_amount; itemplate++) {
        templates.assign_direct(itemplate, new tom::Volume<T>(dims[0], dims[1], dims[2], NULL, NULL));
    }

    CORR3D_VERBOSE_DEBUG( std::vector<std::string> templates_loaded_errormsg(ntemplates_amount); );

    std::vector<double> variance_particles(nparticles_amount);


    std::vector<const tom::Volume<std::complex<T> > *> s_particle_fs(nparticles_amount);
    std::vector<const tom::Volume<T> *> s_templates(ntemplates_amount);


    double variance_template_orig, mean;

    const std::vector<tom::Wedge<T> *> *pwedge_templates = NULL;
    const std::vector<tom::Wedge<T> *> *pwedge_particles = NULL;
    std::auto_ptr<std::vector<tom::Wedge<T> *> > wedge_templates_auto;
    std::auto_ptr<std::vector<tom::Wedge<T> *> > wedge_particles_auto;


    for (itemplates_amount=0; itemplates_amount<(ntemplates+ntemplates_amount-1)/ntemplates_amount; itemplates_amount++) {
        const std::size_t itemplates_amount_length=std::min<int>(ntemplates_amount, (int)ntemplates - (int)itemplates_amount*(int)ntemplates_amount);
        s_templates.resize(itemplates_amount_length);
        CORR3D_VERBOSE_DEBUG( std::cout << "process " << (itemplates_amount+1) << "th part of the " << ntemplates << " templates (" << (itemplates_amount*ntemplates_amount) << ".." << (itemplates_amount*ntemplates_amount+itemplates_amount_length-1) << ")" << std::endl; );

        ////////////////////////////////////////////////////////////////////
        // Get the current template-set
        ////////////////////////////////////////////////////////////////////
        for (itemplate=0, itemplate_all=itemplates_amount*ntemplates_amount; itemplate<itemplates_amount_length; itemplate++, itemplate_all++) {
            s_templates.at(itemplate) = NULL;
            CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = ""; );
            if (!intemplates.at(itemplate_all)) {
                CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = "not loaded"; );
                continue;
            }
            try {
                templates.at(itemplate)->setValues(intemplates.at(itemplate_all)->getVolume());
                intemplates.at(itemplate_all)->clearCache();
            } catch (std::exception &e) {
                CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = e.what(); );
                intemplates.at(itemplate_all)->clearCache();
                continue;
            }

            // Normalise to mean 0 and std 1.
            templates.at(itemplate)->stat(mean, variance_template_orig, false);
            if (!variance_template_orig) {
                CORR3D_VERBOSE_DEBUG( templates_loaded_errormsg.at(itemplate) = "has variance zero"; );
                continue;
            }
            templates.at(itemplate)->template shift_scale<T>(-mean, 1);
            s_templates.at(itemplate) = templates.at(itemplate);
        }

        // Copy the references to the current wedge.
        if (wedge_templates) {
            wedge_templates_auto.reset(new std::vector<tom::Wedge<T> *>(itemplates_amount_length));
            for (std::size_t i=0; i<itemplates_amount_length; i++) {
                wedge_templates_auto.get()->at(i) = wedge_templates->at(itemplates_amount*ntemplates_amount+i);
            }
            pwedge_templates = wedge_templates_auto.get();
        }

        for (iparticles_amount=0; iparticles_amount<(nparticles+nparticles_amount-1)/nparticles_amount; iparticles_amount++) {
            const std::size_t iparticles_amount_length =std::min<int>(nparticles_amount, (int)nparticles - (int)iparticles_amount*(int)nparticles_amount);
            s_particle_fs.resize(iparticles_amount_length);

            CORR3D_VERBOSE_DEBUG( std::cout << "  process " << (iparticles_amount+1) << "th part of the " << nparticles << " particles (" << (iparticles_amount*nparticles_amount) << ".." << (iparticles_amount*nparticles_amount+iparticles_amount_length-1) << ")" << std::endl; );

            ////////////////////////////////////////////////////////////////////
            // Get the current particle set.
            ////////////////////////////////////////////////////////////////////
            for (iparticle=0, iparticle_all=iparticles_amount*nparticles_amount; iparticle<iparticles_amount_length; iparticle++, iparticle_all++) {
                CORR3D_VERBOSE_DEBUG( std::cout << "    particle " << iparticle_all; );
                s_particle_fs.at(iparticle) = NULL;
                if (!inparticles.at(iparticle_all)) {
                    CORR3D_VERBOSE_DEBUG( std::cout << ": NULL" << std::endl; );
                    continue;
                }
                CORR3D_VERBOSE_DEBUG( std::cout << " (\"" << inparticles.at(iparticle_all)->getName() << "\")"; );
                try {
                    voltmp.setValues(inparticles.at(iparticle_all)->getVolume());
                    inparticles.at(iparticle_all)->clearCache();
                } catch (std::exception &e) {
                    CORR3D_VERBOSE_DEBUG( std::cout << ": ERROR: " << e.what() << std::endl; );
                    inparticles.at(iparticle_all)->clearCache();
                    continue;
                }

                // Applay mask to particle and/or normalise
                if (mask) {
                    tom::norm_mask<T, T2, double>(voltmp, *mask, tom::norm::NORM_NO_NORM, &variance_particles.at(iparticle), false);
                } else {
                    voltmp.stat(mean, variance_particles.at(iparticle), false);
                    voltmp.template shift_scale<double>(-mean, 1);
                }
                if (!variance_particles.at(iparticle)) {
                    CORR3D_VERBOSE_DEBUG( std::cout << ": ERROR: has zero variance" << std::endl; );
                    continue;
                }
                plan.at(iparticle)->execute(voltmp, *particles_fs.at(iparticle));
                particles_fs.at(iparticle)->get(0,0,0) = std::complex<T>(0,0);
                s_particle_fs.at(iparticle) = particles_fs.at(iparticle);
                CORR3D_VERBOSE_DEBUG( std::cout << std::endl; );
            }


            CORR3D_VERBOSE_DEBUG(
                for (itemplate=0, itemplate_all=itemplates_amount*ntemplates_amount; itemplate<itemplates_amount_length; itemplate++, itemplate_all++) {
                    std::cout << "    template " << itemplate_all;
                    if (!s_templates.at(itemplate)) {
                        CORR3D_VERBOSE_DEBUG( std::cout << " Not loaded: " << templates_loaded_errormsg.at(itemplate) << ": SKIP" << std::endl; );
                        continue;
                    }
                    std::cout << ": " << intemplates.at(itemplate_all)->getName() << std::endl;
                }
            );

            // Cut out only the relevant wedges.
            if (wedge_particles) {
                wedge_particles_auto.reset(new std::vector<tom::Wedge<T> *>(iparticles_amount_length));
                for (std::size_t i=0; i<iparticles_amount_length; i++) {
                    wedge_particles_auto.get()->at(i) = wedge_particles->at(iparticles_amount*nparticles_amount+i);
                }
                pwedge_particles = wedge_particles_auto.get();
            }

            corr3d_single(s_templates, s_particle_fs, variance_particles, vangles, itemplates_amount*ntemplates_amount, iparticles_amount*nparticles_amount, NULL, mask,
                        pwedge_templates, pwedge_particles, corrhdl, fftw_flags, true);
        }
    }
}













template void tom::transf_particle<float , float >(const tom::Volume<float > &src, tom::Volume<float > &dst, const tom::Volume<float > *mask, double phi, double psi, double theta, double shiftx, double shifty, double shiftz);
template void tom::transf_particle<double, double>(const tom::Volume<double> &src, tom::Volume<double> &dst, const tom::Volume<double> *mask, double phi, double psi, double theta, double shiftx, double shifty, double shiftz);


template
void tom::align<float , float >(const std::vector<std::string> &fparticles,
                                const std::vector<std::string> &ftemplates,
                                const tom::Volume<double> &vangles,
                                const tom::Volume<float > *mask,
                                const tom::Volume<float > &mask_cc,
                                const std::vector<tom::Wedge<float > *> wedge_particles,
                                const std::vector<tom::Wedge<float > *> wedge_templates,
                                std::vector<std::vector<tom::cc_peak<float > > > &list_peak, unsigned fftw_flags);

template
void tom::align<double, double>(const std::vector<std::string> &fparticles,
                                const std::vector<std::string> &ftemplates,
                                const tom::Volume<double> &vangles,
                                const tom::Volume<double> *mask,
                                const tom::Volume<double> &mask_cc,
                                const std::vector<tom::Wedge<double> *> wedge_particles,
                                const std::vector<tom::Wedge<double> *> wedge_templates,
                                std::vector<std::vector<tom::cc_peak<double> > > &list_peak, unsigned fftw_flags);


template
void tom::corr3d<float , float >(   std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                                    const std::vector<tom::VolumeContainer<float > *> &fparticles,
                                    std::size_t nparticles_amount,
                                    const std::vector<tom::VolumeContainer<float > *> &ftemplates,
                                    std::size_t ntemplates_amount,
                                    const tom::Volume<double> &vangles,
                                    const tom::Volume<float > *mask,
                                    const std::vector<tom::Wedge<float > *> *wedge_particles,
                                    const std::vector<tom::Wedge<float > *> *wedge_templates,
                                    tom::CorrelationHandler<float > *corrhdl,
                                    unsigned fftw_flags);

template
void tom::corr3d<double, double>(   std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                                    const std::vector<tom::VolumeContainer<double> *> &fparticles,
                                    std::size_t nparticles_amount,
                                    const std::vector<tom::VolumeContainer<double> *> &ftemplates,
                                    std::size_t ntemplates_amount,
                                    const tom::Volume<double> &vangles,
                                    const tom::Volume<double> *mask,
                                    const std::vector<tom::Wedge<double> *> *wedge_particles,
                                    const std::vector<tom::Wedge<double> *> *wedge_templates,
                                    tom::CorrelationHandler<double> *corrhdl,
                                    unsigned fftw_flags);


template
void tom::average<float , float >(  const std::vector<std::string> &fparticles,
                                    std::vector<double3> angles,
                                    std::vector<double3> shifts,
                                    const tom::Volume<float > *mask,
                                    tom::Volume<float > &sum);
template
void tom::average<double, double>(  const std::vector<std::string> &fparticles,
                                    std::vector<double3> angles,
                                    std::vector<double3> shifts,
                                    const tom::Volume<double> *mask,
                                    tom::Volume<double> &sum);


template
bool tom::average_wedges<float >(   const std::vector<tom::Wedge<float > *> &wedge_particles,
                                    std::vector<double3> angles,
                                    tom::Volume<float > &sum,
                                    std::size_t sizez);

template
bool tom::average_wedges<double>(   const std::vector<tom::Wedge<double> *> &wedge_particles,
                                    std::vector<double3> angles,
                                    tom::Volume<double> &sum,
                                    std::size_t sizez);



template
void tom::corr3d_batch<float , float >(   std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                                    const std::vector<tom::VolumeContainer<float > *> &fparticles,
                                    std::size_t nparticles_amount,
                                    const std::vector<tom::VolumeContainer<float > *> &ftemplates,
                                    std::size_t ntemplates_amount,
                                    const tom::Volume<double> &vangles,
                                    const tom::Volume<float > *mask,
                                    const std::vector<tom::Wedge<float > *> *wedge_particles,
                                    const std::vector<tom::Wedge<float > *> *wedge_templates,
                                    tom::CorrelationHandler<float > *corrhdl,
                                    unsigned fftw_flags);
template
void tom::corr3d_batch<double, double>(   std::size_t sizex__, std::size_t sizey__, std::size_t sizez__,
                                    const std::vector<tom::VolumeContainer<double> *> &fparticles,
                                    std::size_t nparticles_amount,
                                    const std::vector<tom::VolumeContainer<double> *> &ftemplates,
                                    std::size_t ntemplates_amount,
                                    const tom::Volume<double> &vangles,
                                    const tom::Volume<double> *mask,
                                    const std::vector<tom::Wedge<double> *> *wedge_particles,
                                    const std::vector<tom::Wedge<double> *> *wedge_templates,
                                    tom::CorrelationHandler<double> *corrhdl,
                                    unsigned fftw_flags);



