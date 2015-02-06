/****************************************************************************//**
 * \file corrsin.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    18.01.2008
 *******************************************************************************/


#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <vector>
#include <iostream>
#include <fftw3.h>

#define __M_PI 3.141592653589793238512808959406186204433

#include "tom/core/fftw_plan.hpp"
#include "tom/core/volume_fcn.hpp"
#include "tom/core/wedge.hpp"

#include "tom/corr/corr.hpp"

#include "helper/auto_vector.hpp"


#define TIC { ____global_start = clock(); }
#define TOC { ____global_end = clock(); }
#define cpu_time_used (((double) (____global_end - ____global_start)) / CLOCKS_PER_SEC)
clock_t ____global_start, ____global_end;



const char *wisdom_file = "/fs/pool/pool-bmsan/haller/DA/tomc/.libs/wisdom_file.txt";
#if 1
    const char *filename_read = "/fs/pool/pool-nickell/thermoplasma/em/tomography/data/052305/vols/052305_big.vol";
#elif 1
    const char *filename_read = "/fs/pool/pool-bmsan/haller/DA/data/vols/v1.em";
#elif 0
    const char *filename_read = "/fs/pool/pool-nickell/thermoplasma/em/tomography/template_matching/sub_volumes/job/template/Ribo_0_68_nm_fil_0_6.em";
#else
    const char *filename_read = "/fs/pool/pool-bmsan/haller/DA/data/vols/testvol.em";
    const char *filename_read = "/fs/pool/pool-bmsan/haller/DA/data/vols/testvol_double.em";
#endif
const char *filename_write_re = "/fs/pool/pool-bmsan/haller/DA/data/vols/data_re.em";
const char *filename_write_im = "/fs/pool/pool-bmsan/haller/DA/data/vols/data_im.em";


#define Str(x) #x
#define Xstr(x) Str(x)


#define FLOATTYPE double


#include <sys/time.h>


template<typename T>
T own_rand() {
    return drand48();
}


int main(int argc, char *argv[]) {

    tom::Volume<double> v1(10,10,10, NULL,NULL);
    tom::Volume<double> v2(10,10,10, NULL,NULL);
    tom::element_wise_add(v1, v2);
    return 0;



    std::cout.precision(20);

    unsigned const fftw_flags = FFTW_MEASURE;

    typedef FLOATTYPE T;

    printf("%% set wisdom %s: (%d)\n", wisdom_file, tom::fftw::set_wisdom_name(wisdom_file, 1));
    printf("%% work in %s modus\n", Xstr(FLOATTYPE));


    std::vector<std::string> fparticles, ftemplates;
    std::auto_ptr<tom::Volume<double> > angles;

    for (int i=1; i<=153; i++) {
        std::stringstream s;
        s << "/fs/home/haller/haller/DA/data/particles/particle_" << i << ".em";
        //s << "/fs/home/haller/haller/DA/data/gen_part/avg_" << i << ".em";
        fparticles.push_back(s.str());
    }

    ftemplates.push_back("/fs/home/haller/haller/DA/data/avg.em");
    ftemplates.push_back("/fs/home/haller/haller/DA/data/avg_mickey.em");


    {
        tom::Volume<double> *pangles;
        tom_io_em_header header;
        const char *f = NULL;
        f = "/fs/home/haller/haller/DA/data/gen_part/avg_angles.em";
        f = "/fs/home/haller/haller/DA/data/angles.em";
        f = "/fs/home/haller/haller/DA/data/angles_rand2000.em";
        try {
            tom::read_from_em<double>(pangles, f, NULL, NULL, NULL, &header, NULL, NULL);
            angles.reset(pangles);
        } catch (int &res) {
            std::cout << "exception reading the angles-emfile: " << f << " (" << res << ")" << std::endl;
            throw;
        } catch (std::exception &e) {
            std::cout << "exception reading the angles-emfile: " << f << " (" << e.what() << ")" << std::endl;
            throw;
        }
    }

    std::vector<std::vector<tom::cc_peak<T> > > list_peak;


    std::size_t dims[3] = {32,32,32};

    tom::Volume<T> mask(dims[0],dims[1],dims[2],NULL,NULL);
    #define min(a,b) ((a) < (b) ? (a) : (b))
    #define PI 3.141592653589793238512808959406186204433
    //tom::init_spheremask(mask, ceil( min(dims[0], min(dims[1],dims[2]))/2.)-1. -2, 2., ceil( min(dims[0], min(dims[1],dims[2]))/2.));
    tom::init_spheremask(mask, 15, 0, 0);

    tom::Volume<T> mask_cc(dims[0],dims[1],dims[2],NULL,NULL);
    mask_cc.setValues(mask);
    mask_cc.setValues(1);

    std::vector<tom::Wedge<T> *> wedge_particles(fparticles.size());
    std::vector<tom::Wedge<T> *> wedge_templates(ftemplates.size());
    std::auto_ptr<tom::Wedge<T> > wedge_particle(new tom::SimpleWedge<T>(30./180.*PI, (min(dims[0], min(dims[1],dims[2]))-1) / 2 ));
    std::auto_ptr<tom::Wedge<T> > wedge_template(new tom::NoWedge<T>());

    //wedge_particle.reset(new tom::NoWedge<T>());
    //wedge_template.reset(new tom::SimpleWedge<T>(30./180.*PI, (min(dims[0], min(dims[1],dims[2]))-1) / 2 ));

    for (size_t i=0; i<fparticles.size(); i++) {
        wedge_particles.at(i) = wedge_particle.get();
    }
    for (size_t i=0; i<ftemplates.size(); i++) {
        wedge_templates.at(i) = wedge_template.get();
    }

    #if 0
    tom::align<T, T>(fparticles, ftemplates, *angles, &mask, mask_cc, wedge_particles, wedge_templates, list_peak, fftw_flags);
    #else
    {
        tom::CorrelationHandlerPeaks<T> corrhdl(ftemplates.size(), fparticles.size(), mask_cc, true);
        std::vector<tom::VolumeContainer<T> *> intemplates(ftemplates.size());
        std::vector<tom::VolumeContainer<T> *> inparticles(fparticles.size());
        auto_vector<tom::VolumeContainer<T> > auto_inparticles(fparticles.size());
        auto_vector<tom::VolumeContainer<T> > auto_intemplates(ftemplates.size());
        std::auto_ptr<tom::VolumeContainer<T> > ptr;
        for (std::size_t i=0; i<ftemplates.size(); i++) {
            ptr.reset(new tom::VolumeContainerEM<T>(ftemplates.at(i)));
            intemplates.at(i) = ptr.get();
            auto_intemplates.assign(i, ptr);
        }
        for (std::size_t i=0; i<fparticles.size(); i++) {
            ptr.reset(new tom::VolumeContainerEM<T>(fparticles.at(i)));
            inparticles.at(i) = ptr.get();
            auto_inparticles.assign(i, ptr);
        }
        //tom::corr3d<T, T>(dims[0], dims[1], dims[2], inparticles, 20, intemplates, 10, *angles, &mask, &wedge_particles, &wedge_templates, &corrhdl, fftw_flags);
        tom::corr3d_batch<T, T>(dims[0], dims[1], dims[2], inparticles, 20, intemplates, 10, *angles, &mask, &wedge_particles, &wedge_templates, &corrhdl, fftw_flags);

        list_peak = corrhdl.getPeakList();
    }
    #endif

    std::cout << "------------------------------------------" << std::endl << std::endl;
    #define PI 3.141592653589793238512808959406186204433
    std::cout << "FINISHED: " << std::endl;
    for (int i=0; (size_t)i<ftemplates.size(); i++) {
        std::cout << "  template " << i << " (\"" << ftemplates[i] << "\"): " << std::endl;
        for (int j=0; (size_t)j<fparticles.size(); j++) {
            std::cout << "    particle " << j << " (\"" << fparticles[j] << "\"): ";
            tom::cc_peak<T> &peak = list_peak.at(i).at(j);
            if (peak.angle_idx>=0) {
                std::cout << "[" << peak.x << "," << peak.y << "," << peak.z << "]: " << peak.val << " (" << peak.angle_idx << " = [" << (angles->get(0,peak.angle_idx,0)*180./PI) << "," << (angles->get(1,peak.angle_idx,0)*180./PI) << "," << (angles->get(2,peak.angle_idx,0)*180./PI) << "] deg)" << std::endl;
            } else {
                std::cout << "ERROR" << std::endl;
            }
        }
    }


    {
        tom::Volume<T> avg(dims[0], dims[1], dims[2], NULL, NULL);
        tom::Volume<T> avg_wedge(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
        std::cout << "------------------------------------------" << std::endl << std::endl;
        for (std::size_t itemplate=0; itemplate<ftemplates.size(); itemplate++) {
            std::vector<std::string> fparticles_clean;
            std::vector<tom::double3> angles_clean;
            std::vector<tom::double3> shifts_clean;
            std::vector<tom::Wedge<T> *> wedge_particles_clean;
            for (std::size_t iparticle=0; iparticle<fparticles.size(); iparticle++) {
                if (list_peak.at(itemplate).at(iparticle).angle_idx >= 0) {
                    fparticles_clean.push_back(fparticles.at(iparticle));
                    tom::double3 t;
                    t.first  = - angles->get(1,list_peak.at(itemplate).at(iparticle).angle_idx,0);
                    t.second = - angles->get(0,list_peak.at(itemplate).at(iparticle).angle_idx,0);
                    t.third  = - angles->get(2,list_peak.at(itemplate).at(iparticle).angle_idx,0);
                    angles_clean.push_back(t);
                    t.x = - list_peak.at(itemplate).at(iparticle).x;
                    t.y = - list_peak.at(itemplate).at(iparticle).y;
                    t.z = - list_peak.at(itemplate).at(iparticle).z;
                    shifts_clean.push_back(t);
                    wedge_particles_clean.push_back(wedge_particles.at(iparticle));
                }
            }
            std::cout << "average for template " << itemplate << " (\"" << ftemplates.at(itemplate) << "\")";
            if (ftemplates.empty()) {
                std::cout << ": NO particles aligned. skip." << std::endl;
            } else {
                std::cout << " with " << fparticles_clean.size() << " out of " << fparticles.size() << " particles" << std::endl;
                tom::Volume<std::complex<T> > avg_fs(dims[0], dims[1], dims[2]/2+1, NULL, NULL);
                tom::fftw::Plan<T> plan_r2c(avg, avg_fs, fftw_flags);
                tom::fftw::Plan<T> plan_c2r(avg_fs, avg, fftw_flags);
                tom::average<T,T>(fparticles_clean, angles_clean, shifts_clean, &mask, avg);

                avg.write_to_em("/fs/home/haller/haller/DA/data/vols/run_avg.em", NULL);

                if (tom::average_wedges<T>(wedge_particles_clean, angles_clean, avg_wedge, dims[2])) {
                    plan_r2c.execute(avg, avg_fs);

                    avg_wedge.write_to_em("/fs/home/haller/haller/DA/data/vols/run_avg_wedge.em", NULL);
                    tom::element_wise_div<std::complex<T>, T>(avg_fs, avg_wedge, std::complex<T>(0,0));

                    plan_c2r.execute(avg_fs, avg);
                    avg.shift_scale<double>(0., 1./(dims[0]*dims[1]*dims[2]));
                    avg.write_to_em("/fs/home/haller/haller/DA/data/vols/run_avg_with_wedge.em", NULL);
                }
            }
        }
    }



    printf("save wisdom %s: (%d)\n", wisdom_file, tom::fftw::save_wisdom());

    tom::fftw::clear_wisdom_name();
    fftw_forget_wisdom();

    return 0;

}
