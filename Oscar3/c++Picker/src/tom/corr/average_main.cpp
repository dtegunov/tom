/****************************************************************************//**
 * \file average_main.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    15.02.2008
 *******************************************************************************/





#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <set>

#include "boost/tuple/tuple.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/shared_ptr.hpp"

#include "helper/snippets.hpp"
#include "helper/ostream_swapper.hpp"


#include "tom/core/volume_fcn.hpp"
#include "tom/core/volume_loop.hpp"
#include "tom/corr/config_files.hpp"
#include "tom/corr/filename_generators.hpp"
#include "tom/corr/average.hpp"



namespace {
template<typename T, typename TMASK>
struct st_binary_mask {
    T v_;
    st_binary_mask(T v_): v_(v_) { }
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (!mask) {
            v = v_;
        }
    }
};
/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename T, typename TMASK>
void set_where_binary_mask_not_set(tom::Volume<T> &v, const tom::Volume<TMASK> &mask, T zero_value) {
    tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::st_binary_mask<T, TMASK> >(v, mask, ::st_binary_mask<T, TMASK>(zero_value));
}
} // namespace





/****************************************************************************//**
 *
 *
 *******************************************************************************/
int main(int argc, char *argv[]) {

    typedef double TFLOAT;

    if (argc<5 || argc>8) {
        std::cerr << "Usage: " << (argc?argv[0]:"average") << " CONFIG_FILE FTEMPLATES FPARTICLES [LOGFILE] [DO_FSC] [WEDGE_THRESHOLD]"               "\n"
                     "compiled at " __DATE__ " " __TIME__ " from " __FILE__                                                         "\n"
                     "  CONFIG_FILE     : filename of the parallel correlation program."                                                 "\n"
                     "  FTEMPLATES      : filename of the textfile with the names of the templates."                                     "\n"
                     "  FPARTICLES      : filename of the textfile with the names of the particles."                                     "\n"
                     "  ANGLES          : filename of the textfile with the angles."                                                     "\n"
                     "  LOGFILE         : optional parameter for the log file. '-' for STDOUT."                                          "\n"
                     "  DO_FSC          : optional parameter to decide whether to compute the FSC (yes or no)."                          "\n"
                     "  WEDGE_THRESHOLD : Voxels of the averaged wedge are set to zero if they are SMALLER then WEDGE_THRESHOLD% of the number of particles (default 0 == no cutoff)."  << std::endl;
        return -1;
    }

    double use_threshold = 0.005;
    if (argc >= 8) {
        try {
            use_threshold = boost::lexical_cast<double>(argv[7]);
        } catch (boost::bad_lexical_cast &) {
            throw std::invalid_argument("The WEDGE_THRESHOLD must be a floating point number.");
        }
        if (use_threshold > 1 && use_threshold <= 100) {
            use_threshold = use_threshold/100.;
        }
        if (use_threshold > 1 || use_threshold <= 0) {
            throw std::invalid_argument("WEDGE_THRESHOLD must be either a prozentual value ]1, 100] or relative between ]0.0, 1.0].");
        }
    }

    bool do_fourier_shell = true;
    if (argc >= 7) {
        std::string p(argv[6]);
        std::transform(p.begin(), p.end(), p.begin(), &toupper);
        if (p=="YES" || p=="1" || p=="FSC" || p=="TRUE") {
            do_fourier_shell = true;
        } else if (p=="NO" || p=="0" || p=="FALSE") {
            do_fourier_shell = false;
        } else {
            std::cerr << "Argument 6 (do_fsc) must be either \"yes\" or \"no\"." << std::endl;
            return -2;
        }
    }


    std::size_t itemplate, iparticle, ntemplates, nparticles, i;

    std::vector<std::string> ftemplates, fparticles;
    std::vector<boost::shared_ptr<tom::Volume<double> > > angles;
    std::vector<helper::triple<double, double, double> > anglesv;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > wedge_templates, wedge_particles;
    std::vector<tom::cc_peak<TFLOAT> > peak_list;
    tom::corr::ConfigDataServer c;
    std::string outputdir;

    std::stringstream clog;
    std::auto_ptr<std::ofstream> f;
    std::auto_ptr<helper::ostream_swapper> stream_swapper;
    clog.precision(20);
    clog << "# start " << (argc>0?argv[0]:"average") << " at " << helper::now2str() << " (compiled at " __DATE__ " " __TIME__ ")." << std::endl;

    // Parse the input files.
    clog << "# parse configuration (" << argv[1] << ")" << std::endl;
    c.read_from_conf(argv[1]);
    clog << "# parse filelist for templates (" << argv[2] << ")" << std::endl;
    tom::corr::parse_filelist(argv[2], ftemplates, wedge_templates, 1);
    ntemplates = ftemplates.size();
    clog << "#   " << ntemplates << " templates read.\n"
            "# parse filelist for particles (" << argv[3] << ")" << std::endl;
    tom::corr::parse_filelist(argv[3], fparticles, wedge_particles, 1);
    nparticles = fparticles.size();
    clog << "#   " << nparticles << " particles read.\n"
            "# parse the list of angles (" << argv[4] << ")" << std::endl;
    tom::corr::parse_angleslist(argv[4], ntemplates, nparticles, angles);
    clog << "# parse peak list (" << c.peakfilename << ")" << std::endl;
    tom::corr::parse_peaklist(c.peakfilename, peak_list, anglesv, ntemplates, nparticles);

    clog << "#\n# use a threshold for the wedge of " << use_threshold << "\n# (values in the summed wedge are set to zero if they are smaller then the relative threshold compared to the number of particles).\n" << std::endl;

    tom::avg::filename_generator f_gen(c.outputdir);


    // open the logfile.
    {
        std::ostream *use_log = &std::clog;

        std::string logfilename(argc>=6 ? argv[5] : "-");
        if (logfilename != "-") {
            f.reset(new std::ofstream(logfilename.c_str(), std::ios::out | std::ios::trunc));
            if (!f->good()) {
                f.reset(NULL);
                clog << "# ERROR: could not open log-file \"" << logfilename << "\" for writing. Use CERR instead." << std::endl;
            } else {
                use_log = f.get();
            }
        }
        *use_log << clog.rdbuf() << std::flush;
        clog.clear();
        stream_swapper.reset(new helper::ostream_swapper(clog, *use_log));
    }



    // Find out, which particles average to which template.
    std::vector<std::vector<std::size_t> > use_particles(ntemplates), use_particles1(ntemplates), use_particles2(ntemplates), use_particles3(ntemplates);
    for (iparticle=0; iparticle<nparticles; iparticle++) {
        std::size_t max_index = ntemplates;
        double max_ccvalue = - std::numeric_limits<double>::infinity();
        for (itemplate=0; itemplate<ntemplates; itemplate++) {
            i = itemplate * nparticles + iparticle;

            if (angles[i].get()) {
                assert(angles[i]->getSizeX()==3 && angles[i]->getSizeZ()==1);

                if (peak_list[i].angle_idx < 0) {
                    clog << "# WARNING: no peak found for template #" << itemplate << ", particle #" << iparticle << "." << std::endl;
                    continue;
                } else {
                    if (static_cast<std::size_t>(peak_list[i].angle_idx) >= angles[i]->getSizeY()) {
                        clog << "# WARNING: peak-index for template #" << itemplate << ", particle #" << iparticle << " is larger then the entire angle-list (" << peak_list[i].angle_idx << " vs. " << angles[i]->getSizeY() << "). Using peak from file." << std::endl;
                    } else {
                        boost::tuples::tuple<double &, double &, double &> a = boost::tuples::tie(angles[i]->get(0, peak_list[i].angle_idx, 0), angles[i]->get(0, peak_list[i].angle_idx, 0), angles[i]->get(0, peak_list[i].angle_idx, 0));
                        if (fabs(a.get<0>()-anglesv[i].x)>1e-7 || fabs(a.get<0>()-anglesv[i].x)>1e-7 || fabs(a.get<0>()-anglesv[i].x)>1e-7) {
                            clog << "# WARNING: peak-index for template #" << itemplate << ", particle #" << iparticle << " differs from angles list and peakfile: (" << a.get<0>() << " vs. " << anglesv[i].x << "; " << a.get<0>() << " vs. " << anglesv[i].x << "; " << a.get<0>() << " vs. " << anglesv[i].x << "). Take from peak-file." << std::endl;
                        }
                    }
                }
            } else {
                if (peak_list[i].angle_idx >= 0) {
                    clog << "# WARNING: peak-index for template #" << itemplate << ", particle #" << iparticle << " given, but the angles list is empty. Using peak from file." << std::endl;
                } else {
                    // everything is fine.
                    continue;
                }
            }
            const double shiftx = - peak_list[i].x;
            const double shifty = - peak_list[i].y;
            const double shiftz = - peak_list[i].z;
            if (fabs(shiftx) > c.volsize/2 || fabs(shifty) > c.volsize/2 || fabs(shiftz) > c.volsize/2) {
                clog << "# WARNING: the shift of template #" << itemplate << ", particle #" << iparticle << " is too large for the volume of size " << c.volsize << " (" << shiftx << "," << shifty << "," << shiftz << ")." << std::endl;
                continue;
            }

            if (max_index==ntemplates || max_ccvalue<peak_list[i].val) {
                max_ccvalue = peak_list[i].val;
                max_index = itemplate;
            }
        }
        if (max_index != ntemplates) {
            use_particles[max_index].push_back(iparticle);
        }
    }

    // Split the list of particles in two sets for the FSC.
    if (do_fourier_shell) {
        for (itemplate=0; itemplate<ntemplates; itemplate++) {
            std::size_t size = use_particles[itemplate].size();
            if (size >= 2) {
                const std::vector<std::size_t> &v_ = use_particles [itemplate];
                std::vector<std::size_t> &v1 = use_particles1[itemplate];
                std::vector<std::size_t> &v2 = use_particles2[itemplate];
                std::vector<std::size_t> &v3 = use_particles3[itemplate];
                v1.resize(0);
                v2.resize(0);
                v3.resize(0);
                if (size%2) {
                    size--;
                    v3.push_back(v_[size]);
                }
                v1.reserve(size/2);
                v2.reserve(size/2);
                for (iparticle=0; iparticle<size; iparticle+=2) { v1.push_back(v_[iparticle]); }
                for (iparticle=1; iparticle<size; iparticle+=2) { v2.push_back(v_[iparticle]); }
            }
        }
    } else {
        clog << "# don't compute the Fourier Shell Correlation." << std::endl;
    }


    // Init the spheremask
    std::auto_ptr<const tom::Volume<TFLOAT> > mask_sphere;
    if (c.sphere_mask_inner_radius>0 || c.sphere_mask_sigma>0) {
        tom::Volume<TFLOAT> *p;
        mask_sphere.reset(p = new tom::Volume<TFLOAT>(c.volsize, c.volsize, c.volsize, NULL,NULL));
        tom::init_spheremask(*p, c.sphere_mask_inner_radius, c.sphere_mask_sigma, c.sphere_mask_cutoff_radius);
    }



    tom::Volume<TFLOAT> avg(c.volsize, c.volsize, c.volsize, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>());
    tom::Volume<TFLOAT> avgwedge(c.volsize, c.volsize, c.volsize, NULL,NULL);
    tom::Volume<TFLOAT> avgwedge_half(avgwedge, NULL, c.volsize, c.volsize, c.volsize/2+1, avgwedge.getStrideX(), avgwedge.getStrideY(), avgwedge.getStrideZ());
    tom::Volume<TFLOAT> avgwedge_sqrtn(c.volsize, c.volsize, c.volsize, NULL,NULL);
    tom::Volume<TFLOAT> avgwedge_sqrtn_half(avgwedge_sqrtn, NULL, c.volsize, c.volsize, c.volsize/2+1, avgwedge.getStrideX(), avgwedge.getStrideY(), avgwedge.getStrideZ());
    tom::Volume<TFLOAT> avg2(c.volsize, c.volsize, c.volsize, NULL,NULL);
    tom::Volume<TFLOAT> avgwedge2(c.volsize, c.volsize, c.volsize, NULL,NULL);
    tom::Volume<TFLOAT> avgwedge2_half(avgwedge2, NULL, c.volsize, c.volsize, c.volsize/2+1, avgwedge2.getStrideX(), avgwedge2.getStrideY(), avgwedge2.getStrideZ());
    tom::Volume<std::complex<TFLOAT> > avg_ft(c.volsize, c.volsize, c.volsize/2+1, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>());
    tom::Volume<std::complex<TFLOAT> > avg_ft2(c.volsize, c.volsize, c.volsize/2+1, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>());

    tom::fftw::Plan<TFLOAT> plan_fwd(avg2, avg_ft, c.fftw_flag | FFTW_DESTROY_INPUT);
    tom::fftw::Plan<TFLOAT> plan_bwd(avg_ft, avg2, c.fftw_flag | FFTW_DESTROY_INPUT);

    std::auto_ptr<tom::Volume<TFLOAT> > vtemplate;
    tom::Volume<std::complex<TFLOAT> > vtemplate_ft(c.volsize, c.volsize, c.volsize/2+1, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>());
    tom::fftw::Plan<TFLOAT> plan_vtemplate(avg, vtemplate_ft, c.fftw_flag | FFTW_DESTROY_INPUT);


    // Volumes for the fourier shell correlation.
    std::auto_ptr<tom::Volume<TFLOAT> > FSC_v, FSC_w2;
    std::auto_ptr<tom::Volume<std::complex<TFLOAT> > > FSC_f1, FSC_f2, FSC_f1_half, FSC_f2_half;
    std::auto_ptr<tom::fftw::Plan<TFLOAT> > FSC_plan1, FSC_plan2;
    if (do_fourier_shell) {
        FSC_v.reset(        new tom::Volume<TFLOAT>(c.volsize, c.volsize, c.volsize, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>()));
        FSC_f1.reset(       new tom::Volume<std::complex<TFLOAT> >(c.volsize, c.volsize, c.volsize, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>()));
        FSC_f2.reset(       new tom::Volume<std::complex<TFLOAT> >(c.volsize, c.volsize, c.volsize, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>()));
        FSC_f1_half.reset(  new tom::Volume<std::complex<TFLOAT> >(*FSC_f1, NULL, c.volsize, c.volsize, c.volsize/2+1, FSC_f1->getStrideX(), FSC_f1->getStrideY(), FSC_f1->getStrideZ()));
        FSC_f2_half.reset(  new tom::Volume<std::complex<TFLOAT> >(*FSC_f2, NULL, c.volsize, c.volsize, c.volsize/2+1, FSC_f2->getStrideX(), FSC_f2->getStrideY(), FSC_f2->getStrideZ()));
        FSC_plan1.reset(    new tom::fftw::Plan<TFLOAT>(*FSC_v, *FSC_f1_half, c.fftw_flag | FFTW_DESTROY_INPUT));
        FSC_plan2.reset(    new tom::fftw::Plan<TFLOAT>(*FSC_v, *FSC_f2_half, c.fftw_flag | FFTW_DESTROY_INPUT));
        FSC_w2.reset(       new tom::Volume<TFLOAT>(c.volsize, c.volsize, c.volsize/2+1, tom::fftw::fftw_malloc<TFLOAT>(), tom::fftw::fftw_free<TFLOAT>()));
    }


    double mean, variance;

    std::vector<std::size_t>::const_iterator culvit; // Const Unsigned Long Vector ITerator

    std::vector<tom::avg::st_average_result> av(ntemplates);

    double threshold = 1;

    // Now compute the average, the FSC, the CC and save it to file.
    for (itemplate=0; itemplate<ntemplates; itemplate++) {
        const std::vector<std::size_t> &v_ = use_particles [itemplate];
        const std::vector<std::size_t> &v1 = use_particles1[itemplate];
        const std::vector<std::size_t> &v2 = use_particles2[itemplate];
        const std::vector<std::size_t> &v3 = use_particles3[itemplate];

        threshold = v_.size()*use_threshold; // cut off.

        assert(std::set<std::size_t>(v_.begin(), v_.end()).size() == v_.size() &&
               (!do_fourier_shell || (v_.size()<2 || v_.size()==(v1.size()+v2.size()+v3.size()))));

        if (v_.empty()) {
            continue;
        }

        tom::avg::st_average_result &ar = av[itemplate];
        ar.use_idx = v_;

        // Read the original template from file.
        try {
            tom::Volume<TFLOAT> *p;
            tom::read_from_em(p, ftemplates[itemplate], NULL,NULL,NULL, NULL, NULL,NULL);
            vtemplate.reset(p);
        } catch (std::exception &e) {
            std::stringstream ss;
            ss << "# The template #" << itemplate << " (" << ftemplates[itemplate] << ") does not exist.";
            if (c.force_files_exist) {
                clog << ss.str() << " ABORT";
                throw;
            }
            clog << ss << " SKIP CC correlation.";
            vtemplate.reset(NULL);
        }
        if (vtemplate.get()) {
            // normalise it and compute its fourier transformation.
            vtemplate->stat(mean, variance, false);
            vtemplate->shift_scale<TFLOAT>(-mean, 1./sqrt(variance));
            avg.setValues(*vtemplate);
            plan_vtemplate.execute(avg, vtemplate_ft);
        }


        // Average the particles.

        if (!v1.empty() && v1.size()==v2.size()) {
            assert(do_fourier_shell);
            tom::avg::average(c.volsize, fparticles, wedge_particles, v1, itemplate, mask_sphere.get(), peak_list, anglesv, c.force_files_exist, &clog, avg, avgwedge_half);
            tom::avg::average(c.volsize, fparticles, wedge_particles, v2, itemplate, mask_sphere.get(), peak_list, anglesv, c.force_files_exist, &clog, avg2, avgwedge2_half);

            std::vector<double> res;
            std::size_t n_shells = c.volsize/2;

            FSC_w2->setValues(avgwedge_half);
            tom::element_wise_set_below_threshold<TFLOAT>(*FSC_w2, threshold*v1.size()/static_cast<double>(v_.size()), 0);
            tom::avg::apply_avgwedge(avg, *FSC_w2, mask_sphere.get(), c.fftw_flag, *FSC_v);
            //FSC_v->setValues(avg);
            //FSC_v->write_to_em("FSC_v1.em", NULL);
            FSC_plan1->execute(*FSC_v, *FSC_f1_half);

            FSC_w2->setValues(avgwedge2_half);
            tom::element_wise_set_below_threshold<TFLOAT>(*FSC_w2, threshold*v2.size()/static_cast<double>(v_.size()), 0);
            tom::avg::apply_avgwedge(avg2, *FSC_w2, mask_sphere.get(), c.fftw_flag, *FSC_v);
            //FSC_v->setValues(avg2);
            //FSC_v->write_to_em("FSC_v2.em", NULL);
            FSC_plan2->execute(*FSC_v, *FSC_f2_half);

            tom::fourier_shell_correlation(*FSC_f1_half, false, *FSC_f2_half, false, c.volsize, n_shells, res);

            ar.fsc.resize(n_shells);
            for (std::size_t i=0; i<n_shells; i++) {
                ar.fsc[i].first  = res[i+n_shells*1];
                ar.fsc[i].second = res[i+n_shells*7];
            }

            //tom::Volume<double>(&res[0], n_shells, 11, 1, 0,0,0, false, NULL).write_to_em("FSC_cc.em", NULL);

            tom::element_wise_add<TFLOAT, TFLOAT>(avg, avg2);
            tom::element_wise_add<TFLOAT, TFLOAT>(avgwedge_half, avgwedge2_half);
            if (!v3.empty()) {
                // Add the missing partivle to the average.
                tom::avg::average(c.volsize, fparticles, wedge_particles, v3, itemplate, mask_sphere.get(), peak_list, anglesv, c.force_files_exist, &clog, avg2, avgwedge2_half);
                tom::element_wise_add<TFLOAT, TFLOAT>(avg, avg2);
                tom::element_wise_add<TFLOAT, TFLOAT>(avgwedge_half, avgwedge2_half);
            }
        } else {
            tom::avg::average(c.volsize, fparticles, wedge_particles, v_, itemplate, mask_sphere.get(), peak_list, anglesv, c.force_files_exist, &clog, avg, avgwedge_half);
        }

        ////////////////////////////////////////////////////////////////////////
        // Save the average volume (the pure sum): avg_s
        avg.write_to_em(f_gen.get_avg_s(itemplate), NULL);
        tom::hermitian_symmetry_to_full(avgwedge);
        tom::fftshift(avgwedge, avgwedge2, true);
        avgwedge2.write_to_em(f_gen.get_avgwedge_s(itemplate), NULL);
        if (vtemplate.get()) {
            // Compute the correlation with the original template.
            avg.stat(mean, variance, false);
            if (variance) {
                avg2.setValues(avg); // make a copy for normalising.
                avg2.shift_scale<TFLOAT>(-mean, 1./sqrt(variance));
                plan_fwd.execute(avg2, avg_ft);
                tom::element_wise_conj_multiply<std::complex<TFLOAT>, std::complex<TFLOAT> >(avg_ft, vtemplate_ft);
                plan_bwd.execute(avg_ft, avg2);
                ar.ccval_s.reset(new double((avg2.get(0,0,0)/avg2.numel())/avg2.numel()));
            }
        }
        ////////////////////////////////////////////////////////////////////////

        // Set values below threshold in volume to 0
        tom::element_wise_set_below_threshold<TFLOAT>(avgwedge_half, threshold, 0);
        tom::hermitian_symmetry_to_full(avgwedge);

        avg2.setValues(avg);
        plan_fwd.execute(avg2, avg_ft);
        avg_ft2.setValues(avg_ft);

        ////////////////////////////////////////////////////////////////////////
        // Compute the average with frequencies weighted by N
        set_where_binary_mask_not_set<std::complex<TFLOAT>, TFLOAT>(avg_ft, avgwedge_half, 0);
        avg_ft.get(0,0,0) = 0;
        plan_bwd.execute(avg_ft, avg2);
        if (mask_sphere.get()) {
            tom::norm_mask<TFLOAT, TFLOAT, double>(avg2, *mask_sphere, tom::norm::NORM_STDDEV_1, NULL, true);
        } else {
            if ((variance=avg2.variance_mean_free(false))) {
                avg2.shift_scale<TFLOAT>(0, 1./sqrt(variance));
            }
        }
        avg2.write_to_em(f_gen.get_avg_n(itemplate), NULL);
        tom::fftshift(avgwedge, avgwedge2, true);
        avgwedge2.write_to_em(f_gen.get_avgwedge_n(itemplate), NULL);
        if (vtemplate.get()) {
            // Compute the correlation with the original template.
            avg2.stat(mean, variance, false);
            if (variance) {
                avg2.shift_scale<TFLOAT>(-mean, 1./sqrt(variance));
                plan_fwd.execute(avg2, avg_ft);
                tom::element_wise_conj_multiply<std::complex<TFLOAT>, std::complex<TFLOAT> >(avg_ft, vtemplate_ft);
                plan_bwd.execute(avg_ft, avg2);
                ar.ccval_n.reset(new double((avg2.get(0,0,0)/avg2.numel())/avg2.numel()));
            }
        }
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Compute the average with frequencies weighted by sqrt(N)
        avgwedge_sqrtn_half.setValues(avgwedge_half);
        tom::element_wise_operation(avgwedge_sqrtn_half, &tom::math::sqrt<TFLOAT>);
        avg_ft.setValues(avg_ft2);
        tom::element_wise_div<std::complex<TFLOAT>, TFLOAT>(avg_ft, avgwedge_sqrtn_half, std::complex<TFLOAT>(0,0));
        avg_ft.get(0,0,0) = 0;
        plan_bwd.execute(avg_ft, avg2);
        if (mask_sphere.get()) {
            tom::norm_mask<TFLOAT, TFLOAT, double>(avg2, *mask_sphere, tom::norm::NORM_STDDEV_1, NULL, true);
        } else {
            if ((variance=avg2.variance_mean_free(false))) {
                avg2.shift_scale<TFLOAT>(0, 1./sqrt(variance));
            }
        }
        avg2.write_to_em(f_gen.get_avg_sqrtn(itemplate), NULL);
        tom::hermitian_symmetry_to_full(avgwedge_sqrtn);
        tom::fftshift(avgwedge_sqrtn, avgwedge2, true);
        avgwedge2.write_to_em(f_gen.get_avgwedge_sqrtn(itemplate), NULL);
        if (vtemplate.get()) {
            // Compute the correlation with the original template.
            avg2.stat(mean, variance, false);
            if (variance) {
                avg2.shift_scale<TFLOAT>(-mean, 1./sqrt(variance));
                plan_fwd.execute(avg2, avg_ft);
                tom::element_wise_conj_multiply<std::complex<TFLOAT>, std::complex<TFLOAT> >(avg_ft, vtemplate_ft);
                plan_bwd.execute(avg_ft, avg2);
                ar.ccval_sqrtn.reset(new double((avg2.get(0,0,0)/avg2.numel())/avg2.numel()));
            }
        }
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Compute the average with frequencies weighted by 1
        avgwedge2_half.setValues(avgwedge_half);
        avg_ft.setValues(avg_ft2);
        tom::element_wise_div<std::complex<TFLOAT>, TFLOAT>(avg_ft, avgwedge2_half, std::complex<TFLOAT>(0,0));
        avg_ft.get(0,0,0) = 0;
        plan_bwd.execute(avg_ft, avg2);
        if (mask_sphere.get()) {
            tom::norm_mask<TFLOAT, TFLOAT, double>(avg2, *mask_sphere, tom::norm::NORM_STDDEV_1, NULL, true);
        } else {
            if ((variance=avg2.variance_mean_free(false))) {
                avg2.shift_scale<TFLOAT>(0, 1./sqrt(variance));
            }
        }
        avg2.write_to_em(f_gen.get_avg_1(itemplate), NULL);
        tom::make_binary(avgwedge2);
        tom::hermitian_symmetry_to_full(avgwedge2);
        tom::fftshift(avgwedge2, avgwedge, true);
        avgwedge.write_to_em(f_gen.get_avgwedge_1(itemplate), NULL);
        if (vtemplate.get()) {
            // Compute the correlation with the original template.
            avg2.stat(mean, variance, false);
            if (variance) {
                avg2.shift_scale<TFLOAT>(-mean, 1./sqrt(variance));
                plan_fwd.execute(avg2, avg_ft);
                tom::element_wise_conj_multiply<std::complex<TFLOAT>, std::complex<TFLOAT> >(avg_ft, vtemplate_ft);
                plan_bwd.execute(avg_ft, avg2);
                ar.ccval_1.reset(new double((avg2.get(0,0,0)/avg2.numel())/avg2.numel()));
            }
        }
        ////////////////////////////////////////////////////////////////////////

    }

    tom::avg::write_average_to_log(clog, ftemplates, fparticles, av, c.outputdir);

    clog << "# finished at " << helper::now2str() << std::endl;

    return 0;
}








