/****************************************************************************//**
 * \file jobmanager_client.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    13.01.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__JOBMANAGER_CLIENT_HPP__
#define ___INCLUDE_CORR__JOBMANAGER_CLIENT_HPP__


#include <boost/shared_ptr.hpp>


#include "tom/core/wedge.hpp"
#include "tom/corr/correlation_handler.hpp"





namespace tom {
namespace corr {






class ConfigDataClient;
class ConfigDataServer;



/****************************************************************************//**
 *
 *******************************************************************************/
class ConfigDataClient {

public:
    unsigned long volsize;
    unsigned binning;
    double sphere_mask_inner_radius;
    double sphere_mask_sigma;
    double sphere_mask_cutoff_radius;
    double cc_mask_radius;
    unsigned fftw_flag;
    bool reapply_mask_after_wedge;
    std::string fftw_wisdom_dir;
    std::string outputdir;
    bool return_peak_to_root;
    bool saveccvols;
    bool resume;
    bool force_files_exist;
    int nice;

    ConfigDataClient() { }
    ConfigDataClient(const tom::corr::ConfigDataServer &c);

    int assert_status() const {
        return this->volsize > 0 &&
               this->binning >= 1 &&
               this->volsize / this->binning > 0 &&
               (!this->return_peak_to_root || this->cc_mask_radius > 0) &&
               !(this->outputdir.empty() && (this->saveccvols || !this->resume));
    }

    std::size_t get_volsize_effective() const { return this->volsize / this->binning; }
    void write_to_log(std::ostream &stream) const;
};






/****************************************************************************//**
 *
 * The job manager has two modes. \n
 * The configuration mode where it is possible to access the underlying vectors
 * directly as a non const reference. (through the methods access_ftemplates
 * etc.). It can be used, to configure (set) the job manager without assigning
 * a new vector directly and thus reducing the amount of copies to be done.\n
 * In working mode the call of these functions throws an exception. To switch
 * into working mode (by calling the method set_working_mode) the parameters
 * must be in a consistent state (otherwise an exceptio is thrown).
 * Executing the working functions in configuration mode on the other hand throws
 * an exception, too.
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
class JobManagerClient {

public:

    struct filename_triple {
        std::string filename_ccvol;
        std::string filename_ccidx;
        std::string filename_angle;
    };

    enum {
        STATUS_CONFIG = 1,
        STATUS_PROCESS = 2,
        STATUS_FINISHED = 4
    };


    JobManagerClient(const tom::corr::ConfigDataClient &c);

    void set_config_mode() { this->status = STATUS_CONFIG; }
    void set_working_mode();

    bool is_in_config_mode()  { return this->status==STATUS_CONFIG; }
    bool is_in_working_mode() { return this->status!=STATUS_CONFIG; }

    const std::vector<std::string> &get_ftemplates() const { if (this->status==STATUS_CONFIG) { throw std::runtime_error("Access to ftemplates in configuration mode."); } return this->ftemplates; }
    const std::vector<std::string> &get_fparticles() const { if (this->status==STATUS_CONFIG) { throw std::runtime_error("Access to fparticles in configuration mode."); } return this->fparticles; }
    const tom::Volume<double> *get_angles() const { if (this->status==STATUS_CONFIG) { throw std::runtime_error("Access to angles in configuration mode."); } return this->angles.get(); }
    std::vector<std::string>                             &access_ftemplates()      { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to ftemplates in working mode.");      } return this->ftemplates;      }
    std::vector<std::string>                             &access_fparticles()      { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to fparticles in working mode.");      } return this->fparticles;      }
    std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > &access_wedge_templates() { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to wedge_templates in working mode."); } return this->wedge_templates; }
    std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > &access_wedge_particles() { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to wedge_particles in working mode."); } return this->wedge_particles; }
    std::vector<std::size_t>                             &access_templates_idx()   { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to templates_idx in working mode.");   } return this->templates_idx;   }
    std::vector<std::size_t>                             &access_particles_idx()   { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to particles_idx in working mode.");   } return this->particles_idx;   }
    std::auto_ptr<tom::Volume<double> >                  &access_angles()          { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to angles in working mode.");          } return this->angles;          }
    std::vector<char>                                    &access_process_pair()    { if (this->status!=STATUS_CONFIG) { throw std::runtime_error("Access to process_pair in working mode.");    } return this->process_pair;    }

    int assert_status() const;

    const std::vector<tom::cc_peak<TFLOAT> > &get_peak_list() const;

    bool client_has_work() const    { return this->status!=STATUS_CONFIG   && this->ntemplates>0; }
    bool client_has_results() const { return this->status==STATUS_FINISHED && !this->peak_list.empty(); }

    std::size_t get_volsize_effective() const { return this->volsize_effective; }

    void process();

private:

    void init_effective_jobs();



    void do_corr3d(tom::CorrelationHandler<TFLOAT> *corrhdl) const;


    bool check_existing(std::size_t itemplate, std::size_t iparticle, tom::cc_peak<TFLOAT> *ipeak) const;

    std::string outputdir;
    bool resume;
    bool saveccvols;
    bool get_peaks;
    int status;
    std::size_t volsize;
    std::size_t binning;
    std::size_t volsize_effective;
    std::auto_ptr<const tom::Volume<TMASKCC> > maskcc;
    unsigned fftw_flag;
    bool reapply_mask_after_wedge;
    bool force_files_exist;
    std::auto_ptr<const tom::Volume<TFLOAT> > mask_sphere;

    std::vector<std::string> ftemplates;
    std::vector<std::string> fparticles;
    std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > wedge_templates;
    std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > wedge_particles;
    std::vector<std::size_t> templates_idx;
    std::vector<std::size_t> particles_idx;
    std::auto_ptr<tom::Volume<double> > angles;
    std::vector<char> process_pair;

    std::size_t ntemplates, nparticles;


    std::vector<std::string> job_ftemplates;
    std::vector<std::string> job_fparticles;
    std::vector<tom::Wedge<TFLOAT> *> job_wedge_templates;
    std::vector<tom::Wedge<TFLOAT> *> job_wedge_particles;
    std::vector<std::size_t> job_templates_idx;
    std::vector<std::size_t> job_particles_idx;
    std::vector<char> job_process_pair;

    std::size_t job_ntemplates, job_nparticles;



    std::vector<typename tom::corr::JobManagerClient<TFLOAT, TMASKCC>::filename_triple > output_filenames;
    std::vector<tom::cc_peak<TFLOAT> > peak_list;

};



}
}



// INLINE FUNCTIONS


#include "tom/corr/jobmanager_server.hpp"


/****************************************************************************//**
 *
 *******************************************************************************/
inline tom::corr::ConfigDataClient::ConfigDataClient(const tom::corr::ConfigDataServer &c)
    : volsize(c.volsize),
        binning(c.binning),
        sphere_mask_inner_radius(c.sphere_mask_inner_radius),
        sphere_mask_sigma(c.sphere_mask_sigma),
        sphere_mask_cutoff_radius(c.sphere_mask_cutoff_radius),
        cc_mask_radius(c.cc_mask_radius),
        fftw_flag(c.fftw_flag),
        reapply_mask_after_wedge(c.reapply_mask_after_wedge),
        fftw_wisdom_dir(c.fftw_wisdom_dir),
        outputdir(c.outputdir),
        return_peak_to_root(!c.peakfilename.empty()),
        saveccvols(c.saveccvols),
        resume(c.resume),
        force_files_exist(c.force_files_exist),
        nice(c.nice) {
    if (this->binning < 1) {
        this->binning = 1;
    }
    assert(c.assert_status());
    assert(this->assert_status());
}





#endif


