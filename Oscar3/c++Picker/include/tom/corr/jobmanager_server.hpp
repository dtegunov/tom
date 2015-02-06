/****************************************************************************//**
 * \file jobmanager_server.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    11.01.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__JOBMANAGER_SERVER_HPP__
#define ___INCLUDE_CORR__JOBMANAGER_SERVER_HPP__



#include <map>
#include <string>
#include <vector>
#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "tom/core/wedge.hpp"


#define TOM_MPI_MAIN__JOB_NOT_ASSIGNED    int(1)
#define TOM_MPI_MAIN__JOB_ASSIGNED        int(2)
#define TOM_MPI_MAIN__JOB_FINISHED        int(4)



namespace tom {
namespace corr {


class ConfigDataClient;
class ConfigDataServer;


/****************************************************************************//**
 *
 *
 *******************************************************************************/
class ConfigDataServer {

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
    unsigned ntemplates_amount;
    unsigned nparticles_amount;
    std::string peakfilename;
    std::string outputdir;
    std::string logfile;
    bool saveccvols;
    bool resume;
    bool force_files_exist;
    int nice;

    int assert_status() const {
        return this->volsize > 0 &&
               this->binning >= 1 &&
               this->volsize / this->binning > 0 &&
               (this->peakfilename.empty() || this->cc_mask_radius > 0) &&
               !(this->outputdir.empty() && (this->saveccvols || !this->resume));
    }

    void write_to_log(std::ostream &stream) const;
    void write_to_conf(const std::string &filename) const;
    void read_from_conf(const std::string &filename);
};





/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
class JobManagerServer {

public:
    void init_job(const std::vector<std::string> &ftemplates,
                const std::vector<std::string> &fparticles,
                const std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > &wedge_templates,
                const std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > &wedge_particles,
                const std::vector<const tom::Volume<double> *> &angles);

    std::vector<std::pair<std::size_t, std::size_t> > get_jobs_in_state(int type) const;
    std::size_t get_number_of_jobs_in_state(int type) const;

    bool get_assigned_jobs(int rank, std::vector<std::size_t> &templates_idx, std::vector<std::size_t> &particles_idx, std::vector<char> &process_pair) const ;
    bool get_assigned_jobs(int rank,
                        std::vector<std::string> &ftemplates,
                        std::vector<std::string> &fparticles,
                        std::vector<std::size_t> &templates_idx,
                        std::vector<std::size_t> &particles_idx,
                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_templates,
                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_particles,
                        const tom::Volume<double> *&angle,
                        std::vector<char> &process_pair) const;



    std::pair<int,int> get_job_state(std::size_t itemplate, std::size_t iparticle) const;

    std::size_t get_number_pending_jobs() const { return this->ntemplates*this->nparticles - this->number_finished_jobs; };
    std::size_t get_number_assigned_jobs() const { return this->number_assigned_jobs; };
    std::size_t get_number_finished_jobs() const { return this->number_finished_jobs; };
    std::size_t get_number_jobs_not_yet_assigned() const { return this->ntemplates*this->nparticles - this->number_finished_jobs - this->number_assigned_jobs; };
    std::size_t get_number_remaining_angles_tags() const { return this->angles_not_yet_assigned.size(); }

    void mark_jobs_as_done(int rank);
    bool assign_new_job(int rank,
                        std::size_t ntemplates_amount, std::size_t nparticles_amount,
                        std::vector<std::string> &ftemplates,
                        std::vector<std::string> &fparticles,
                        std::vector<std::size_t> &templates_idx,
                        std::vector<std::size_t> &particles_idx,
                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_templates,
                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_particles,
                        const tom::Volume<double> *&angle,
                        std::vector<char> &process_pair);

    int assert_status() const;
    void print() const;
    void print_jobs(int rank) const;

    const std::vector<std::string> &get_ftemplates() const { return this->ftemplates; }
    const std::vector<std::string> &get_fparticles() const { return this->fparticles; }
    const std::vector<std::pair<int, boost::shared_ptr<tom::Volume<double> > > > &get_angles() const { return this->angles; };


private:
    void mark_jobs_as_done(int rank, bool erase_value);

    std::size_t ntemplates, nparticles;
    std::vector<std::string> ftemplates;
    std::vector<std::string> fparticles;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > wedge_templates;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > wedge_particles;
    std::vector<std::pair<int, boost::shared_ptr<tom::Volume<double> > > > angles;
    std::size_t number_assigned_jobs;
    std::size_t number_finished_jobs;
    std::map<int, std::vector<std::size_t> > angles_not_yet_assigned;
    std::vector<std::pair<int, int> > job_assignment;
    std::map<int, std::vector<std::size_t> > job_currently_assigned;
};


void tag_angles(const std::vector<const tom::Volume<double> *> &input, std::vector<std::pair<int, boost::shared_ptr<tom::Volume<double> > > > &angles);


}
}




// INLINE FUNCTIONS...


#include "tom/corr/jobmanager_client.hpp"




/****************************************************************************//**
 * \brief Marks all jobs of a process as done.
 *
 * \param[in] rank The id of of the process.
 *
 * Marks the jobs of the process as finished.
 *******************************************************************************/
template<typename TFLOAT>
inline void tom::corr::JobManagerServer<TFLOAT>::mark_jobs_as_done(int rank) {
    this->mark_jobs_as_done(rank, true);
}






#endif




