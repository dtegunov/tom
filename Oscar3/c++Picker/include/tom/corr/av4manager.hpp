/****************************************************************************//**
 * \file av4manager.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    14.02.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__AV4MANAGER_HPP__
#define ___INCLUDE_CORR__AV4MANAGER_HPP__


#include <string>
#include <vector>
#include <map>
#include <assert.h>

#include <boost/shared_ptr.hpp>

#include "tom/core/volume.hpp"
#include "tom/core/wedge.hpp"





class GetPot; // forward declaration.



namespace tom {


namespace avg {
struct st_average_result; // forward declaration.
} // namespace avg


namespace mpi {
class exec_params; // forward declaration.
} // namespace mpi




namespace av4 {

namespace nt {
    enum e_next_template { AVG_N, AVG_1, AVG_SQRTN, AVG_NONE, AVG_S };

    e_next_template str2e_next_template(const std::string &e);

    std::string e_next_template2str(e_next_template e);

} // namespace nt





/****************************************************************************//**
 *
 *******************************************************************************/
class manager {
public:
    typedef double TFLOAT;

    static manager *create(const std::string &f, std::ostream &clog);
    virtual ~manager();
    virtual std::string getConfigStr() const;

    void process();
    static std::map<std::string, std::string> get_parameter_description();


protected:
    manager(GetPot &f, std::ostream &clog);

    void run_average();
    void run_correlation();
    void create_correlation_input();
    void create_average_input();
    virtual bool generateNextRun() = 0;



    // Parameters in config file.
    std::string c_outputdir;
    std::string c_templates;
    std::string c_particles;
    std::string c_prefix;
    std::size_t c_skip_iterations;
    std::size_t c_max_iteration;
    std::size_t c_nice;
    std::string c_fftw_wisdom_dir;
    std::size_t c_np;
    std::size_t c_sleep_before_alignment;
    unsigned    c_fftw_flag;
    bool        c_force_files_exist;
    double      c_sphere_mask_inner_radius;
    double      c_sphere_mask_sigma;
    double      c_sphere_mask_cutoff_radius;
    double      c_drop_template_wedge_sphere_radius; // If the template is a binary wedge, with radius at least of c_drop_template_wedge_sphere_radius it is ignored. Zero, means never drop.
    double      c_average_threshold;
    bool        c_saveccvols;
    bool        c_resume;
    std::string c_correlation_exe;
    std::string c_average_exe;
    std::string c_mpi_command;

    // Paramaters in configfile, but they may change for each iteration.
    // Although they may change depending on the implementation of
    // generateNextRun they appear in every implementation. Thus the prefix "c0_"
    double      c0_cc_mask_radius;
    bool        c0_reapply_mask_after_wedge;

    std::ostream &clog;

    std::size_t volsize; // The volume size stays constant. Thus it is read at startup from the first particle.

    std::vector<std::string> ftemplates0;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_templates0;
    std::vector<std::string> fparticles0;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_particles0;

    tom::mpi::exec_params *mpi_command;

    struct run_config {
        unsigned binning;
        double cc_mask_radius;
        bool reapply_mask_after_wedge;
        unsigned ntemplates_amount;
        unsigned nparticles_amount;
        std::vector<std::string> ftemplates;
        std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_templates;
        std::vector<std::string> fparticles;
        // There is a different wedge for the metric and for the averaging of the particles. (one to used in corrpar, one in average).
        std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_particles_metric;
        std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_particles_avg;
        std::vector<boost::shared_ptr<tom::Volume<double> > > angles;
        std::vector<std::pair<int, int> > particle_shift;
        bool do_fsc;

        void log(std::ostream &clog) const;
    };
    std::vector<boost::shared_ptr<run_config> > rc;

private:
    manager(const manager &m): clog(m.clog) { assert(0); }; // HIDE copy constructor.
    manager &operator=(const manager &m)    { assert(0); return *this; }; // HIDE assignemnt operator.
}; // class tom::av4::manager






/****************************************************************************//**
 *
 *******************************************************************************/
class manager_static: public manager {
public:
    manager_static(GetPot &f, std::ostream &clog);
    virtual ~manager_static();
    virtual std::string getConfigStr() const;
    static std::map<std::string, std::string> get_parameter_description();
protected:

private:
    virtual bool generateNextRun();

    // Parameters in config file.
    std::size_t c_binning;
    std::string c_angles_list;
    std::string c_wedge_particles_metric;
    bool c_ignore_particle_wedge;
    tom::av4::nt::e_next_template c_next_template;


    std::auto_ptr<tom::Volume<double> > angles_list;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_particles_metric;

}; // tom::av4::manager_static






/****************************************************************************//**
 *
 *******************************************************************************/
class manager_angular_refinement: public manager {
public:
    manager_angular_refinement(GetPot &f, std::ostream &clog);
    virtual ~manager_angular_refinement();
    virtual std::string getConfigStr() const;
    static std::map<std::string, std::string> get_parameter_description();
protected:

private:
    virtual bool generateNextRun();

    // Parameters in config file.
    std::string c_angles_list;
    std::string c_wedge_particles_metric;
    bool c_ignore_particle_wedge;
    tom::av4::nt::e_next_template c_next_template;

    std::size_t c_break_max_iteration_list;
    std::size_t c_break_max_iteration_refinement;

    double c_break_max_cc_list;
    double c_break_max_cc_refinement;

    double c_break_min_cc_change_list;
    double c_break_min_cc_change_refinement;

    double      c_sampling_step_list;
    std::size_t c_sampling_step_refinement;

    std::size_t c_binning;



    std::auto_ptr<tom::Volume<double> > angles_list;
    std::vector<boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> > > wedge_particles_metric;

    enum e_ang_ref_status { GENERAL_SAMPLING, REFINEMENT };
    double cc_min_1, cc_min_2;
    e_ang_ref_status ang_ref_status;
    double angular_sampling_step;
    std::size_t current_binning;

    std::size_t iter_rc_idx0;
    bool check_iteration_terminated() const;


}; // tom::av4::manager_static





} // namespace av4
} // namespace tom






#endif



