/****************************************************************************//**
 * \file av4manager.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    14.02.2008
 *******************************************************************************/
#include "tom/corr/av4manager.hpp"


#include <fstream>
#include <sstream>
#include <stdexcept>
#include <fftw3.h>
#include <algorithm>
#include <cctype>
#include <iomanip>

#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>




#include <helper/filesystem.hpp>
#include <helper/GetPot>
#include <helper/snippets.hpp>
#include <helper/exec.hpp>

#include <tom/core/io.h>
#include <tom/core/volume.hpp>

#include <tom/corr/config_files.hpp>
#include <tom/corr/jobmanager_server.hpp>
#include <tom/corr/filename_generators.hpp>

#include <tom/mpi/exec_params.hpp>







/****************************************************************************//**
 *
 *******************************************************************************/
tom::av4::manager *tom::av4::manager::create(const std::string &filename, std::ostream &clog) {

    if (!helper::fs::is_regular(filename)) {
        throw std::runtime_error("Parsing error: Could not open \"" + filename + "\".");
    }

    GetPot f(filename.c_str(), 0x0, 0x0, " \t");

    std::auto_ptr<tom::av4::manager> m;

    std::string s = f("refinement", "");
    std::transform(s.begin(), s.end(), s.begin(), toupper);

    if (s == "STATIC") {
        m.reset(new manager_static(f, clog));
    } else if (s == "ANGULAR_REFINEMENT") {
        m.reset(new manager_angular_refinement(f, clog));
    } else {
        if (s.empty()) {
            throw std::runtime_error("The parameter refinement is missing.");
        } else {
            throw std::runtime_error("The refinement \"" + s + "\" is not valid.");
        }
    }

    // If there are other fields, throw an error.
    std::vector<std::string> ufos = f.unidentified_variables();
    if (!ufos.empty()) {
        std::ostringstream ss; ss << "Parsing error: " << ufos.size() << " unrecognized options in configfile: ";
        std::size_t i, max_ufos = ufos.size() - 1;
        for (i=0; i<max_ufos; i++) {
            ss << "\"" << ufos.at(i) << "\", ";
        }
        ss << "\"" << ufos.at(i) << "\".";
        throw std::runtime_error(ss.str());
    }

    return m.release();
}


namespace {
/****************************************************************************//**
 *
 *******************************************************************************/
inline std::string get_description(const std::map<std::string, std::string> &m, const std::string &key) {
    const std::map<std::string, std::string>::const_iterator mit = m.find(key);
    if (mit != m.end()) {
        return mit->second;
    }
    return "(NO DESCRIPTION FOR " + key + ")";
}
}

/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::av4::manager_static::getConfigStr() const {

    std::map<std::string, std::string> parameter_description;
    this->get_parameter_description().swap(parameter_description);

    const std::size_t nlength = 80;
    const std::string prefix("# ");
    std::ostringstream ss;

    ss << this->manager::getConfigStr() <<
          ""                                                                    "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "refinement"), prefix, nlength);
    ss << "refinement=STATIC"                                                   "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "ignore_particle_wedge"), prefix, nlength);
    ss << "ignore_particle_wedge=" << c_ignore_particle_wedge <<                "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "next_template"), prefix, nlength);
    ss << "next_template=" << e_next_template2str(c_next_template) <<           "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "binning"), prefix, nlength);
    ss << "binning=" << c_binning <<                                            "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "wedge_particles_metric"), prefix, nlength);
    ss << "wedge_particles_metric=" << c_wedge_particles_metric <<              "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "angles_list"), prefix, nlength);
    ss << "angles_list=" << c_angles_list <<                                    std::endl;
    return ss.str();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::av4::manager_angular_refinement::getConfigStr() const {

    std::map<std::string, std::string> parameter_description;
    this->get_parameter_description().swap(parameter_description);

    const std::size_t nlength = 80;
    const std::string prefix("# ");
    std::ostringstream ss;
    ss << this->manager::getConfigStr() <<
          ""                                                                    "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "refinement"), prefix, nlength);
    ss << "refinement=ANGULAR_REFINEMENT"                                       "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "angles_list"), prefix, nlength);
    ss << "angles_list=" << c_angles_list <<                                    "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "ignore_particle_wedge"), prefix, nlength);
    ss << "ignore_particle_wedge=" << c_ignore_particle_wedge <<                "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "wedge_particles_metric"), prefix, nlength);
    ss << "wedge_particles_metric=" << c_wedge_particles_metric <<              "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "next_template"), prefix, nlength);
    ss << "next_template=" << e_next_template2str(c_next_template) <<           "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "break_max_iteration_list"), prefix, nlength);
    ss << "break_max_iteration_list=" << c_break_max_iteration_list <<          "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "break_max_iteration_refinement"), prefix, nlength);
    ss << "break_max_iteration_refinement=" << c_break_max_iteration_refinement << "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "break_max_cc_list"), prefix, nlength);
    ss << "break_max_cc_list=" << c_break_max_cc_list <<                        "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "break_max_cc_refinement"), prefix, nlength);
    ss << "break_max_cc_refinement=" << c_break_max_cc_refinement <<            "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "break_min_cc_change_list"), prefix, nlength);
    ss << "break_min_cc_change_list=" << c_break_min_cc_change_list <<          "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "break_min_cc_change_refinement"), prefix, nlength);
    ss << "break_min_cc_change_refinement=" << c_break_min_cc_change_refinement << "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "sampling_step_list"), prefix, nlength);
    ss << "sampling_step_list=" << c_sampling_step_list <<                      "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "sampling_step_refinement"), prefix, nlength);
    ss << "sampling_step_refinement=" << c_sampling_step_refinement<<           "\n" "\n";
          helper::print_parameters(ss, "", ::get_description(parameter_description, "binning"), prefix, nlength);
    ss << "binning=" << c_binning <<                                            "\n";
    return ss.str();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::av4::manager::getConfigStr() const {

    std::map<std::string, std::string> parameter_description;
    this->get_parameter_description().swap(parameter_description);

    const std::size_t nlength = 80;
    const std::string prefix("# ");

    std::ostringstream ss;
    helper::print_parameters(ss, "", ::get_description(parameter_description, "outputdir"), prefix, nlength);
    ss << "outputdir=" << c_outputdir <<                                      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "templates"), prefix, nlength);
    ss << "templates=" << c_templates <<                                      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "particles"), prefix, nlength);
    ss << "particles=" << c_particles <<                                      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "prefix"), prefix, nlength);
    ss << "prefix=" << c_prefix <<                                            "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "skip_iterations"), prefix, nlength);
    ss << "skip_iterations=" << c_skip_iterations <<                           "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "max_iteration"), prefix, nlength);
    ss << "max_iteration=" << c_max_iteration <<                              "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "nice"), prefix, nlength);
    ss << "nice=" << c_nice <<                                                "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "fftw_wisdom_dir"), prefix, nlength);
    ss << "fftw_wisdom_dir=" << c_fftw_wisdom_dir <<                          "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "np"), prefix, nlength);
    ss << "np=" << c_np <<                                                    "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "fftw_flag"), prefix, nlength);
    ss << "fftw_flag=" << tom::fftw::flag2str(c_fftw_flag) <<                 "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "force_files_exist"), prefix, nlength);
    ss << "force_files_exist=" << (c_force_files_exist?1:0) <<                "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "sphere_mask_inner_radius"), prefix, nlength);
    ss << "sphere_mask_inner_radius=" << c_sphere_mask_inner_radius <<        "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "sphere_mask_sigma"), prefix, nlength);
    ss << "sphere_mask_sigma=" << c_sphere_mask_sigma <<                      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "sphere_mask_cutoff_radius"), prefix, nlength);
    ss << "sphere_mask_cutoff_radius=" << c_sphere_mask_cutoff_radius <<      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "drop_template_wedge_sphere_radius"), prefix, nlength);
    ss << "drop_template_wedge_sphere_radius=" << c_drop_template_wedge_sphere_radius << "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "saveccvols"), prefix, nlength);
    ss << "saveccvols=" << (c_saveccvols?1:0) <<                              "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "resume"), prefix, nlength);
    ss << "resume=" << (c_resume?1:0) <<                                      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "cc_mask_radius"), prefix, nlength);
    ss << "cc_mask_radius=" << c0_cc_mask_radius <<                           "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "reapply_mask_after_wedge"), prefix, nlength);
    ss << "reapply_mask_after_wedge=" << (c0_reapply_mask_after_wedge?1:0) << "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "average_threshold"), prefix, nlength);
    ss << "average_threshold=" << c_average_threshold <<                      "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "correlation_exe"), prefix, nlength);
    ss << "correlation_exe=" << c_correlation_exe <<                          "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "average_exe"), prefix, nlength);
    ss << "average_exe=" << c_average_exe <<                                  "\n" "\n";
    helper::print_parameters(ss, "", ::get_description(parameter_description, "mpi_command"), prefix, nlength);
    ss << "mpi_command=" << c_mpi_command <<                                  "\n";
            ;
    return ss.str();
}


/****************************************************************************//**
 *\brief Destructor.
 *******************************************************************************/
tom::av4::manager::~manager() {

    delete mpi_command;

}

/****************************************************************************//**
 *\brief Destructor.
 *******************************************************************************/
tom::av4::manager_angular_refinement::~manager_angular_refinement() {

}



/****************************************************************************//**
 *\brief Destructor.
 *******************************************************************************/
tom::av4::manager_static::~manager_static() {

}



/****************************************************************************//**
 *
 *******************************************************************************/
void tom::av4::manager::create_average_input() {

    // just check that the files are still there. Could also be omitted...

    const std::size_t iit = rc.size() - 1;

    tom::av4::filename_generator f_gen(c_outputdir, c_prefix);

    if (!helper::fs::is_regular(f_gen.get_config_correlation(iit) )) { throw std::runtime_error("The config file \""   + f_gen.get_config_correlation(iit)  + "\" disappeared unexpectedly."); }
    if (!helper::fs::is_regular(f_gen.get_templates_list(iit)     )) { throw std::runtime_error("The template file \"" + f_gen.get_templates_list(iit)      + "\" disappeared unexpectedly."); }
    if (!helper::fs::is_regular(f_gen.get_angles_list(iit)        )) { throw std::runtime_error("The angles file \""   + f_gen.get_angles_list(iit)         + "\" disappeared unexpectedly."); }

    const manager::run_config &rci = *rc.back();

    const std::size_t W = 31;

    // Create the particles list.
    std::string s = f_gen.get_particles_list_average(iit);
    clog << "# " << std::setw(W) << "particles_list_average = " << "\"" << s << "\"." << std::endl;
    tom::corr::write_filelist(s, rci.fparticles, rci.wedge_particles_avg);

}



/****************************************************************************//**
 * \brief Creates all the necessary files and directories for the next run.
 *******************************************************************************/
void tom::av4::manager::create_correlation_input() {

    assert(!rc.empty());
    assert(!rc.at(rc.size()-1)->ftemplates.empty() && rc.at(rc.size()-1)->ftemplates.size()==rc.at(rc.size()-1)->wedge_templates.size());
    assert(!rc.at(rc.size()-1)->fparticles.empty() && rc.at(rc.size()-1)->fparticles.size()==rc.at(rc.size()-1)->wedge_particles_avg.size() && rc.at(rc.size()-1)->fparticles.size()==rc.at(rc.size()-1)->wedge_particles_metric.size());

    tom::av4::filename_generator f_gen(c_outputdir, c_prefix);
    const std::size_t iit = rc.size() - 1;

    const std::size_t W = 31;
    std::string s;

    // Check and create the output directory.
    const std::string outputdir_iit = f_gen.get_outputdir(iit);
    bool created = false;
    if (!helper::fs::is_directory(outputdir_iit)) {
        if (!helper::fs::create_directory(outputdir_iit)) {
            clog << "# ERROR CREATING OUTPUTDIR: \"" << outputdir_iit << "\"." << std::endl;
            throw std::runtime_error("can not create output directory.");
        }
        created = true;
    }
    clog << "# " << std::setw(W) << "outputdir = " << "\"" << outputdir_iit << "\"" << (created?" (created)":"") << std::endl;

    const manager::run_config &rci = *rc.back();

    // Create the templates list.
    s = f_gen.get_templates_list(iit);
    clog << "# " << std::setw(W) << "templates_list = " << "\"" << s << "\"." << std::endl;
    tom::corr::write_filelist(s, rci.ftemplates, rci.wedge_templates);

    // Create the particles list.
    s = f_gen.get_particles_list_correlation(iit);
    clog << "# " << std::setw(W) << "particles_list = " << "\"" << s << "\"." << std::endl;
    tom::corr::write_filelist(s, rci.fparticles, rci.wedge_particles_metric);

    // Create the angles list.
    s = f_gen.get_angles_list(iit);
    clog << "# " << std::setw(W) << "angles_list = " << "\"" << s << "\"." << std::endl;
    tom::corr::write_angleslist(s, rci.ftemplates.size(), rci.fparticles.size(), rci.angles, f_gen, iit);

    // write the configuration file.
    tom::corr::ConfigDataServer cs;
    cs.volsize = volsize;
    cs.binning = rci.binning;
    cs.sphere_mask_inner_radius = c_sphere_mask_inner_radius;
    cs.sphere_mask_sigma = c_sphere_mask_sigma;
    cs.sphere_mask_cutoff_radius = c_sphere_mask_cutoff_radius;
    cs.cc_mask_radius = rci.cc_mask_radius;
    cs.fftw_flag = c_fftw_flag;
    cs.reapply_mask_after_wedge = rci.reapply_mask_after_wedge;
    cs.fftw_wisdom_dir = c_fftw_wisdom_dir;
    cs.ntemplates_amount = rci.ntemplates_amount;
    cs.nparticles_amount = rci.nparticles_amount;
    cs.peakfilename = f_gen.get_peakfilename(iit);
    cs.outputdir = outputdir_iit;
    cs.logfile = f_gen.get_logfilename_correlation(iit);
    cs.saveccvols = c_saveccvols;
    cs.resume = c_resume;
    cs.force_files_exist = c_force_files_exist;
    cs.nice = c_nice;


    s = f_gen.get_config_correlation(iit);
    clog << "# " << std::setw(W) << "corr_config = " << "\"" << s << "\"." << std::endl;
    cs.write_to_conf(s);

    //helper::fs::remove(cs.peakfilename);
}



namespace {
void log_exec_results(std::ostream &clog, int status, const std::vector<char> &stdout_, const std::vector<char> &stderr_, const std::string &prefix) {

    int ndata = (stdout_.empty()?0:1) + (stderr_.empty()?0:1);

    if (WIFEXITED(status)) {
        clog << prefix << "child process terminated normally with status " << WEXITSTATUS(status) << std::endl;
    } else if (WIFSIGNALED(status)) {
        clog << prefix << "child process  terminated due to signal " << WTERMSIG(status) << std::endl;
    } else {
        clog << prefix << "wait for child returned status " << status << ". Strange: this is neighter WIFEXITED nor WIFSIGNALED... ignore it..." << std::endl;
    }
    if (ndata) {
        clog << prefix << "output stream" << (ndata>1?"s":"") << " contain" << (ndata==1?"s":"") << " the following data..." << std::endl;
        if (!stdout_.empty()) {
            clog << prefix << "STDOUT(" << stdout_.size() << " bytes) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            clog.write(&stdout_[0], stdout_.size());
            clog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        }
        if (!stderr_.empty()) {
            clog << prefix << "STDERR(" << stderr_.size() << " bytes) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            clog.write(&stderr_[0], stderr_.size());
            clog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        }
    }
}

} // namespace


/****************************************************************************//**
 * \brief calls mpirun to start the parallel correlation program.
 *******************************************************************************/
void tom::av4::manager::run_correlation() {

    const std::size_t iit = rc.size() - 1;

    tom::av4::filename_generator f_gen(c_outputdir, c_prefix);

    // Create the command line arguments.
    std::string file = c_correlation_exe;
    std::vector<std::string> argv;
    argv.push_back(file);
    argv.push_back(f_gen.get_config_correlation(iit));
    argv.push_back(f_gen.get_templates_list(iit));
    argv.push_back(f_gen.get_particles_list_correlation(iit));
    argv.push_back(f_gen.get_angles_list(iit));



    // warp the call of mpirun around the parameters.
    mpi_command->execvp(file, argv);

    #if 0
    file = "ls";
    argv.resize(0);
    argv.push_back(file);
    argv.push_back("-l");
    argv.push_back("-a");
    #endif

    int status;
    helper::c_execv call_exec(file, argv, true);

    std::vector<char> stdout_, stderr_;

    clog << "# $ " << file;
    {
        std::vector<std::string>::const_iterator it = argv.begin();
        if (it != argv.end()) {
            it++;
            for (; it!=argv.end(); it++) {
                clog << " '" << *it << "'";
            }
        }
    }
    clog << std::endl;

    if (iit < c_skip_iterations) {
        clog << "#     skip the actual call of the correlation because of skip_iterations=" << c_skip_iterations << std::endl;
    } else {
        if (iit>0 && c_sleep_before_alignment) {
            sleep(c_sleep_before_alignment);
        }
        clog << "#     correlation started at " << helper::now2str();
        if (iit>0 && c_sleep_before_alignment) {
            clog << " (after sleeping for " << c_sleep_before_alignment << " second" << (c_sleep_before_alignment==1?"":"s") << ")";
        }
        clog << std::endl;
        try {
            helper::exec<const helper::c_execv &>(call_exec, stdout_, stderr_, &status);
        } catch (std::exception &e) {
            clog << "#     caught exception: " << e.what() << std::endl;
            ::log_exec_results(clog, status, stdout_, stderr_, "#     ");
            throw;
        } catch (...) {
            clog << "#     caught unknown exception" << std::endl;
            ::log_exec_results(clog, status, stdout_, stderr_, "#     ");
            throw;
        }
        clog << "#     finished at " << helper::now2str() << std::endl;
        ::log_exec_results(clog, status, stdout_, stderr_, "#     ");

        if (!WIFEXITED(status) ||WEXITSTATUS(status)!=0) {
            throw std::runtime_error("The child process terminated not with error status 0");
        }
    }

}


/****************************************************************************//**
 * \brief calls mpirun to start the parallel correlation program.
 *******************************************************************************/
void tom::av4::manager::run_average() {

    const std::size_t iit = rc.size() - 1;

    tom::av4::filename_generator f_gen(c_outputdir, c_prefix);


    // Create the command line arguments.
    std::string file = c_average_exe;

    // Wait for some time, because as the peakfile were created by an other process
    // on maybe an other time, the file, may not yet be available in the network.
    const std::string peakfilename = f_gen.get_peakfilename(iit);
    {
        bool peakfile_exist = false;
        const std::size_t nwait = 30;
        for (std::size_t i=0; i<nwait; i++) {
            if (!helper::fs::is_regular(peakfilename)) {
                if (i==0) {
                    clog << "# WAITING for peakfile " << peakfilename << std::endl;
                }
                sleep(1);
            } else {
                peakfile_exist = true;
                break;
            }
        }
        if (!peakfile_exist) {
            clog << "# Peakfile even after " << nwait << " seconds does not exist. Resume with the averaging, though it will probably fail..." << std::endl;
        }
    }

    std::vector<std::string> argv;
    argv.push_back(file);
    argv.push_back(f_gen.get_config_correlation(iit));
    argv.push_back(f_gen.get_templates_list(iit));
    argv.push_back(f_gen.get_particles_list_average(iit));
    argv.push_back(f_gen.get_angles_list(iit));
    argv.push_back(f_gen.get_logfilename_average(iit));
    argv.push_back(rc[iit]->do_fsc ? "yes" : "no");
    argv.push_back(boost::lexical_cast<std::string>(c_average_threshold));



    int status;
    helper::c_execv call_exec(file, argv, true);

    std::vector<char> stdout_, stderr_;

    clog << "# $ " << file;
    {
        std::vector<std::string>::const_iterator it = argv.begin();
        if (it != argv.end()) {
            it++;
            for (; it!=argv.end(); it++) {
                clog << " '" << *it << "'";
            }
        }
    }
    clog << std::endl;


    if (iit < c_skip_iterations) {
        clog << "#     skip the actual call of the averaging because of skip_iterations=" << c_skip_iterations << std::endl;
    } else {
        clog << "#     average started at " << helper::now2str() << std::endl;
        try {
            helper::exec<const helper::c_execv &>(call_exec, stdout_, stderr_, &status);
        } catch (std::exception &e) {
            clog << "#     caught exception: " << e.what() << std::endl;
            ::log_exec_results(clog, status, stdout_, stderr_, "#     ");
            throw;
        } catch (...) {
            clog << "#     caught unknown exception" << std::endl;
            ::log_exec_results(clog, status, stdout_, stderr_, "#     ");
            throw;
        }
        clog << "#     finished at " << helper::now2str() << std::endl;
        ::log_exec_results(clog, status, stdout_, stderr_, "#     ");

        if (!WIFEXITED(status) ||WEXITSTATUS(status)!=0) {
            throw std::runtime_error("The child process terminated not with error status 0");
        }
    }
}







/****************************************************************************//**
 *
 *******************************************************************************/
void tom::av4::manager::run_config::log(std::ostream &clog) const {

    const std::size_t W = 31;
    clog <<                                                                          "# " << std::setw(W) <<
            "binning = "                    << binning <<                       "\n" "# " << std::setw(W) <<
            "cc_mask_radius = "             << cc_mask_radius <<                "\n" "# " << std::setw(W) <<
            "reapply_mask_after_wedge = "   << reapply_mask_after_wedge <<      "\n" "# " << std::setw(W) <<
            "ntemplates_amount = "          << ntemplates_amount <<             "\n" "# " << std::setw(W) <<
            "nparticles_amount = "          << nparticles_amount <<             "\n" "# " << std::setw(W) <<
            "do_fsc = "                     << do_fsc <<                        "\n" "# " << std::setw(W) <<
            "#templates = "                 << ftemplates.size() <<             "\n" "# " << std::setw(W) <<
            "#particles = "                 << fparticles.size() <<             std::endl;

}


/****************************************************************************//**
 *
 *******************************************************************************/
void tom::av4::manager::process() {

    do {
        clog << "###############################################################\n"
                "# generate iteration #" << rc.size() << std::endl;

        if (!generateNextRun()) {
            clog << "###############################################################" << std::endl;
            clog << "# generateNextRun did not create a new run. END." << std::endl;
            break;
        }

        // Print current configuration to log.
        rc.back()->log(clog);

        create_correlation_input();
        run_correlation();

        create_average_input();
        run_average();

        if (c_max_iteration>0 && rc.size()>=c_max_iteration) {
            clog << "###############################################################" << std::endl;
            clog << "# Reached maximum number of iterations. END" << std::endl;
            break;
        }

    } while (true);


}

/****************************************************************************//**
 *
 *******************************************************************************/
std::map<std::string, std::string> tom::av4::manager::get_parameter_description() {
    std::map<std::string, std::string> parameter_description;

    parameter_description["cc_mask_radius"] = "Parameter cc_mask_radius for alignment: shifts larger then cc_mask_radius (relative to the original volume size) are ignored when searching the peak.";
    parameter_description["reapply_mask_after_wedge"] = "Parameter reapply_mask_after_wedge for alignment: after applying the spectral weighting function (missing wedge) to the volumes, the previously applied sphere mask (sphere_mask_inner_radius) is maybe slightly distorted. Reapply the mask as a binary mask. Is maybe more accurate, but also more time consuming.";
    parameter_description["outputdir"] = "Name of the directory where to save the ouput. Existing files may be overwritten or taken to resume.";
    parameter_description["correlation_exe"] = "Filename of the program for alignment (correlation).";
    parameter_description["average_exe"] = "Filename of the program for averaging.";
    parameter_description["mpi_command"] = "Select how to call the parallel environment. Currently supported is 'mpirun' and 'poe'.";
    parameter_description["templates"] = "Filename of the list of template with which the program starts (may change during iterating).";
    parameter_description["particles"] = "Filename of the list of particles with which the program starts (may change during iterating).";
    parameter_description["prefix"] = "Prefix to prepend to the name of the output files.";
    parameter_description["skip_iterations"] = "Do not call the alignment nor the averaging program for the first iterations. Used to resume a previously started run with the same parameters.";
    parameter_description["sleep_before_alignment"] = "Gives the number of seconds to wait before calling the parallel alignment program. Because of the network it can happen, that the parallel processes do not yet see the files. In that case, they would fail. Thus wait for sleep_before_alignment seconds every time before the alignment is started.";
    parameter_description["max_iteration"] = "Maximum number of iterations (0 means unlimited).";
    parameter_description["nice"] = "Renice the parallel processes for the alignment.";
    parameter_description["fftw_wisdom_dir"] = "Name of the directory where to save and read the fftw-wisdom from.";
    parameter_description["np"] = "Number of parallel processes for the alignment.";
    parameter_description["fftw_flag"] = "Flag for creation of the fftw-plan.";
    parameter_description["force_files_exist"] = "If true, error during reading of a file results in a fatal error.";
    parameter_description["sphere_mask_inner_radius"] = "The inner radius of the sphere mask which is used for alignment and averaging (relative to the original volume size).";
    parameter_description["sphere_mask_sigma"] = "The make a soft edge of the sphere mask used for alignment and averaging (relative to the original volume size, 0 makes a binary sphere).";
    parameter_description["sphere_mask_cutoff_radius"] = "The sphere mask used for alignment and averaging is set to zero outside of this radius (relative to the original volume size, in conjunction sphere_mask_sigma).";
    parameter_description["average_threshold"] = "The parameter WEDGE_THRESHOLD of the averaging: Frequencies corresponding to the summed wedge of the particles which are smaller then this relative threshold, are set to zero.";
    parameter_description["drop_template_wedge_sphere_radius"] = "Check that the wedge of the new average (template) is all equally 1 inside this radius. In that case, start the next iteration with NOWEDGE (because it would only be a lowpass filter; speeds up alignment) (NOT YET IMPLEMENTED!).";
    parameter_description["saveccvols"] = "Parameter saveccvols for averaging. Correlation volumes are saved to the output directory, for resuming or manual examination.";
    parameter_description["resume"] = "Parameter resume for averaging. If true, the averaging process(es) look if the correlation volume is already saved from a previous run and in case, read the peak from there.";

    return parameter_description;
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::map<std::string, std::string> tom::av4::manager_static::get_parameter_description() {
    std::map<std::string, std::string> parameter_description;

    //tom::av4::manager::get_parameter_description().swap(parameter_description);

    parameter_description["refinement"] = "STATIC: Iterative alignment with a unchanging list of rotations.";

    parameter_description["next_template"] = "During averaging several averages (different weighted in fourier space) are computed. Which one to take for the next iteration?.";
    parameter_description["binning"] = "Binning faktor. Says unchanged for each iteration (ATTENTION. for example a value of 4 means 'combine 4 voxels', i.e. use a binning of 2 = log2(4)).";
    parameter_description["wedge_particles_metric"] = "Use for correlation/alignment an other wedge than for averaging (where the missing wedge from 'particles' is taken).";
    parameter_description["angles_list"] = "The list of angles which are used for alignment. (em-file with 3xNx1 angles; phi, psi, theta, in radians)";
    parameter_description["ignore_particle_wedge"] = "If true, the particles are treated as they have no missing wedge (regardless what the lists in 'particles' and 'wedge_particles_metric' says).";


    return parameter_description;
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::map<std::string, std::string> tom::av4::manager_angular_refinement::get_parameter_description() {
    std::map<std::string, std::string> parameter_description;

    //tom::av4::manager::get_parameter_description().swap(parameter_description);

    parameter_description["sampling_step_list"] = "This is the angluar increment to which the angles in the em-file 'angles_list' correspond (in degrees). It is the starting angle for the refinement.";
    parameter_description["binning"] = "Use this binning in the first step (0 or 1 means no binning). For each refinement steps the binning is reduced (ATTENTION. for example a value of 4 means 'combine 4 voxels', i.e. use a binning of 2 = log2(4)).";
    parameter_description["sampling_step_refinement"] = "For each refinement step for the rotations, the next angle-increment is the old one, divided by this faktor, where the first refinement step starts with 'sampling_step_list'.";
    parameter_description["break_min_cc_change_refinement"] = "If the correlation of all the current averages with their averages from the previous iteration does change less then this value (in percent), break the iteration for the current refinement (-1 means never break upon this condition).";
    parameter_description["break_min_cc_change_list"] = "If the correlation of all the current averages with their averages from the previous iteration does change less then this value (in percent), break the iteration for the first step (-1 means never break upon this condition).";
    parameter_description["break_max_cc_refinement"] = "If the correlation of all the new averages with their average from the previous  run, is larger, break the iteration for the current refinement iteration.";
    parameter_description["break_max_cc_list"] = "If the correlation of all the new averages with their average from the privious run, is larger, break the iteration for the first step.";
    parameter_description["break_max_iteration_refinement"] = "During the next steps, where the angles are refined, this is the maximum number of iterations (0 means endlessly).";
    parameter_description["break_max_iteration_list"] = "During the first step, where the angles list is used for sampling the alignment, this is the maximum number of iterations (0 means endlessly).";
    parameter_description["next_template"] = "During averaging several averages (different weighted in fourier space) are computed. Which one to take for the next iteration?.";
    parameter_description["ignore_particle_wedge"] = "If true, the particles are treated as they have no missing wedge (regardless what the lists in 'particles' and 'wedge_particles_metric' says).";
    parameter_description["wedge_particles_metric"] = "Use for correlation/alignment an other wedge than for averaging (where the missing wedge from 'particles' is taken; empty defaults to the ones given in 'particle').";
    parameter_description["refinement"] = "ANGULAR_REFINEMENT: First make an angular sampling as given by a input list. Then make angular refinement (Each step is repeated several times until one of the termination criterias is fullfiled).";
    parameter_description["angles_list"] = "The list of angles which are used for alignment in the first step. (em-file with 3xNx1 angles; phi, psi, theta, in radians)";

    return parameter_description;
}


/****************************************************************************//**
 *
 *******************************************************************************/
tom::av4::manager::manager(GetPot &f, std::ostream &clog)
    : clog(clog),
      rc() {

    std::ostringstream ss;
    std::string s;
    int i;

    c_outputdir = f("outputdir", "");
        if (c_outputdir.empty()) { throw std::runtime_error("Mandatory parameter \"outputdir\" is missing or empty."); }
        // Check that the output directory exists.
        helper::fs::path p(this->c_outputdir);
        if (helper::fs::exists(p)) {
            if (!helper::fs::is_directory(p)) {
                ss << "The output dir \"" << c_outputdir << "\" does not name a directory."; throw std::runtime_error(ss.str());
            }
        } else {
            ss << "The output dir \"" << c_outputdir << "\" does not exist."; throw std::runtime_error(ss.str());
        }
    c_templates = f("templates", "");
        if (c_templates.empty()) { throw std::runtime_error("Mandatory parameter \"templates\" is missing or empty."); }
        tom::corr::parse_filelist<manager::TFLOAT>(c_templates, ftemplates0, wedge_templates0, 0);
        if (ftemplates0.empty()) {
            throw std::runtime_error("The template list \"" + c_templates + "\" does not contain any templates.");
        }
    c_particles = f("particles", "");
        if (c_particles.empty()) { throw std::runtime_error("Mandatory parameter \"particles\" is missing or empty."); }
        tom::corr::parse_filelist<manager::TFLOAT>(c_particles, fparticles0, wedge_particles0, 0);
        if (fparticles0.empty()) {
            throw std::runtime_error("The particle list \"" + c_particles + "\" does not contain any particles.");
        } else {
            tom_io_em_header header;
            int i = ::tom_io_em_read_header(fparticles0[0].c_str(), NULL, &header, NULL);
            if (i != TOM_ERR_OK) {
                throw std::runtime_error("Opening the first particle \"" + fparticles0[0] + "\" to get the volume size failed.");
            }
            if (header.dims[0]!=header.dims[1] || header.dims[0]!=header.dims[2]) {
                ss << "The first particle \"" << fparticles0[0] << "\" is not cubic (has size " << header.dims[0] << "x" << header.dims[1] << "x" << header.dims[2] << ").";
                throw std::runtime_error(ss.str());
            }
            volsize = header.dims[0];
        }
    c_prefix = f("prefix", "");
        s = c_prefix;
        std::transform(s.begin(), s.end(), s.begin(), toupper);
        if (s.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") != std::string::npos) {
            throw std::runtime_error("The parameter prefix should contain only alphanumeric symbols.");
        }
    c_skip_iterations = ((i=f("skip_iterations", -1)) > 0) ? i : 0;
    c_max_iteration = ((i=f("max_iteration", -1))<0) ? 0 : i;
    c_sleep_before_alignment = (i=f("sleep_before_alignment", 1))<0 ? 0 : i;
    c_nice = f("nice", 0);
    c_fftw_wisdom_dir = f("fftw_wisdom_dir", c_outputdir.c_str());
    c_np = (i=f("np", 32));
        if (i < 2) { throw std::runtime_error("The parameter \"np\" must be >= 2."); }
    c_fftw_flag                         = tom::fftw::str2flag(f("fftw_flag", "FFTW_MEASURE"));
    c_force_files_exist                 = 0!=f("force_files_exist", 1);
    c_sphere_mask_inner_radius          = f("sphere_mask_inner_radius", ceil(volsize/2.) - 1.1);
    c_sphere_mask_sigma                 = f("sphere_mask_sigma", 0);
    c_sphere_mask_cutoff_radius         = f("sphere_mask_cutoff_radius", 0);
    c_average_threshold                 = f("average_threshold", 0.01);
    c_drop_template_wedge_sphere_radius = std::max(1., f("drop_template_wedge_sphere_radius", static_cast<double>(volsize/2)));
    c_saveccvols                        = 0!=f("saveccvols", 1);
    c_resume                            = 0!=f("resume", 1);

    c_correlation_exe = f("correlation_exe", "");
        if (c_correlation_exe.empty()) { throw std::runtime_error("Mandatory parameter \"correlation_exe\" is missing or empty."); }
        if (!helper::fs::is_regular(c_correlation_exe)) { throw std::runtime_error("\"" + c_correlation_exe + "\" does not exist."); }

    c_average_exe = f("average_exe", "");
        if (c_average_exe.empty()) { throw std::runtime_error("Mandatory parameter \"average_exe\" is missing or empty."); }
        if (!helper::fs::is_regular(c_average_exe)) { throw std::runtime_error("\"" + c_average_exe + "\" does not exist."); }

    c0_cc_mask_radius                   = f("cc_mask_radius", (ceil(volsize/2.) - 1.1) * 0.5);
    c0_reapply_mask_after_wedge         = 0!=f("reapply_mask_after_wedge", 1);

    #ifndef _AIX
    c_mpi_command = f("mpi_command", "mpirun");
    #else
    // defaults to poe on AIX
    c_mpi_command = f("mpi_command", "poe");
    #endif
        std::auto_ptr<tom::mpi::exec_params> mpi_command_ = create_from_config_parameter(c_mpi_command, c_np);
        assert(mpi_command_.get());



    mpi_command = mpi_command_.release();

}




/****************************************************************************//**
 *
 *******************************************************************************/
tom::av4::nt::e_next_template tom::av4::nt::str2e_next_template(const std::string &e) {
    std::string s(e);
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    if (s == "AVG_N") {
        return tom::av4::nt::AVG_N;
    } else if (s == "AVG_1") {
        return tom::av4::nt::AVG_1;
    } else if (s == "AVG_SQRTN") {
        return tom::av4::nt::AVG_SQRTN;
    } else if (s == "AVG_S") {
        return tom::av4::nt::AVG_S;
    } else if (s == "AVG_NONE") {
        return tom::av4::nt::AVG_NONE;
    } else {
        throw std::domain_error("str2e_next_template gets an unexpected string value.");
    }
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::av4::nt::e_next_template2str(tom::av4::nt::e_next_template  e) {
    switch (e) {
        case tom::av4::nt::AVG_N:
            return "AVG_N";
        case tom::av4::nt::AVG_1:
            return "AVG_1";
        case tom::av4::nt::AVG_SQRTN:
            return "AVG_SQRTN";
        case tom::av4::nt::AVG_S:
            return "AVG_S";
        case tom::av4::nt::AVG_NONE:
            return "AVG_NONE";
        default:
            throw std::domain_error("e_next_template2str gets an unexpected enum value.");
    }
}


/****************************************************************************//**
 * \brief The constructor of manager_static
 *
 * \param[in] f The config parameters to initialise the class.
 * \param[in,out] clog The output stream for the logging information of the
 *    class.
 *******************************************************************************/
tom::av4::manager_static::manager_static(GetPot &f, std::ostream &clog):
    manager(f, clog),
    c_angles_list(),
    c_ignore_particle_wedge(false),
    c_next_template(tom::av4::nt::AVG_1),
    angles_list(NULL) {



    std::string s;
    c_angles_list          = f("angles_list", "");
        if (c_angles_list.empty()) {
            throw std::runtime_error("the refinement STATIC needs the parameter \"angles_list\" with the name of an em-file of angles (phi,psi,theta in rad).");
        } else if (!helper::fs::is_regular(c_angles_list)) {
            throw std::runtime_error("the angle list \"" + c_angles_list + "\" must be the name of and em-file of angles (phi,psi,theta in rad), but the file does not exist.");
        }
        try {
            tom::Volume<double> *p;
            tom::read_from_em(p, c_angles_list, NULL,NULL,NULL, NULL,NULL);
            angles_list.reset(p);
        } catch (...) {
            throw std::runtime_error("the angle list \"" + c_angles_list + "\" could not be opened. It must be an em-file of angles (phi,psi,theta in rad).");
        }

        if (angles_list->getSizeX()!=3 || angles_list->getSizeZ()!=1) {
            throw std::runtime_error("the angle list \"" + c_angles_list + "\" must be an em-file of angles (phi,psi,theta in rad, 3xNx1).");
        }

    c_ignore_particle_wedge = f("ignore_particle_wedge", 0) != 0;

    c_wedge_particles_metric = f("wedge_particles_metric", "");
        if (c_wedge_particles_metric.empty()) {
            wedge_particles_metric = wedge_particles0;
        } else {
            std::vector<std::string> f2;
            tom::corr::parse_filelist(c_wedge_particles_metric, f2, wedge_particles_metric, 0);
            if (f2.size() != fparticles0.size()) {
                throw std::runtime_error("The wedge-metric list \"" + c_wedge_particles_metric + "\" does not contain the same number of particles.");
            }
            std::vector<std::string>::const_iterator it0 = fparticles0.begin();
            std::vector<std::string>::const_iterator it1 = f2.begin();
            std::size_t i=0;
            for (; it0!=fparticles0.end(); it0++,it1++,i++) {
                if (*it0 != *it1) {
                    clog << "# WARNING: metric file \"" << c_wedge_particles_metric << "\" contains different filename for particle #" << (i+1) << ": (\"" << *it1 << "\" but use \"" << *it0 << "\")." << std::endl;
                }
            }
        }

    c_binning = f("binning", 0);


    s = f("next_template", e_next_template2str(tom::av4::nt::AVG_1).c_str());
        try {
            c_next_template = tom::av4::nt::str2e_next_template(s);
        } catch (std::domain_error &e) {
            throw std::runtime_error("The parameter 'next_template' has an unrecognized value \"" + s + "\".");
        }

}


/****************************************************************************//**
 * \brief The constructor of manager_angular_refinement
 *
 * \param[in] f The config parameters to initialise the class.
 * \param[in,out] clog The output stream for the logging information of the
 *    class.
 *******************************************************************************/
tom::av4::manager_angular_refinement::manager_angular_refinement(GetPot &f, std::ostream &clog):
    manager(f, clog),
    c_angles_list(),
    c_wedge_particles_metric(),
    c_ignore_particle_wedge(false),
    c_next_template(tom::av4::nt::AVG_1),
    c_break_max_iteration_list(0),
    c_break_max_iteration_refinement(0),
    c_break_max_cc_list(2.),
    c_break_max_cc_refinement(2.),
    c_break_min_cc_change_list(-1.),
    c_break_min_cc_change_refinement(-1.),
    c_sampling_step_list(0.),
    c_sampling_step_refinement(2),
    c_binning(-1),
    angles_list(NULL),
    wedge_particles_metric(),
    cc_min_1(1.),
    cc_min_2(1.),
    ang_ref_status(tom::av4::manager_angular_refinement::GENERAL_SAMPLING),
    iter_rc_idx0(0) {

    std::string s;
    int i;


    c_angles_list = f("angles_list", c_angles_list.c_str());
        if (c_angles_list.empty()) {
            throw std::runtime_error("the refinement ANGULAR_REFINEMENT needs the parameter \"angles_list\" with the name of an em-file of angles (phi,psi,theta in rad) for the first iteration.");
        } else if (!helper::fs::is_regular(c_angles_list)) {
            throw std::runtime_error("the angle list \"" + c_angles_list + "\" must be the name of and em-file of angles (phi,psi,theta in rad), but the file does not exist.");
        }
        try {
            tom::Volume<double> *p;
            tom::read_from_em(p, c_angles_list, NULL,NULL,NULL, NULL,NULL);
            angles_list.reset(p);
        } catch (...) {
            throw std::runtime_error("the angle list \"" + c_angles_list + "\" could not be opened. It must be an em-file of angles (phi,psi,theta in rad).");
        }
        if (angles_list->getSizeX()!=3 || angles_list->getSizeZ()!=1) {
            throw std::runtime_error("the angle list \"" + c_angles_list + "\" must be an em-file of angles (phi,psi,theta in rad, 3xNx1).");
        }
        {
            double min, max;
            angles_list->minmax(min, max);
            if (min < -6.283185307179586477 || max > 6.283185307179586477) {
                clog << "# WARNING! The angles list \"" << c_angles_list << "\" contains values smaller then 0 or larger than 2*PI. There is nothing wrong with it, just be sure that these angles are in radians!" << std::endl;
            }
        }

    c_wedge_particles_metric = f("wedge_particles_metric", c_wedge_particles_metric.c_str());
        if (c_wedge_particles_metric.empty()) {
            wedge_particles_metric = wedge_particles0;
        } else {
            std::vector<std::string> f2;
            tom::corr::parse_filelist(c_wedge_particles_metric, f2, wedge_particles_metric, 0);
            if (f2.size() != fparticles0.size()) {
                throw std::runtime_error("The wedge-metric list \"" + c_wedge_particles_metric + "\" does not contain the same number of particles.");
            }
            std::vector<std::string>::const_iterator it0 = fparticles0.begin();
            std::vector<std::string>::const_iterator it1 = f2.begin();
            std::size_t i=0;
            for (; it0!=fparticles0.end(); it0++,it1++,i++) {
                if (*it0 != *it1) {
                    clog << "# WARNING: metric file \"" << c_wedge_particles_metric << "\" contains different filename for particle #" << (i+1) << ": (\"" << *it1 << "\" but use \"" << *it0 << "\")." << std::endl;
                }
            }
        }

    c_ignore_particle_wedge = f("ignore_particle_wedge", int(c_ignore_particle_wedge)) != 0;

    s = f("next_template", e_next_template2str(c_next_template).c_str());
        try {
            c_next_template = tom::av4::nt::str2e_next_template(s);
        } catch (std::domain_error &e) {
            throw std::runtime_error("The parameter 'next_template' has an unrecognized value \"" + s + "\".");
        }

    i = f("break_max_iteration_list", int(c_break_max_iteration_list));
        if (i < 0) {
            c_break_max_iteration_list = 0;
        } else {
            c_break_max_iteration_list = i;
        }

    i = f("break_max_iteration_refinement", int(c_break_max_iteration_refinement));
        if (i < 0) {
            c_break_max_iteration_refinement = 0;
        } else {
            c_break_max_iteration_refinement = i;
        }

    c_break_max_cc_list = f("break_max_cc_list", c_break_max_cc_list);

    c_break_max_cc_refinement = f("break_max_cc_refinement", c_break_max_cc_list);

    c_break_min_cc_change_list = f("break_min_cc_change_list", c_break_min_cc_change_list);
        if (c_break_min_cc_change_list < 0.) {
            c_break_min_cc_change_list = -1.;
        }

    c_break_min_cc_change_refinement = f("break_min_cc_change_refinement", c_break_min_cc_change_list);
        if (c_break_min_cc_change_refinement < 0.) {
            c_break_min_cc_change_refinement = -1.;
        }


    c_binning = (i=f("binning", int(c_binning))) < 0 ? 0 : i;

    c_sampling_step_list = f("sampling_step_list", -1.);
        if (c_sampling_step_list <= 0) {
            throw std::runtime_error("The parameter 'sampling_step_list' must be a positive floating point (angle in degree).");
        }

    i = f("sampling_step_refinement", 0);
        if (i < 0) {
            throw std::runtime_error("The parameter 'sampling_step_refinement' can not be negative.");
        } else if (i == 1) {
            throw std::runtime_error("The parameter 'sampling_step_refinement' can not be 1 (larger then 1 or zero for no angular refinement).");
        }
        c_sampling_step_refinement = i;

    current_binning = c_binning;

}









namespace {
/****************************************************************************//**
 * init the names of the next run, according to the value of next_template.
 *******************************************************************************/
void init_filenames_for_next_template(const tom::avg::filename_generator &f_gen, tom::av4::nt::e_next_template nt, std::string &ftemplate, std::string &fwedge, bool &has_wedge, std::size_t itemplate) {
    switch (nt) {
        case tom::av4::nt::AVG_N:
            f_gen.get_avg_n(itemplate).swap(ftemplate);
            f_gen.get_avgwedge_n(itemplate).swap(fwedge);
            has_wedge = true;
            break;
        case tom::av4::nt::AVG_1:
            f_gen.get_avg_1(itemplate).swap(ftemplate);
            f_gen.get_avgwedge_1(itemplate).swap(fwedge);
            has_wedge = true;
            break;
        case tom::av4::nt::AVG_SQRTN:
            f_gen.get_avg_sqrtn(itemplate).swap(ftemplate);
            f_gen.get_avgwedge_sqrtn(itemplate).swap(fwedge);
            has_wedge = true;
            break;
        case tom::av4::nt::AVG_S:
            f_gen.get_avg_s(itemplate).swap(ftemplate);
            f_gen.get_avgwedge_s(itemplate).swap(fwedge);
            has_wedge = true;
            break;
        case tom::av4::nt::AVG_NONE:
            f_gen.get_avg_1(itemplate).swap(ftemplate);
            fwedge = "";
            has_wedge = false;
            break;
        default:
            assert(0);
    };
}

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TFLOAT>
bool init_next_templates(   const tom::avg::filename_generator &f_gen,
                            const std::vector<tom::avg::st_average_result> &av,
                            std::ostream &clog,
                            tom::av4::nt::e_next_template nt,
                            std::vector<std::string> &ftemplates,
                            std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > &wedge_templates,
                            std::vector<std::size_t> &av_index,
                            std::size_t binning) {
    std::string f_t, f_w;
    ftemplates.resize(0);
    wedge_templates.resize(0);
    av_index.resize(0);
    bool has_wedge;
    const std::size_t avsize = av.size();
    for (std::size_t i=0; i<avsize; i++) {
        if (av[i].use_idx.empty()) {
            // No average computed.
            clog << "#        no average computed for template #" << i << " in the last run. REMOVE THIS CLASS." << std::endl;
        } else {
            ::init_filenames_for_next_template(f_gen, nt, f_t, f_w, has_wedge, i);
            if (!helper::fs::is_regular(f_t)) {
                clog << "#        average for template #" << i << " should be named \"" << f_t << "\", but the file does not exist." << std::endl;
                throw std::runtime_error("generateNextRun: file with average from previous run does not exist.");
            }
            ftemplates.push_back(f_t);

            if (has_wedge) {
                if (!helper::fs::is_regular(f_w)) {
                    clog << "#        wedge for average of template #" << i << " should be named \"" << f_w << "\", but the file does not exist." << std::endl;
                    throw std::runtime_error("generateNextRun: file with average from previous run does not exist.");
                }
                wedge_templates.push_back(boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> >(new tom::WedgeDescriptor_EMWedge<TFLOAT>(f_w, binning)));
            } else {
                wedge_templates.push_back(boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> >(new tom::WedgeDescriptor_NoWedge<TFLOAT>()));
            }

            av_index.push_back(i);
        }
    }
    return !ftemplates.empty();
}

/****************************************************************************//**
 *
 *******************************************************************************/
double get_minimal_cc_value(const std::vector<tom::avg::st_average_result> &av, tom::av4::nt::e_next_template next_template) {
    double cc_min = 1.;
    std::vector<tom::avg::st_average_result>::const_iterator avit;
    for (avit=av.begin(); avit!=av.end(); avit++) {
        switch (next_template) {
            case tom::av4::nt::AVG_NONE:
            case tom::av4::nt::AVG_1:
                assert((avit->ccval_1.get() && !avit->use_idx.empty()) || (!avit->ccval_1.get() && avit->use_idx.empty()));
                if (avit->ccval_1.get() && *avit->ccval_1<cc_min) { cc_min = *avit->ccval_1; }
                break;
            case tom::av4::nt::AVG_SQRTN:
                assert((avit->ccval_sqrtn.get() && !avit->use_idx.empty()) || (!avit->ccval_sqrtn.get() && avit->use_idx.empty()));
                if (avit->ccval_sqrtn.get() && *avit->ccval_sqrtn<cc_min) { cc_min = *avit->ccval_sqrtn; }
                break;
            case tom::av4::nt::AVG_N:
                assert((avit->ccval_n.get() && !avit->use_idx.empty()) || (!avit->ccval_n.get() && avit->use_idx.empty()));
                if (avit->ccval_n.get() && *avit->ccval_n<cc_min) { cc_min = *avit->ccval_n; }
                break;
            case tom::av4::nt::AVG_S:
                assert((avit->ccval_s.get() && !avit->use_idx.empty()) || (!avit->ccval_s.get() && avit->use_idx.empty()));
                if (avit->ccval_s.get() && *avit->ccval_s<cc_min) { cc_min = *avit->ccval_s; }
                break;
            default:
                assert(0);
        }
    }
    return cc_min;
}

/****************************************************************************//**
 *
 *******************************************************************************/
std::auto_ptr<tom::Volume<double> > angles_refinement(std::size_t refinement_steps, double angular_sampling_step, double phi, double psi, double theta) {

    assert(refinement_steps>=2 && angular_sampling_step>0.);

    // convert to radians...
    angular_sampling_step *= 0.01745329251994329576913914624236578987393;

    std::size_t j;
    signed long int j_phi, j_psi, j_theta;
    j = 2*refinement_steps-1;
    std::auto_ptr<tom::Volume<double> > res(new tom::Volume<double>(3, j*j*j, 1, NULL,NULL));

    assert(res->isContiguous());
    double *pres = &res->get();

    const signed long int rs_a = -static_cast<long signed int>(refinement_steps-1);
    const signed long int rs_z =  static_cast<long signed int>(refinement_steps  );

    for (j = 0, j_phi=rs_a; j_phi<rs_z; j_phi++) {
        for (j_psi=rs_a; j_psi<rs_z; j_psi++) {
            for (j_theta=rs_a; j_theta<rs_z; j_theta++, j++) {
                *pres++ = j_phi   * angular_sampling_step + phi;
                *pres++ = j_psi   * angular_sampling_step + psi;
                *pres++ = j_theta * angular_sampling_step + theta;
            }
        }
    }
    return res;
}

} // namespace

/****************************************************************************//**
 * creates all parameters for the next run (i.e. \c c.rc).
 *******************************************************************************/
bool tom::av4::manager_static::generateNextRun() {
    boost::shared_ptr<manager::run_config> rc_new_(new manager::run_config());
    manager::run_config &rc_new = *rc_new_;

    assert(wedge_particles0.size()==fparticles0.size() && wedge_templates0.size()==ftemplates0.size());

    std::size_t iit = rc.size();
    std::size_t ntemplates, nparticles;

    rc_new.binning                  = c_binning;
    rc_new.cc_mask_radius           = c0_cc_mask_radius;
    rc_new.reapply_mask_after_wedge = c0_reapply_mask_after_wedge;
    rc_new.ntemplates_amount        = 0;
    rc_new.nparticles_amount        = 0;
    rc_new.do_fsc                   = true;


    if (iit == 0) {
        rc_new.ftemplates           = ftemplates0;
        if (c_next_template == tom::av4::nt::AVG_NONE) {
            rc_new.wedge_templates.assign(rc_new.ftemplates.size(), boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> >(new tom::WedgeDescriptor_NoWedge<manager::TFLOAT>()));
        } else {
            rc_new.wedge_templates      = wedge_templates0;
        }
    } else {
        tom::av4::filename_generator f_gen(c_outputdir, c_prefix);
        tom::avg::filename_generator f_gen_avg(f_gen.get_outputdir(iit-1));

        const std::string filename = f_gen.get_logfilename_average(iit-1);
        std::vector<tom::avg::st_average_result> av;
        std::vector<std::size_t> av_index;

        tom::avg::parse_average(filename, av);
        if (rc[iit-1]->ftemplates.size() != av.size()) {
            throw std::runtime_error("ERROR parsing average log: number of parsed classes does not correspond to expected one.");
        }

        if (!::init_next_templates(f_gen_avg, av, clog, c_next_template, rc_new.ftemplates, rc_new.wedge_templates, av_index, rc_new.binning)) {
            clog << "# There are no classes left. BREAK." << std::endl;
            return false;
        }
    }
    rc_new.fparticles = fparticles0;
    rc_new.wedge_particles_avg = wedge_particles0;
    if (c_ignore_particle_wedge) {
        rc_new.wedge_particles_metric.assign(rc_new.fparticles.size(), boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> >(new tom::WedgeDescriptor_NoWedge<manager::TFLOAT>()));
    } else {
        rc_new.wedge_particles_metric = wedge_particles_metric;
    }

    ntemplates = rc_new.ftemplates.size();
    nparticles = rc_new.fparticles.size();

    rc_new.angles.assign(ntemplates*nparticles, boost::shared_ptr<tom::Volume<double> >(new tom::Volume<double>(*angles_list)));
    rc_new.particle_shift.assign(ntemplates*nparticles, std::pair<int, int>(0,0));

    rc.push_back(rc_new_);
    return true;
}



#define PI 3.141592653589793238512808959406186204433

namespace {
std::size_t compute_binning_faktor(std::size_t volsize, double increment) {
    # if 0
    // TODO: not yet implemented.
    if (increment <= 0) {
        return 0;
    }
    double threshold = 2. / tan(increment*(PI/180.));
    std::size_t i = 1;
    while (volsize/(i+1) > threshold) {
        i++;
    }
    #endif
    return 0;
}
}


/****************************************************************************//**
 * creates all parameters for the next run (i.e. \c c.rc).
 *******************************************************************************/
bool tom::av4::manager_angular_refinement::generateNextRun() {


    boost::shared_ptr<manager::run_config> rc_new_(new manager::run_config());
    manager::run_config &rc_new = *rc_new_;

    const std::size_t iit = rc.size();

    if (iit < 1) {

        // Initialise the first run-configuration.
        clog << "# initialise the overall sampling..." << std::endl;

        rc_new.binning = c_binning;

        rc_new.cc_mask_radius           = c0_cc_mask_radius;
        rc_new.reapply_mask_after_wedge = c0_reapply_mask_after_wedge;
        rc_new.ntemplates_amount        = 0;
        rc_new.nparticles_amount        = 0;
        rc_new.do_fsc                   = true;

        rc_new.ftemplates = ftemplates0;
        if (c_next_template == tom::av4::nt::AVG_NONE) {
            rc_new.wedge_templates.assign(rc_new.ftemplates.size(), boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> >(new tom::WedgeDescriptor_NoWedge<manager::TFLOAT>()));
        } else {
            rc_new.wedge_templates = wedge_templates0;
        }

        rc_new.fparticles = fparticles0;
        rc_new.wedge_particles_avg = wedge_particles0;
        if (c_ignore_particle_wedge) {
            rc_new.wedge_particles_metric.assign(rc_new.fparticles.size(), boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> >(new tom::WedgeDescriptor_NoWedge<manager::TFLOAT>()));
        } else {
            rc_new.wedge_particles_metric = wedge_particles_metric;
        }

        const std::size_t ntemplates = rc_new.ftemplates.size();
        const std::size_t nparticles = rc_new.fparticles.size();

        rc_new.angles.assign(ntemplates*nparticles, boost::shared_ptr<tom::Volume<double> >(new tom::Volume<double>(*angles_list)));
        rc_new.particle_shift.assign(ntemplates*nparticles, std::pair<int, int>(0,0));

        rc_new.angles.assign(ntemplates*nparticles, boost::shared_ptr<tom::Volume<double> >(new tom::Volume<double>(*angles_list)));
        rc_new.particle_shift.assign(ntemplates*nparticles, std::pair<int, int>(0,0));

        ang_ref_status = tom::av4::manager_angular_refinement::GENERAL_SAMPLING;
        iter_rc_idx0 = 0;
    } else {
        const manager::run_config &rc_old = *rc.at(iit-1);


        std::vector<tom::avg::st_average_result> av;
        std::vector<tom::avg::st_average_result>::const_iterator avit;
        tom::av4::filename_generator f_gen(c_outputdir, c_prefix);
        tom::avg::filename_generator f_gen_avg(f_gen.get_outputdir(iit-1));

        tom::avg::parse_average(f_gen.get_logfilename_average(iit-1), av);
        if (rc.at(iit-1)->ftemplates.size() != av.size()) {
            throw std::runtime_error("ERROR parsing average log: number of parsed classes does not correspond to expected one.");
        }

        std::vector<std::string> rc_new_ftemplates;
        std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > rc_new_wedge_templates;
        std::vector<std::size_t> av_index;
        if (!::init_next_templates(f_gen_avg, av, clog, c_next_template, rc_new_ftemplates, rc_new_wedge_templates, av_index, 0)) {
            clog << "# There are no classes left. BREAK." << std::endl;
            return false;
        }

        bool next_run_created = false;
        {
            // Check if the iteration with the last parameters has to be repeated...
            const std::size_t citeration = iit - iter_rc_idx0;
            assert(iter_rc_idx0<iit && iit>0);

            const double      &c_break_max_cc_        = ang_ref_status==tom::av4::manager_angular_refinement::GENERAL_SAMPLING ? c_break_max_cc_list        : c_break_max_cc_refinement       ;
            const std::size_t &c_break_max_iteration_ = ang_ref_status==tom::av4::manager_angular_refinement::GENERAL_SAMPLING ? c_break_max_iteration_list : c_break_max_iteration_refinement;
            const double      &c_break_min_cc_change_ = ang_ref_status==tom::av4::manager_angular_refinement::GENERAL_SAMPLING ? c_break_min_cc_change_list : c_break_min_cc_change_refinement;

            if (c_break_max_iteration_>0 && citeration>=c_break_max_iteration_) {
                clog << "# maximum number of iterations reached... proceed with angular refinement..." << std::endl;
            }  else {
                // Find the minimal correlation coefficient.
                cc_min_2 = cc_min_1;
                cc_min_1 = ::get_minimal_cc_value(av, c_next_template);
                const double cc_min_change_percentage = std::fabs((cc_min_1-cc_min_2)/cc_min_1*100.);

                if (cc_min_1 > c_break_max_cc_) {
                    clog << "# correlation of the computed average with their templates are all larger than " << c_break_max_cc_ << "... proceed with angular refinement..." << std::endl;
                } else if (citeration>=2 && c_break_min_cc_change_>0. && cc_min_change_percentage<c_break_min_cc_change_) {
                    clog << "# difference between the correlations of the computed average with their templates over the last two iterations was less then " << c_break_min_cc_change_ << "%... proceed with angular refinement..." << std::endl;
                } else {
                    {
                        std::streamsize p = clog.precision();
                        clog.precision(5);
                        clog << "# repeat iteration with updated template... (current iteration == " << citeration;
                        if (c_break_max_iteration_ > 0) {
                            clog << "/" << c_break_max_iteration_;
                        }
                        clog << ", CC=" << cc_min_1 << "/" << c_break_max_cc_;
                        if (iit>=2) {
                            clog << ", CCPREV=" << cc_min_2 << "=" << cc_min_change_percentage << "%/" << c_break_min_cc_change_;
                        }
                        clog << ")" << std::endl;
                        clog.precision(p);
                    }
                    // Copy the new parameters from the start of the current iteration.
                    const manager::run_config &rc_0 = *rc.at(iter_rc_idx0);
                    rc_new.binning                  = rc_0.binning;
                    rc_new.cc_mask_radius           = rc_0.cc_mask_radius;
                    rc_new.reapply_mask_after_wedge = rc_0.reapply_mask_after_wedge;
                    rc_new.ntemplates_amount        = rc_0.ntemplates_amount;
                    rc_new.nparticles_amount        = rc_0.nparticles_amount;
                    rc_new.do_fsc                   = rc_0.do_fsc;

                    rc_new.ftemplates.swap(rc_new_ftemplates);
                    rc_new.wedge_templates.swap(rc_new_wedge_templates);

                    rc_new.fparticles = rc_0.fparticles;
                    rc_new.wedge_particles_avg = rc_0.wedge_particles_avg;
                    rc_new.wedge_particles_metric = rc_0.wedge_particles_metric;

                    rc_new.angles = rc_0.angles;
                    rc_new.particle_shift = rc_0.particle_shift;

                    next_run_created = true;
                }
            }
        }
        if (!next_run_created) {
            assert(rc.at(iit-1)->ftemplates.size() == av.size());

            // current lokal iteration is finished. start a new one...

            if (ang_ref_status == tom::av4::manager_angular_refinement::GENERAL_SAMPLING) {
                // up to now, it did GENERAL_SAMPLING. Now switch to REFINEMENT...
                if (c_sampling_step_list <= 0) {
                    clog << "# The parameter is sampling_step_list set to 0. Thus do no angular refinement... BREAK." << std::endl;
                    return 0;
                }
                ang_ref_status = tom::av4::manager_angular_refinement::REFINEMENT;
                angular_sampling_step = c_sampling_step_list;
                clog << "# switch from general angular sampling to angular refinement..." << std::endl;
            }

            clog << "# Use angle refinement: maximum step = " << angular_sampling_step << "deg, with " << c_sampling_step_refinement << " steps." << std::endl;
            angular_sampling_step /= c_sampling_step_refinement;

            {
                // Is the angluar_sampling too low for the volume size?
                const double max_length_vector = volsize * 0.86602540378443864676372317075293618347; // The half length of the diagonal of the volume.
                const double max_length_vector_displacement = max_length_vector * sqrt(2*(1-cos(angular_sampling_step*0.01745329251994329577)));
                if (max_length_vector_displacement < 0.5) {
                    clog << "# Angular sampling is now at " << angular_sampling_step << "deg. For the volume of size " << volsize << " it makes no sense to refine further. BREAK" << std::endl;
                    return false;
                }
            }

            current_binning = current_binning / c_sampling_step_refinement;

            rc_new.cc_mask_radius = c0_cc_mask_radius;
            rc_new.reapply_mask_after_wedge = c0_reapply_mask_after_wedge;
            rc_new.ntemplates_amount = 0;
            rc_new.nparticles_amount = 0;
            rc_new.ftemplates = rc_new_ftemplates;
            rc_new.wedge_templates = rc_new_wedge_templates;
            rc_new.fparticles = fparticles0;
            rc_new.wedge_particles_avg = wedge_particles0;
            if (c_ignore_particle_wedge) {
                rc_new.wedge_particles_metric.assign(rc_new.fparticles.size(), boost::shared_ptr<tom::WedgeDescriptor<manager::TFLOAT> >(new tom::WedgeDescriptor_NoWedge<manager::TFLOAT>()));
            } else {
                rc_new.wedge_particles_metric = wedge_particles_metric;
            }

            const std::size_t ntemplates = rc_new.ftemplates.size();
            const std::size_t nparticles = rc_new.fparticles.size();
            rc_new.particle_shift.assign(ntemplates*nparticles, std::pair<int, int>(0,0));
            rc_new.do_fsc = true;

            { // Initialise the new angles list...
                const std::size_t ntemplates_old = rc_old.fparticles.size();
                const std::size_t nparticles_old = rc_old.fparticles.size();

                // Parse the peaklist to get the angles where to refine...
                std::vector<helper::triple<double, double, double> > vangles;
                std::vector<tom::cc_peak<tom::av4::manager::TFLOAT> > peak_list;
                tom::corr::parse_peaklist(f_gen.get_peakfilename(iit-1), peak_list, vangles, ntemplates_old, nparticles_old);

                assert(peak_list.size() == ntemplates_old*nparticles_old && vangles.size() == ntemplates_old*nparticles_old);
                assert(nparticles == nparticles_old);

                std::vector<std::size_t>::const_iterator vsit;
                std::size_t itemplate, itemplate_old, iparticle, iparticle_old;
                rc_new.angles.assign(ntemplates*nparticles, boost::shared_ptr<tom::Volume<double> >());
                for (itemplate=0; itemplate < av_index.size(); itemplate++) {
                    itemplate_old = av_index[itemplate];
                    assert(itemplate_old < rc_old.ftemplates.size());
                    for (iparticle=0; iparticle<nparticles; iparticle++) {
                        iparticle_old = iparticle; // There is not change in the particle list. thus... the index stays the same.
                        const tom::cc_peak<tom::av4::manager::TFLOAT> &peak_listi = peak_list[itemplate_old*nparticles_old + iparticle_old];
                        if (peak_listi.angle_idx >= 0) {
                            const helper::triple<double, double, double> &vanglesi = vangles[itemplate_old*nparticles_old + iparticle_old];
                            rc_new.angles[itemplate*nparticles + iparticle] = angles_refinement(c_sampling_step_refinement, angular_sampling_step, vanglesi.x, vanglesi.y, vanglesi.z);
                        }
                    }
                }
            }

            // Remember the current index as the start of a new iteration...
            iter_rc_idx0 = iit;
        }
    }


    rc.push_back(rc_new_);
    return true;
}















