/****************************************************************************//**
 * \file corrpar.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    07.01.2008
 *******************************************************************************/


#include <iostream>
#include <fstream>

#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <fftw3.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <unistd.h>
#include <limits>
#include <sys/times.h>

#include <boost/shared_ptr.hpp>

#include "helper/auto_vector.hpp"
#include "helper/GetPot"
#include "helper/ostream_swapper.hpp"
#include "helper/snippets.hpp"

#include "tom/core/wedge.hpp"
#include "tom/core/volume_fcn.hpp"

#include "tom/corr/corr.hpp"
#include "tom/corr/correlation_handler.hpp"
#include "tom/corr/jobmanager_client.hpp"
#include "tom/corr/jobmanager_server.hpp"
#include "tom/corr/config_files.hpp"
#include "tom/mpi/mpi_fcn.hpp"



#define DBGMSG { std::cerr << "I AM HERE: " __FILE__ << ":" << __LINE__ << std::endl; }



#ifndef THREAD_SAFE
#  error Define THREAD_SAFE.
#endif







/****************************************************************************//**
 * \brief Loads the configuration.
 *
 * \param[in] argv The command line arguments from the main program.
 * \param[out] volsize The size of the originale volumes in the em-files.
 *    All volumes must have the same size and be cubic.
 * \param[out] ftemplates A list of all the templates (their em-filenames).
 * \param[out] fparticles A list of all the particles (their em-filenames).
 * \param[out] wedge_templates A vector of the same length as ftemplates containing
 *    a description of the wedge for the corresponding template. This is send the
 *    the reciever processes and they can use it to create the proper tom::Wedge object.
 * \param[out] wedge_particles Analog to \a wedge_templates the wedge for the particles.
 * \param[out] angles Returns a vector with a list of the ntemplates*nparticles
 *    angles to test for a specific combination of template x particle. Each
 *    Volume must be a 3xNx1 matrix with the angles phi, psi, theta in radians.
 * \param[out] binning The binning factor used when reading the volumes. Note that this
 *    is the binning factor as passed to tom_em_read and thus means that every \a binning
 *    voxel is taken (instaead of taking every 2^binning voxel).
 * \param[out] sphere_mask_inner_radius  For correlation a sphere mask is used. The 3 parameters
 *    \a sphere_mask_inner_radius \a sphere_mask_sigma \a sphere_mask_cutoff_radius correspond
 *    to the parameters of tom::init_spheremask which is used to initialise the mask.
 *    Take care that the 3 parameters are passed unchanged to tom::init_spheremask for a
 *    cubic mask of size (\a volsize / \a binning ).
 * \param[out] sphere_mask_sigma See \a sphere_mask_inner_radius.
 * \param[out] sphere_mask_cutoff_radius See \a sphere_mask_inner_radius.
 * \param[out] maskcc A boolean mask. If the peak is sought (i.e. \a peakfilename not empty)
 *    a correlation volume is only considered as eventual peak, if the corresponding
 *    voxel of the mask as a value different from 0. Note that mask is fft-shifted (i.e. the
 *    zero frequency element is in the center of the mask. A common value would be initialising
 *    the mask with a spheremask. maskcc must have the size \a volsize / \a binning.
 *    *maskcc == NULL is valid if the peak is not sought at all, or any! voxel will be considered.
 * \param[out] fftw_flags The fftw-flag used to create all fftw-plans (e.g. FFTW_MEASURE)
 * \param[out] fftw_wisdom_dir The directory where the fftw-wisdom is stored and where
 *    it will be saved after running the program. Every process passes it to ::load_fftw_wisdom()
 *    to load the wisdom.
 * \param[out] ntemplates_amount The number of templates grouped together as a work assignment
 *    for the processes. Every single work assignment consists of
 *    \a ntemplates_amount * \a nparticles_amount combinations of templates and particles
 *    and all their angles as saved in \a angles. A larger number means, that the local
 *    process needs more memory and more time to handle the job. But it will be also faster
 *    because certain operations can be saved. This parameter corresponds to the parameter
 *    ntemplates_amount of tom::corr3d_batch
 * \param[out] nparticles_amount Analog to \a ntemplates_amount.
 * \param[out] peakfilename The filename of the textfile where to save the peaks of all
 *    operations. Empty means that the final peak will not be computed and thus not saved.
 * \param[out] outputdir The directory where the single resulting correlation volumes
 *    will be saved. If empty, no correlation volumes will be saved or leaded upon
 *    resume.
 * \param[out] saveccvols (boolean value) If true, the computed maximum correlation volumes, its
 *    angles index volume and the list of angles will be saved in \a outputdir.
 *    The filename will be chosen as X_Y_ccvol.em, X_Y_ccidx.em X_Y_angle.em where
 *    X is the zero based index of the template in \a ftemplates and Y the zero based
 *    index in \a fparticles. X_Y_angle.em corresponds to the content in \a angles
 *    and the indexes saved in X_Y_ccidx.em are the index of the angle in that file.
 * \param[out] resume (boolean value) If true, first it will be looked in the file \a outputdir whether
 *    a part of the work is already done. It will be loaded from there instead of
 *    computing it new.
 * \param[out] force_files_exist Boolean parameter, wheter a non existing em-file (or invalid filesize)
 *    cases an exception of it is ignored tacitly.
 *******************************************************************************/
template<typename T, typename Tmask>
void parse_input_args(const std::vector<std::string> &argv, tom::corr::JobManagerServer<T> &jobmanager, tom::corr::ConfigDataServer &c) {


    if (argv.size() != 5) {
        std::clog <<    "# ERROR: " << (argv.empty()?"The program":argv[0]) << " needs 4 input arguments:\n"
                        "#   1: config_file.txt\n"
                        "#   2: filenames_templates.txt\n"
                        "#   3: filenames_particles.txt\n"
                        "#   4: angles.txt" << std::endl;
        throw std::invalid_argument("Wrong number of command line arguments (should be 4).");
    }

    try {

        std::vector<std::string> ftemplates;
        std::vector<std::string> fparticles;
        std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > > wedge_templates;
        std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > > wedge_particles;
        std::vector<boost::shared_ptr<tom::Volume<double> > > angles;

                std::clog <<    "# parse configuration (" << argv[1] << ")" << std::endl;
        c.read_from_conf(argv[1]);

                std::clog <<    "# parse filelist for templates (" << argv[2] << ")" << std::endl;
        tom::corr::parse_filelist(argv[2], ftemplates, wedge_templates, c.binning);
                std::clog <<    "#   " << ftemplates.size() << " templates read.\n"
                                "# parse filelist for particles (" << argv[3] << ")" << std::endl;
        tom::corr::parse_filelist(argv[3], fparticles, wedge_particles, c.binning);
                std::clog <<    "#   " << fparticles.size() << " particles read.\n"
                                "# parse the list of angles (" << argv[4] << ")" << std::endl;
        tom::corr::parse_angleslist(argv[4], ftemplates.size(), fparticles.size(), angles);



        std::vector<const tom::Volume<double> *> angles_(angles.size());
        std::size_t i = 0;
        std::vector<boost::shared_ptr<tom::Volume<double> > >::const_iterator it = angles.begin();
        for (; it!=angles.end(); it++, i++) {
            angles_[i] = it->get();
        }

                std::clog <<    "# initialise the job manager on server" << std::endl;
        jobmanager.init_job(ftemplates, fparticles, wedge_templates, wedge_particles, angles_);
                std::clog <<    "#   " << jobmanager.get_number_pending_jobs() << " jobs with " << jobmanager.get_number_remaining_angles_tags() << " different angles" << std::endl;


    } catch (std::exception &e) {
        std::clog << "# EXCEPTION in parse_input_args: " << e.what() << std::endl;
        throw;
    }
}




/****************************************************************************//**
 *
 *******************************************************************************/
MPI_Datatype get_mpi_type_client_config_data(const tom::corr::ConfigDataClient &c_const) {

    tom::corr::ConfigDataClient &c = const_cast<tom::corr::ConfigDataClient &>(c_const);

    MPI_Datatype type[8] = {    MPI_UNSIGNED_LONG,         // volsize
                                MPI_UNSIGNED,              // binning
                                MPI_DOUBLE,                // sphere_mask_inner_radius
                                MPI_DOUBLE,                // sphere_mask_sigma
                                MPI_DOUBLE,                // sphere_mask_cutoff_radius
                                MPI_DOUBLE,                // cc_mask_radius
                                MPI_UNSIGNED,              // fftw_flag
                                MPI_INT };                 // nice

    int blocklen[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
    int i;
    MPI_Aint base;
    MPI_Aint disp[8];
    MPI_Address(&c.volsize,                         &disp[0]);
    MPI_Address(&c.binning,                         &disp[1]);
    MPI_Address(&c.sphere_mask_inner_radius,        &disp[2]);
    MPI_Address(&c.sphere_mask_sigma,               &disp[3]);
    MPI_Address(&c.sphere_mask_cutoff_radius,       &disp[4]);
    MPI_Address(&c.cc_mask_radius,                  &disp[5]);
    MPI_Address(&c.fftw_flag,                       &disp[6]);
    MPI_Address(&c.nice,                            &disp[7]);

    MPI_Address(&c, &base);
    for (i=0; i<8; i++) { disp[i] -= base; }

    MPI_Datatype peak_type;
    MPI_Type_struct(8, blocklen, disp, type, &peak_type);
    MPI_Type_commit(&peak_type);
    return peak_type;
}


/****************************************************************************//**
 * \brief Broadcasts the configuration data of the client processes.
 *******************************************************************************/
int bcast_ConfigDataClient(int root, MPI_Comm comm, tom::corr::ConfigDataClient &c) {

    int res;

    MPI_Datatype c_type = get_mpi_type_client_config_data(c);
    if ((res=MPI_Bcast(&c, 1, c_type, root, comm)) != MPI_SUCCESS) {
        MPI_Type_free(&c_type);
        return res;
    }
    MPI_Type_free(&c_type);

    // fftw_wisdom_dir
    if ((res=tom::mpi::bcast_string(root, comm, c.fftw_wisdom_dir)) != MPI_SUCCESS) {
        return res;
    }

    // outputdir
    if ((res=tom::mpi::bcast_string(root, comm, c.outputdir)) != MPI_SUCCESS) {
        return res;
    }


    char bool2char[5] = { !c.return_peak_to_root, !c.saveccvols, !c.resume, !c.force_files_exist, !c.reapply_mask_after_wedge };

    if ((res=MPI_Bcast(bool2char, 5, MPI_CHAR, root, comm)) != MPI_SUCCESS) {
        return res;
    }

    c.return_peak_to_root = !bool2char[0];
    c.saveccvols = !bool2char[1];
    c.resume = !bool2char[2];
    c.force_files_exist = !bool2char[3];
    c.reapply_mask_after_wedge = !bool2char[4];

    assert(c.assert_status());

    return res;
}









#define CORR_MPI_TAG_SNDRCV_WORK int(5)


/****************************************************************************//**
 * \brief Analog to Send_work when no more work is on the server side.
 *
 *******************************************************************************/
void Send_work_finished(int dst, MPI_Comm comm) {
    const int tag = CORR_MPI_TAG_SNDRCV_WORK;
    unsigned long num__[2] = { 0, 0 };
    MPI_Send(num__, 2, MPI_UNSIGNED_LONG, dst, tag, comm);
}


/****************************************************************************//**
 * \brief Send the work from the root process to the worker process.
 *
 * \WARNING As currently implemented with static variables it is not thread safe!
 *******************************************************************************/
template<typename TFLOAT>
bool Send_work(int dst, MPI_Comm comm, tom::corr::JobManagerServer<TFLOAT> &jobmanager, int ntemplates_amount, int nparticles_amount) {

    const int tag = CORR_MPI_TAG_SNDRCV_WORK;


    #if THREAD_SAFE
    #  define __MAKE_STATIC__
    #else
    #  define __MAKE_STATIC__ static
    #endif
    // Make them static, so that everytime the same vector is used (and overwritten).
    // This way the memory has to be allocated only one (because a std::vector never
    // actually frees memory.
    __MAKE_STATIC__ std::vector<std::string> ftemplates;
    __MAKE_STATIC__ std::vector<std::string> fparticles;
    __MAKE_STATIC__ std::vector<tom::WedgeDescriptor<TFLOAT> *> wedge_templates;
    __MAKE_STATIC__ std::vector<tom::WedgeDescriptor<TFLOAT> *> wedge_particles;
    __MAKE_STATIC__ std::vector<std::size_t> templates_idx;
    __MAKE_STATIC__ std::vector<std::size_t> particles_idx;
    __MAKE_STATIC__ std::vector<char> process_pair;
    #undef __MAKE_STATIC__

    const tom::Volume<double> *angles;

    const bool res = jobmanager.assign_new_job(dst, ntemplates_amount, nparticles_amount, ftemplates, fparticles, templates_idx, particles_idx, wedge_templates, wedge_particles, angles, process_pair);

    unsigned long num__[2];


    #if 1
    {
        std::size_t nn = 0;
        for (long int i=ftemplates.size()*fparticles.size()-1; i>=0; i--) {
            nn += process_pair[i]!=0;
        }
        std::clog << "talk with P " << std::setw(4) << dst << ": ";
        if (res) {
            std::clog << "new job: " << std::setw(3) << ftemplates.size() << "x" << std::setw(4) << fparticles.size() << "x" << std::setw(5) << (angles?angles->getSizeY():0) << " = " << std::setw(5) << nn << "; ";
        } else {
            std::clog << "no more work left;               ";
        }
        std::clog <<    std::setw(5) << jobmanager.get_number_of_jobs_in_state(TOM_MPI_MAIN__JOB_NOT_ASSIGNED) << " n, " <<
                        std::setw(5) << jobmanager.get_number_of_jobs_in_state(TOM_MPI_MAIN__JOB_ASSIGNED) << " a, " <<
                        std::setw(5) << jobmanager.get_number_of_jobs_in_state(TOM_MPI_MAIN__JOB_FINISHED) << " f" <<
                        "  (" << helper::now2str() << ")" << std::endl;
    }
    #endif

    if (res) {

        //jobmanager.print();

        //std::cout << "========================================================" << std::endl;

        unsigned long &ntemplates = num__[0] = ftemplates.size();
        unsigned long &nparticles = num__[1] = fparticles.size();

        MPI_Send(num__, 2, MPI_UNSIGNED_LONG, dst, tag, comm);

        assert((ntemplates!=0&&nparticles!=0) && ntemplates==wedge_templates.size() && nparticles==wedge_particles.size() && ntemplates==templates_idx.size() && nparticles==particles_idx.size());

        tom::mpi::send_vstring(dst, tag, comm, ftemplates);
        tom::mpi::send_vstring(dst, tag, comm, fparticles);

        MPI_Send(&process_pair[0], ntemplates*nparticles, MPI_CHAR, dst, tag, comm);

        particles_idx.resize(ntemplates+nparticles);
        std::size_t i, j;
        for (i=0, j=nparticles; i<ntemplates; i++, j++) { particles_idx[j] = templates_idx[i]; }
        MPI_Send(&particles_idx[0], ntemplates+nparticles, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm);

        tom::mpi::send_wedges<TFLOAT>(dst, tag, comm, wedge_templates);
        tom::mpi::send_wedges<TFLOAT>(dst, tag, comm, wedge_particles);

        tom::mpi::send_volume<double>(dst, tag, comm, *angles);
    } else {
        Send_work_finished(dst, comm);
    }
    return res;
}






/****************************************************************************//**
 * \brief Recieve the work from the server process.
 *
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
bool Recv_work(int src, MPI_Comm comm, tom::corr::JobManagerClient<TFLOAT, TMASKCC> &jobmanager) {

    const int tag = CORR_MPI_TAG_SNDRCV_WORK;
    MPI_Status status;

    unsigned long num__[2];
    unsigned long &ntemplates = num__[0];
    unsigned long &nparticles = num__[1];
    MPI_Recv(num__, 2, MPI_UNSIGNED_LONG, src, tag, comm, &status);

    assert(ntemplates==0&&nparticles==0 || ntemplates!=0&&nparticles!=0);

    jobmanager.set_config_mode();

    std::vector<std::string> &ftemplates = jobmanager.access_ftemplates();
    std::vector<std::string> &fparticles = jobmanager.access_fparticles();
    std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > &wedge_templates = jobmanager.access_wedge_templates();
    std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > &wedge_particles = jobmanager.access_wedge_particles();
    std::vector<std::size_t> &templates_idx = jobmanager.access_templates_idx();
    std::vector<std::size_t> &particles_idx = jobmanager.access_particles_idx();
    std::auto_ptr<tom::Volume<double> > &angles = jobmanager.access_angles();
    std::vector<char> &process_pair = jobmanager.access_process_pair();

    const unsigned long n_prod = ntemplates*nparticles;
    ftemplates.resize(ntemplates);
    fparticles.resize(nparticles);
    wedge_templates.resize(n_prod);
    wedge_particles.resize(n_prod);
    angles.reset(NULL);
    templates_idx.resize(ntemplates);
    particles_idx.resize(nparticles);
    process_pair.resize(ntemplates*nparticles);


    if (ntemplates) {
        tom::mpi::recv_vstring(src, tag, comm, ftemplates);
        tom::mpi::recv_vstring(src, tag, comm, fparticles);

        MPI_Recv(&process_pair[0], ntemplates*nparticles, MPI_CHAR, src, tag, comm, &status);

        particles_idx.resize(ntemplates+nparticles);
        MPI_Recv(&particles_idx[0], ntemplates+nparticles, tom::mpi::get_MPI_Type<std::size_t>(), src, tag, comm, &status);
        std::size_t i, j;
        for (i=0, j=nparticles; i<ntemplates; i++, j++) { templates_idx[i] = particles_idx[j]; }
        particles_idx.resize(nparticles);

        tom::mpi::recv_wedges(src, tag, comm, wedge_templates);
        tom::mpi::recv_wedges(src, tag, comm, wedge_particles);

        tom::mpi::recv_volume<double>(src, tag, comm, angles);
    }

    #if 0
    {
        std::ostream &s = std::cerr;
        int my_rank;
        MPI_Comm_rank(comm, &my_rank);
        s << my_rank << " recieved new job: " << ntemplates << "x" << nparticles << "x" << (angles.get()?angles->getSizeY():0) << std::endl;
        for (int ttt=0; ttt<2; ttt++) {
            const std::vector<std::string> &f_ = (ttt? fparticles : ftemplates);
            const std::vector<boost::shared_ptr<tom::Wedge<TFLOAT> > > &w_ = (ttt ? wedge_particles : wedge_templates);
            const std::vector<std::size_t> &i_ = (ttt ? particles_idx : templates_idx);
            const std::string sname_ = (ttt ? "particle" : "template");
            const std::size_t &n_ = (ttt ? nparticles : ntemplates);
            for (size_t i=0; i<n_; i++) {
                s << " " << sname_ << " " << std::setw(3) << i << ": " << std::setw(3) << i_.at(i) << ": " << f_.at(i) << " - " << (*w_.at(i)).toString() << std::endl;
            }
        }
    }
    #endif


    jobmanager.set_working_mode();
    return ntemplates;
}






/****************************************************************************//**
 * Converts the (zero based) index \a idx which ranges from 0 to (\a s - 1)
 * to the corresponding index after doing an fft-shift and subtracting the center.
 * The index can be a non integer floating point. \n
 *******************************************************************************/
template<typename TFLOAT>
typename tom::cc_peak<TFLOAT>::idx_type ccpeak_to_shift(typename tom::cc_peak<TFLOAT>::idx_type idx, std::size_t s) {

    bool invert_sign = false;
    if (idx < 0) {
        invert_sign = true;
        idx = -idx;
    }

    typename tom::cc_peak<TFLOAT>::idx_type res = idx - tom::math::ceil(idx / static_cast<double>(s))*s;

    if (tom::math::abs(res) > idx) {
        res = idx;
    }

    return invert_sign ? -res : res;

}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
void server(MPI_Comm comm, const std::vector<std::string> &vargv) {

    const int msg_tag = 1;
    std::size_t j, itemplate, iparticle;

    int numprocs;
    int rank_server, rank_client;

    MPI_Comm_size(comm, &numprocs);
    MPI_Comm_rank(comm, &rank_server);

    tom::corr::JobManagerServer<TFLOAT> jobmanager;
    tom::corr::ConfigDataServer config_server;


    // auto_ptrs to the logfile stream and the ostream_swapper for clog.
    // flogss must be decleared after flog, so that hist deconstructor is called
    // first.
    std::auto_ptr<std::ostream> flog;
    std::auto_ptr<helper::ostream_swapper> flogss;


    std::vector<std::pair<pid_t, std::string> > v_hostnames;
    tom::mpi::gather_hostname(comm, v_hostnames, rank_server);
    assert(v_hostnames.size()==static_cast<unsigned>(numprocs));

    {
        // Temporary buffer the log in a stringstream before the logfilename is
        // clear.
        std::stringstream sclog;
        helper::ostream_swapper os(std::clog, sclog);
        try {
            {
                std::clog <<     "# " << helper::now2str() <<                                                                       "\n"
                                 "# " << (vargv.empty()?"corrpar":vargv[0]) << ":  PID=" << getpid() << ", rank=" << rank_server << "\n"
                                 "# Besides the server there are " << (numprocs-1) << " processes in the used communicator: "       "\n";
                for (j=0; j<static_cast<unsigned>(numprocs); j++) {
                    std::clog << "# " << std::setfill(' ') << std::setw(4) << j << ": " << std::setfill(' ') << std::setw(5) << v_hostnames[j].first << " on " << v_hostnames[j].second << '\n';
                }
                std::clog <<     "#"  << std::endl;
            }

            // Parse the input arguments from vargv...
            parse_input_args<TFLOAT, TMASKCC>(vargv, jobmanager, config_server);

            // Create the logfile-stream.
            if (!config_server.logfile.empty()) {
                flog.reset(new std::ofstream(config_server.logfile.c_str(), std::ios::out|std::ios::trunc));
                if (!flog->good()) {
                    std::clog << "# ERROR opening logfile " << config_server.logfile << " for writing: using STDERR" << std::endl;
                    flog.reset(NULL);
                }
            }
        } catch (std::exception &e) {
            // Print the error in STDERR and rethrow the exception.
            std::cerr << sclog.rdbuf() << "\nEXCEPTION: " << e.what() << std::endl;
            throw;
        } catch (...) {
            // Print the error in STDERR and rethrow the exception.
            std::cerr << sclog.rdbuf() << std::flush;
            throw;
        }

        os.swap_back();
        if (flog.get()) {
            flogss.reset(new helper::ostream_swapper(std::clog, *flog));
        }

        // Copy the content from the temporary buffer again into the clog.
        std::clog << sclog.rdbuf() << std::flush;
    }

    config_server.write_to_log(std::clog);


    const std::size_t ntemplates = jobmanager.get_ftemplates().size();
    const std::size_t nparticles = jobmanager.get_fparticles().size();

    tom::corr::ConfigDataClient config_client(config_server);

    // Broadcast parameters to clients...
    bcast_ConfigDataClient(rank_server, comm, config_client);

    tom::mpi::load_fftw_wisdom(config_client.fftw_wisdom_dir, comm);
    if (config_client.nice > 0) {
        // Make the server less nice than the others.
        nice(config_client.nice / 2);
    } else {
        nice(config_client.nice);
    }


    // The list of all the peaks, which is collected by the server (iff needed).
    std::vector<tom::cc_peak<TFLOAT> > peak_list;
    if (config_client.return_peak_to_root) {
        peak_list.assign(ntemplates*nparticles, tom::cc_peak<TFLOAT>());
    }


    // Temporary list of the peaks returned by the client.
    std::vector<tom::cc_peak<TFLOAT> > job_peak_list;

    // Temporary lists of the jobs which are assigned to the current
    // client (the server talks to).
    std::vector<std::size_t> templates_idx;
    std::vector<std::size_t> particles_idx;
    std::vector<char> process_pair;

    MPI_Status status;
    int client_has_results;

    double wtime_var;
    double wtime_sleep = 0;
    double wtime_all = MPI_Wtime();


    assert(jobmanager.assert_status());

    // Init ntemplates_amount and nparticles_amount.
    const double namount_factor = 1.5;
    double namount_ratio = 2. / 8.;      // desired ratio of ntemplates_amount to nparticles_amount.
    if (config_server.ntemplates_amount < 1 || config_server.nparticles_amount < 1) {
        if (config_server.ntemplates_amount >= 1) {
            config_server.nparticles_amount = std::max<unsigned>(1, static_cast<unsigned>(round(config_server.ntemplates_amount / namount_ratio)));
        } else if (config_server.nparticles_amount >= 1) {
            config_server.ntemplates_amount = std::max<unsigned>(1, static_cast<unsigned>(round(config_server.nparticles_amount * namount_ratio)));
        } else {
            const double f = sqrt(jobmanager.get_number_pending_jobs() / namount_factor / numprocs);
            config_server.ntemplates_amount = std::max<unsigned>(1, static_cast<unsigned>(ceil(f)));
            config_server.nparticles_amount = std::max<unsigned>(1, static_cast<unsigned>(ceil(config_server.ntemplates_amount/namount_ratio)));
        }
    } else {
        namount_ratio = static_cast<double>(config_server.ntemplates_amount) / config_server.nparticles_amount;
    }


    if (jobmanager.get_number_pending_jobs() > 0) {
        while (jobmanager.get_number_pending_jobs() > 0) { // Processing loop....

            while ((config_server.ntemplates_amount>1 || config_server.nparticles_amount>1) && config_server.ntemplates_amount*config_server.nparticles_amount*numprocs/2. > jobmanager.get_number_jobs_not_yet_assigned()) {
                if (config_server.ntemplates_amount == 1) {
                    config_server.nparticles_amount--;
                } else if (config_server.nparticles_amount == 1) {
                    config_server.ntemplates_amount--;
                } else if ((config_server.nparticles_amount-1)*namount_ratio > (config_server.ntemplates_amount-1)) {
                    config_server.nparticles_amount--;
                } else {
                    config_server.ntemplates_amount--;
                }
            }

            //std::cout << "HELLO WORLD: " << __FILE__ << ":" << __LINE__ << std::endl;
            // Block until message from client comes.
            wtime_var = MPI_Wtime();
            MPI_Probe(MPI_ANY_SOURCE, msg_tag, comm, &status);
            wtime_sleep += MPI_Wtime() - wtime_var;
            //std::cout << "HELLO WORLD: " << __FILE__ << ":" << __LINE__ << std::endl;

            // Give the requesting process work.
            rank_client = status.MPI_SOURCE;

            // Recieve the message from the client.
            MPI_Recv(&client_has_results, 1, MPI_INT, rank_client, msg_tag, comm, &status);

            //std::cout << "Server recieves request from " << rank_client << ": " << client_has_results << std::endl;
            //jobmanager.print();

            // Get the jobs currently assigned to the client.
            jobmanager.get_assigned_jobs(rank_client, templates_idx, particles_idx, process_pair);
            //std::cout << rank_client << " has a " << templates_idx.size() << "x" << particles_idx.size() << " job" << std::endl;

            assert(jobmanager.assert_status() && process_pair.size()==templates_idx.size()*particles_idx.size() && !(!config_client.return_peak_to_root && client_has_results) && (client_has_results&&!templates_idx.empty() || !client_has_results&&templates_idx.empty()));

            if (client_has_results) {

                // Recieve the peaklist from the client.
                tom::mpi::recv_peaklist<TFLOAT>(job_peak_list, rank_client, msg_tag, comm);

                assert(config_client.return_peak_to_root && job_peak_list.size()==process_pair.size());

                for (j=0, itemplate=0; itemplate<templates_idx.size(); itemplate++) {
                    for (iparticle=0; iparticle<particles_idx.size(); iparticle++, j++) {
                        //std::cout << "SERVER:::: " << (int)process_pair[j] << job_peak_list[j].angle_idx << " " << !(process_pair[j] && job_peak_list[j].angle_idx==-1) << (templates_idx.at(itemplate)<ntemplates) << (particles_idx.at(iparticle)<nparticles) << std::endl;
                        assert(templates_idx.at(itemplate)<ntemplates && particles_idx.at(iparticle)<nparticles);
                        if (process_pair[j]) {
                            peak_list[templates_idx[itemplate]*nparticles + particles_idx[iparticle]] = job_peak_list[j];
                        }
                    }
                }
            }

            // Send new job to the client (or no job - depending on what the jobmanager says.
            Send_work(rank_client, comm, jobmanager, config_server.ntemplates_amount, config_server.nparticles_amount);

            //jobmanager.print();
            //std::cout << "==============" << std::endl << std::endl;
        }
    } else {
        // There are not jobs at all. Send work finished to all clients.
        for (int i=0; i<numprocs; i++) {
            if (i != rank_server) {
                Send_work_finished(i, comm);
            }
        }
    }
    assert(jobmanager.assert_status());

    wtime_all = MPI_Wtime() - wtime_all;
    //std::cout << "server process ends (it slept for " << wtime_sleep << " of " << wtime_all << ")" << std::endl;

    if (config_client.return_peak_to_root) {
        assert(!config_server.peakfilename.empty());
        std::ofstream f(config_server.peakfilename.c_str());
        if (f) {
            f.precision(10);
            const std::vector<std::pair<int, boost::shared_ptr<tom::Volume<double> > > > &angles = jobmanager.get_angles();
            const signed long binning = config_server.binning>1 ? config_server.binning : 1;
            std::size_t volsize_effective = config_server.volsize / binning;
            const double peak_correction = 3;
            //const double peak_correction = binning>1 ? binning/2.-0.5 : 0.;

            f << "# num_template  num_particle  cc_val  shiftx  shifty  shiftz angle_idx angle_phi_rad angle_psi_rad angle_theta_rad" << std::endl;
            for (j=0, itemplate=0; itemplate<ntemplates; itemplate++) {
                for (iparticle=0; iparticle<nparticles; iparticle++, j++) {
                    const tom::cc_peak<TFLOAT> &peak = peak_list[j];
                    const tom::Volume<double> *angle = angles[j].second.get();
                    if (peak.angle_idx >= 0) {
                        assert(angle && static_cast<std::size_t>(peak.angle_idx) < angle->getSizeY());
                        f << itemplate << " " <<
                             iparticle << " " <<
                             (peak.val) << " " <<
                             (ccpeak_to_shift<TFLOAT>(peak.x, volsize_effective)*binning + peak_correction) << " " <<
                             (ccpeak_to_shift<TFLOAT>(peak.y, volsize_effective)*binning + peak_correction) << " " <<
                             (ccpeak_to_shift<TFLOAT>(peak.z, volsize_effective)*binning + peak_correction) << " " <<
                             peak.angle_idx << " " <<
                             angle->get(0, peak.angle_idx, 0) << " " <<
                             angle->get(1, peak.angle_idx, 0) << " " <<
                             angle->get(2, peak.angle_idx, 0) << " " <<
                             std::endl;
                    }
                }
            }
            std::clog << "# write peak-file \"" << config_server.peakfilename << "\"" << std::endl;
        } else {
            std::clog << "# ERROR writing the peak-file \"" << config_server.peakfilename << "\"" << std::endl;
        }
    }

    std::clog <<    "# finished: " << helper::now2str() << std::endl;
}




/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
void client(MPI_Comm comm, int root) {

    const int msg_tag = 1;
    int rank_client, numprocs;
    MPI_Comm_rank(comm, &rank_client);
    MPI_Comm_size(comm, &numprocs);

    tom::corr::ConfigDataClient config_client;

    {
        std::vector<std::pair<pid_t, std::string> > v;
        tom::mpi::gather_hostname(comm, v, root);
    }


    // Share the common data with all processes.
    bcast_ConfigDataClient(root, comm, config_client);

    tom::mpi::load_fftw_wisdom(config_client.fftw_wisdom_dir, comm);
    nice(config_client.nice);

    tom::corr::JobManagerClient<TFLOAT, TMASKCC> jobmanager(config_client);


    int client_has_results = 0;

    do { // Processing loop....

        //sleep(rank_client * 2);


        // say hallo to the server.
        client_has_results = jobmanager.client_has_results();

        MPI_Send(&client_has_results, 1, MPI_INT, root, msg_tag, comm);

        if (client_has_results) {
            // Send the result back to the server.
            //std::cout << rank_client << ": " << msg_cnt++ << ": send peak (" << jobmanager.get_peak_list().size() << ")" << std::endl;
            tom::mpi::send_peaklist<TFLOAT>(jobmanager.get_peak_list(), root, msg_tag, comm);
            client_has_results = false;
        }

        if ((Recv_work<TFLOAT>(root, comm, jobmanager))) {

            try {
                jobmanager.process();
            } catch(std::exception &e) {
                std::cerr << "EXCEPTION " << __FILE__ << ":" << __LINE__ << " at " << rank_client << ": " << e.what() << std::endl;
                throw;
            } catch (...) {
                std::cerr << "EXCEPTION " << __FILE__ << ":" << __LINE__ << " at " << rank_client << ": UNKNOWN." << std::endl;
                throw;
            }

        }
        //sleep((numprocs - rank_client) * 2);

    } while (jobmanager.client_has_work());

}






/****************************************************************************//**
 * \brief Main program.
 *
 *******************************************************************************/
int main(int argc, char *argv[]) {

    typedef double TFLOAT;
    typedef char TMASKCC;

    MPI_Init(&argc,&argv);

    const int root = 0;

    std::cout.precision(20);

    // Set the memory allocation functions.
    tom::fftw::setFcnMemFFTW<TFLOAT>();

    int numprocs_world, my_rank_world;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs_world);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_world);

    // Vector containing the command line arguments.
    std::vector<std::string> vargv;
    tom::mpi::startup(vargv, argc, argv);

    //std::cout << __FILE__ << ": " << __LINE__ << ": " << my_rank_world << std::endl;

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    {
        if (my_rank_world == root) {
            server<TFLOAT, TMASKCC>(comm, vargv);
        } else {
            client<TFLOAT, TMASKCC>(comm, root);
        }
    }
    tom::mpi::save_fftw_wisdom(comm);
    MPI_Comm_free(&comm);



    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}








