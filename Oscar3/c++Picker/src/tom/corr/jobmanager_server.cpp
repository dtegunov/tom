/****************************************************************************//**
 * \file jobmanager_server.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    11.01.2007
 *******************************************************************************/
#include "tom/corr/jobmanager_server.hpp"


#include <set>
#include <map>
#include <assert.h>
#include <iostream>
#include <unistd.h>
#include <stdexcept>



#include <helper/GetPot>


#define TOM_MPI_MAIN__RANK_NONE           int(-1)




/****************************************************************************//**
 *
 *
 *******************************************************************************/
void tom::corr::ConfigDataServer::write_to_conf(const std::string &filename) const {

    std::ofstream f(filename.c_str(), std::ios_base::trunc | std::ios_base::out);
    if (!f) {
        std::ostringstream ss; ss << "write_to_conf: could not open \"" << filename << "\" for writeing."; throw std::runtime_error(ss.str());
    }

    f <<    "volsize=" << this->volsize <<                                      "\n"
            "binning=" << this->binning <<                                      "\n"
            "sphere_mask_inner_radius=" << this->sphere_mask_inner_radius <<    "\n"
            "sphere_mask_sigma=" << this->sphere_mask_sigma <<                  "\n"
            "sphere_mask_cutoff_radius=" << this->sphere_mask_cutoff_radius <<  "\n"
            "fftw_flag=" << tom::fftw::flag2str(this->fftw_flag) <<             "\n"
            "reapply_mask_after_wedge=" << this->reapply_mask_after_wedge <<    "\n"
            "fftw_wisdom_dir=" << this->fftw_wisdom_dir <<                      "\n"
            "nparticles_amount=" << this->nparticles_amount <<                  "\n"
            "ntemplates_amount=" << this->ntemplates_amount <<                  "\n"
            "peakfilename=" << this->peakfilename <<                            "\n"
            "cc_mask_radius=" << this->cc_mask_radius <<                        "\n"
            "outputdir=" << this->outputdir <<                                  "\n"
            "saveccvols=" << this->saveccvols <<                                "\n"
            "resume=" << this->resume <<                                        "\n"
            "force_files_exist=" << this->force_files_exist <<                  "\n"
            "logfile=" << this->logfile <<                                      "\n"
            "nice=" << this->nice <<                                            std::endl;
}







/****************************************************************************//**
 *
 *
 *******************************************************************************/
void tom::corr::ConfigDataServer::read_from_conf(const std::string &filename) {

    {
        if (!std::ifstream(filename.c_str()).good()) {
            std::ostringstream ss; ss << "Parsing error: Could not open \"" << filename << "\".";
            throw std::runtime_error(ss.str());
        }
    }
    GetPot f(filename.c_str(), 0x0, 0x0, " \t");


    this->volsize = f("volsize", 0);
    if (this->volsize <= 0) {
        std::ostringstream ss; ss << "Parsing error: The volsize can not be ommited (and must a positve integer).";
        throw std::runtime_error(ss.str());
    }
    this->binning = f("binning", 0);
    if (this->binning == 0) {
        this->binning = 1;
    }
    if (this->volsize / this->binning < 1) {
        std::ostringstream ss; ss << "Parsing error: The binning factor for the volume is to high (" << this->volsize << " / " << this->binning << ").";
        throw std::runtime_error(ss.str());
    }

    this->sphere_mask_inner_radius  = f("sphere_mask_inner_radius", ceil(this->volsize/2.) - 1.1);
    this->sphere_mask_sigma         = f("sphere_mask_sigma", 0);
    this->sphere_mask_cutoff_radius = f("sphere_mask_cutoff_radius", 0);

    this->fftw_flag = tom::fftw::str2flag(f("fftw_flag", "FFTW_MEASURE"));

    this->reapply_mask_after_wedge = f("reapply_mask_after_wedge", 0);
    this->fftw_wisdom_dir   = f("fftw_wisdom_dir", "");
    this->nparticles_amount = f("nparticles_amount", 0);
    this->ntemplates_amount = f("ntemplates_amount", 0);
    this->peakfilename      = f("peakfilename", "");
    this->cc_mask_radius    = f("cc_mask_radius", this->peakfilename.empty() ? 0. : (ceil(this->volsize/2.) - 1.1) * 0.5);
    this->outputdir         = f("outputdir", "");
    this->saveccvols        = f("saveccvols", 0);
    this->resume            = f("resume", 0);
    this->force_files_exist = f("force_files_exist", 1);
    this->logfile           = f("logfile", "");
    this->nice              = f("nice", 10);

    if (!this->peakfilename.empty() && this->cc_mask_radius <= 0) {
        std::ostringstream ss; ss << "Parsing error: The cc_mask_radius radius can not be <= 0 with a non empty peakfilename.";
        throw std::runtime_error(ss.str());
    }
    if (this->peakfilename.empty() && this->cc_mask_radius > 0) {
        std::ostringstream ss; ss << "Parsing error: Specifying a cc_mask_radius while having an empty peakfilename is not allowed.";
        throw std::runtime_error(ss.str());
    }
    if (this->ntemplates_amount < 0 || this->nparticles_amount < 0) {
        std::ostringstream ss; ss << "Parsing error: n" << (this->ntemplates_amount<0?"templates":"particles") << "_amount can not be less then 0.";
        throw std::runtime_error(ss.str());
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

}






/****************************************************************************//**
 *
 *
 *******************************************************************************/
void tom::corr::tag_angles(const std::vector<const tom::Volume<double> *> &input, std::vector<std::pair<int, boost::shared_ptr<tom::Volume<double> > > > &angles) {

    std::map<const tom::Volume<double> *, std::vector<std::size_t>, tom::volume_less<double> > keymap;
    std::map<const tom::Volume<double> *, std::vector<std::size_t>, tom::volume_less<double> >::const_iterator it;

    std::size_t i;
    std::size_t n = input.size();
    for (i=0; i<n; i++) {
        keymap[input[i]].push_back(i);
    }
    angles.resize(n);

    int tag = 0;
    for (it=keymap.begin(); it!=keymap.end(); it++, tag++) {
        const std::vector<std::size_t> &vidx = it->second;
        boost::shared_ptr<tom::Volume<double> > p(it->first ? new tom::Volume<double>(*it->first) : NULL);
        n = vidx.size();
        for (i=0; i<n; i++) {
            assert(vidx.at(i) < angles.size());
            const std::size_t &j = vidx[i];
            angles[j].first = tag;
            angles[j].second = p;
        }
    }

    #ifndef NDEBUG
    n = input.size();
    for (i=0; i<n; i++) {
        if (input[i]) {
            assert(angles[i].second.get() && *input[i]==*angles[i].second);
        } else {
            assert(!angles[i].second.get());
        }
    }
    #endif
}



/****************************************************************************//**
 * \brief Initialise the jobmanager.
 *
 * \param[in] ftemplates_ Filenames of the templates.
 * \param[in] fparticles_ Filenames of the particles.
 * \param[in] wedge_templates_
 * \param[in] wedge_particles_ The wedges for the particles. Must have the same
 *    size as \a fparticles_
 * \param[in] angles Vector with the angles to try for correlation. It must
 *    have the size of ftemplates_.size() * fparticles_.size(). The particles
 *    are one after the other row wise.
 *******************************************************************************/
template<typename TFLOAT>
void tom::corr::JobManagerServer<TFLOAT>::init_job(const std::vector<std::string> &ftemplates_,
                                            const std::vector<std::string> &fparticles_,
                                            const std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > &wedge_templates_,
                                            const std::vector<boost::shared_ptr<tom::WedgeDescriptor<TFLOAT> > > &wedge_particles_,
                                            const std::vector<const tom::Volume<double> *> &angles_) {
    if (ftemplates_.size()!=wedge_templates_.size() || fparticles_.size()!=wedge_particles_.size() || ftemplates_.size()*fparticles_.size()!=angles_.size()) {
        throw std::invalid_argument("The dimensions of the arguments missmatch");
    }

    this->ntemplates = ftemplates_.size();
    this->nparticles = fparticles_.size();
    this->ftemplates = ftemplates_;
    this->fparticles = fparticles_;
    this->wedge_templates = wedge_templates_;
    this->wedge_particles = wedge_particles_;
    tom::corr::tag_angles(angles_, this->angles);

    job_assignment.assign(this->ntemplates*this->nparticles, std::pair<char, int>(TOM_MPI_MAIN__JOB_NOT_ASSIGNED, TOM_MPI_MAIN__RANK_NONE));
    this->number_assigned_jobs = 0;
    this->number_finished_jobs = 0;

    std::size_t itemplate, iparticle, j=0;
    std::pair<int, boost::shared_ptr<tom::Volume<double> > > *pangles = &this->angles[0];
    for (itemplate=0; itemplate<this->ntemplates; itemplate++) {
        for (iparticle=0; iparticle<this->nparticles; iparticle++, j++, pangles++) {
            if (pangles->second.get() && (pangles->second->getSizeX()!=3 || pangles->second->getSizeZ()!=1)) {
                throw std::invalid_argument("The angles must be either NULL or a list of euler angles (Nx3x1).");
            }
            if (pangles->second.get()) {
                this->angles_not_yet_assigned[pangles->first].push_back(j);
            } else {
                this->job_assignment[j].first = TOM_MPI_MAIN__JOB_FINISHED;
                this->number_finished_jobs++;
            }
        }
    }
}


namespace {
/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TIT>
inline std::size_t choose_subgrid_count_elements(std::size_t rows, std::size_t cols, const std::vector<bool> &mat, TIT row_first, TIT row_last, TIT col_first, TIT col_last) {
    assert(mat.size() == rows*cols);

    //std::cout << "call choose_subgrid_count_elements: ";

    std::size_t count = 0;
    TIT it;
    for (; row_first!=row_last; row_first++) {
        for (it=col_first; it!=col_last; it++) {
            //std::cout << " " << *row_first << "," << *it << " ";
            assert(*row_first<rows && *it<cols);
            if (mat[*row_first * cols + *it]) {
                count++;
            }
        }
    }
    //std::cout << "   " << count << std::endl;
    return count;

}
}

namespace {
/****************************************************************************//**
 * This function takes a rectangular grid of dimension \a rows * \a cols with values
 * set to 1 at the positions \a idx and otherwise 0. (i.e. a sparse matrix)
 * The aim is to choose \a max_rows and \a max_cols so that the sum of the numbers
 * in the submatrix is maximal (or at least good :) ). \n
 * This is needed in assign_new_job with its list of the pending jobs (\a idx )
 * from the whole workload of \a max_rows * \a max_cols jobs. A optimal/good
 * subset of the jobs is needed.
 *******************************************************************************/
void choose_subgrid(std::size_t rows, std::size_t cols, const std::vector<std::size_t> &idx, std::size_t max_rows, std::size_t max_cols, std::vector<std::size_t> &result) {

    // The (naive) algorithm works as follows:
    // I rememenber the columns and rows of the "current", "best" result in the
    // sets result_rows and result_cols. The other are on the depot.
    // Then i try randomly to exchange columns and rows from the result and the depot
    // to improve the result.
    // Not really good. "may" be improved :)

    if (idx.empty()) {
        result.resize(0);
        return;
    }

    const std::size_t nidx = idx.size();

    std::set<std::size_t> result_rows, result_cols, depot_rows, depot_cols;
    std::set<std::size_t>::iterator it, jt;
    std::vector<bool> mat(rows*cols, false);

    std::size_t i, j, irow, icol;

    for (i=0; i<nidx; i++) {
        irow = idx[i] / cols;
        icol = idx[i] % cols;
        if (result_rows.find(irow)==result_rows.end() && depot_rows.find(irow)==depot_rows.end()) {
            if (result_rows.size() < max_rows) {
                result_rows.insert(irow);
            } else {
                depot_rows.insert(irow);
            }
        }
        if (result_cols.find(icol)==result_cols.end() && depot_cols.find(icol)==depot_cols.end()) {
            if (result_cols.size() < max_cols) {
                result_cols.insert(icol);
            } else {
                depot_cols.insert(icol);
            }
        }
        mat[irow*cols + icol] = true;
    }


    bool swapped;
    std::size_t result_linemin, result_linemin_val;
    std::size_t depot_linemax, depot_linemax_val;
    std::size_t cnt = 0;
    do {
        ++cnt;
        assert(cnt < 100);

        swapped = false;
        if ((it=depot_rows.begin()) != depot_rows.end()) {
            //std::cout << "search rows" << std::endl;
            jt = it; jt++;
            depot_linemax = *it;
            depot_linemax_val = choose_subgrid_count_elements(rows, cols, mat, it, jt, result_cols.begin(), result_cols.end());
            for (it++, jt++; it!=depot_rows.end(); it++, jt++) {
                j = choose_subgrid_count_elements(rows, cols, mat, it, jt, result_cols.begin(), result_cols.end());
                if (depot_linemax < j) {
                    depot_linemax_val = j;
                    depot_linemax = *it;
                }
            }
            //std::cout << "search rows result_linemin" << std::endl;
            it=result_rows.begin(); jt = it; jt++;
            result_linemin = *it;
            result_linemin_val = choose_subgrid_count_elements(rows, cols, mat, it, jt, result_cols.begin(), result_cols.end());
            for (it++, jt++; it!=result_rows.end(); it++, jt++) {
                j = choose_subgrid_count_elements(rows, cols, mat, it, jt, result_cols.begin(), result_cols.end());
                if (result_linemin_val > j) {
                    result_linemin_val = j;
                    result_linemin = *it;
                }
            }
            if (depot_linemax_val > result_linemin_val) {
                //std::cout << "swap rows " << result_linemin << "(" << result_linemin_val << ") - " << depot_linemax << "(" << depot_linemax_val << ") " << std::endl;
                depot_rows.erase(depot_linemax);
                result_rows.erase(result_linemin);
                depot_rows.insert(result_linemin);
                result_rows.insert(depot_linemax);
                swapped = true;
            }
        }
        if (!swapped && (it=depot_cols.begin()) != depot_cols.end()) {
            //std::cout << "search cols" << std::endl;
            jt = it; jt++;
            depot_linemax = *it;
            depot_linemax_val = choose_subgrid_count_elements(rows, cols, mat, result_rows.begin(), result_rows.end(), it, jt);
            for (it++, jt++; it!=depot_cols.end(); it++, jt++) {
                j = choose_subgrid_count_elements(rows, cols, mat, result_rows.begin(), result_rows.end(), it, jt);
                if (depot_linemax < j) {
                    depot_linemax_val = j;
                    depot_linemax = *it;
                }
            }
            //std::cout << "search cols result_linemin" << std::endl;
            it=result_cols.begin(); jt = it; jt++;
            result_linemin = *it;
            result_linemin_val = choose_subgrid_count_elements(rows, cols, mat, result_rows.begin(), result_rows.end(), it, jt);
            for (it++, jt++; it!=result_cols.end(); it++, jt++) {
                j = choose_subgrid_count_elements(rows, cols, mat, result_rows.begin(), result_rows.end(), it, jt);
                if (result_linemin_val > j) {
                    result_linemin_val = j;
                    result_linemin = *it;
                }
            }
            if (depot_linemax_val > result_linemin_val) {
                //std::cout << "swap cols " << result_linemin << "(" << result_linemin_val << ") - " << depot_linemax << "(" << depot_linemax_val << ") " << std::endl;
                depot_cols.erase(depot_linemax);
                result_cols.erase(result_linemin);
                depot_cols.insert(result_linemin);
                result_cols.insert(depot_linemax);
                swapped = true;
            }
        }
    } while (swapped);
    //std::cout << "rows: " << result_rows.size() << "vs" << depot_rows.size() << " - cols: " << result_cols.size() << "vs" << depot_cols.size() << std::endl;

    result.resize(0);
    for (it=result_rows.begin(); it!=result_rows.end(); it++) {
        for (jt=result_cols.begin(); jt!=result_cols.end(); jt++) {
            i = *it * cols + *jt;
            if (mat[i]) {
                result.push_back(i);
            }
        }
    }
}
}




/****************************************************************************//**
 * \brief Marks all jobs of a process as done.
 *
 * \param[in] rank The id of of the process.
 * \param[in] erase_value If true, the entrie is deleted by the function.
 *    Otherwise it is just cleared. This brings the list in an inconsistent state,
 *    because the list of currently assigned jobs can not be empty while still
 *    existing. This is usefull if later the list is filled again, so to keep
 *    the list without need to creating a new one.
 *
 * Marks the jobs of the process as finished.
 *******************************************************************************/
template<typename TFLOAT>
void tom::corr::JobManagerServer<TFLOAT>::mark_jobs_as_done(int rank, bool erase_value) {

    // Mark the job which were assigned before as done.
    std::map<int, std::vector<std::size_t> >::iterator iter = this->job_currently_assigned.find(rank);
    if (iter != this->job_currently_assigned.end()) {
        std::vector<std::size_t> &v = iter->second;
        for (std::size_t i=0; i<v.size(); i++) {
            //std::vector<std::pair<int, int> > job_assignment;
            assert(this->job_assignment.at(v.at(i)).first==TOM_MPI_MAIN__JOB_ASSIGNED && this->job_assignment.at(v.at(i)).second==rank &&
                   this->number_assigned_jobs >= 1 && this->number_finished_jobs<this->nparticles*this->ntemplates);
            this->job_assignment[v[i]].first = TOM_MPI_MAIN__JOB_FINISHED;
            this->number_assigned_jobs--;
            this->number_finished_jobs++;
        }
        if (erase_value) { this->job_currently_assigned.erase(rank); }
    }
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
bool tom::corr::JobManagerServer<TFLOAT>::get_assigned_jobs(int rank,
                        std::vector<std::string> &ftemplates,
                        std::vector<std::string> &fparticles,
                        std::vector<std::size_t> &templates_idx,
                        std::vector<std::size_t> &particles_idx,
                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_templates,
                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_particles,
                        const tom::Volume<double> *&angle,
                        std::vector<char> &process_pair) const {

    if (!this->get_assigned_jobs(rank, templates_idx, particles_idx, process_pair)) {
        // there are no jobs left. exit.
        ftemplates.resize(0);
        fparticles.resize(0);
        templates_idx.resize(0);
        particles_idx.resize(0);
        wedge_templates.resize(0);
        wedge_particles.resize(0);
        process_pair.resize(0);
        angle = NULL;
        return false;
    }
    std::size_t i, j, k, n, n2;

    n = templates_idx.size();
    n2 = particles_idx.size();

    ftemplates.resize(n);
    wedge_templates.resize(n);
    angle = NULL;
    for (i=0; i<n; i++) {
        j = templates_idx[i];
        ftemplates[i] = this->ftemplates[j];
        wedge_templates[i] = this->wedge_templates[j].get();
        if (!angle) {
            j = j*this->nparticles;
            for (k=0; k<n2; k++) {
                if (this->job_assignment[j+particles_idx[k]].first == TOM_MPI_MAIN__JOB_ASSIGNED &&
                    this->job_assignment[j+particles_idx[k]].second == rank) {
                    angle = this->angles[j+particles_idx[k]].second.get();
                    break;
                }
            }
        }
    }
    assert(angle);

    fparticles.resize(n2);
    wedge_particles.resize(n2);
    for (i=0; i<n2; i++) {
        j = particles_idx[i];
        fparticles[i] = this->fparticles[j];
        wedge_particles[i] = this->wedge_particles[j].get();
    }

    return true;

}





/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
bool tom::corr::JobManagerServer<TFLOAT>::assign_new_job( int rank,
                                                        std::size_t ntemplates_amount, std::size_t nparticles_amount,
                                                        std::vector<std::string> &ftemplates,
                                                        std::vector<std::string> &fparticles,
                                                        std::vector<std::size_t> &templates_idx,
                                                        std::vector<std::size_t> &particles_idx,
                                                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_templates,
                                                        std::vector<tom::WedgeDescriptor<TFLOAT> *> &wedge_particles,
                                                        const tom::Volume<double> *&angle,
                                                        std::vector<char> &process_pair) {

    assert(this->assert_status());
    assert(ntemplates_amount>=1 && nparticles_amount>=1);


    std::size_t i;
    int angle_id = 0;
    {
        i = 0;
        // Find the angle-list with the most common tuples to compute.
        for(std::map<int, std::vector<std::size_t> >::iterator iter=this->angles_not_yet_assigned.begin();
            iter!=this->angles_not_yet_assigned.end(); iter++ ) {
            if (iter->second.size() > i) {
                i = iter->second.size();
                angle_id =  iter->first;
            }
        }
    }

    if (i < 1) {
        // Mark the job which were assigned before as done.
        this->mark_jobs_as_done(rank);
    } else {
        // Mark the job which were assigned before as done.
        this->mark_jobs_as_done(rank, false);
        std::vector<std::size_t> &v_to_assign_now = this->job_currently_assigned[rank];
        std::vector<std::size_t> &v_angles_not_yet_assigned = this->angles_not_yet_assigned[angle_id];

        assert(v_angles_not_yet_assigned.size() == i);

        // Split a part of the list of not assigned angles in two halfes:
        choose_subgrid(this->ntemplates, this->nparticles, v_angles_not_yet_assigned, ntemplates_amount, nparticles_amount, v_to_assign_now);
        const std::size_t jobsize = v_to_assign_now.size();

        // Remove the indices from the new job (v_to_assign_now) from the non assigned list (v_angles_not_yet_assigned)
        if (jobsize == v_angles_not_yet_assigned.size()) {
            assert(std::set<std::size_t>(v_to_assign_now.begin(), v_to_assign_now.end()) == std::set<std::size_t>(v_angles_not_yet_assigned.begin(), v_angles_not_yet_assigned.end()));
            this->angles_not_yet_assigned.erase(angle_id);
        } else {
            std::set<std::size_t> s_to_assign_now(v_to_assign_now.begin(), v_to_assign_now.end());
            std::vector<std::size_t> vtmp;
            vtmp.reserve(v_angles_not_yet_assigned.size()-v_to_assign_now.size());
            for (std::vector<std::size_t>::iterator it=v_angles_not_yet_assigned.begin(); it!=v_angles_not_yet_assigned.end(); it++) {
                if (s_to_assign_now.find(*it) == s_to_assign_now.end()) {
                    vtmp.push_back(*it);
                }
            }
            v_angles_not_yet_assigned.swap(vtmp);
        }

        // Update the job_assignment list...
        for (i=0; i<jobsize; i++) {
            //std::vector<std::pair<int, int> > job_assignment;
            //std::cout << " assign job " << v_to_assign_now.at(i) << " == [" << (v_to_assign_now.at(i)/this->nparticles) << "," << (v_to_assign_now.at(i)%this->nparticles) << "]  " << this->job_assignment.at(v_to_assign_now.at(i)).first << "," << this->job_assignment.at(v_to_assign_now.at(i)).second << "]" << std::endl;
            assert(this->job_assignment.at(v_to_assign_now.at(i)).first==TOM_MPI_MAIN__JOB_NOT_ASSIGNED && this->job_assignment.at(v_to_assign_now.at(i)).second==TOM_MPI_MAIN__RANK_NONE);
            std::pair<int, int> &ja = this->job_assignment[v_to_assign_now[i]];
            ja.first = TOM_MPI_MAIN__JOB_ASSIGNED;
            ja.second = rank;
        }
        this->number_assigned_jobs += jobsize;
    }

    //std::cout << "print in assign_new_job..... " << std::endl;
    //this->print();

    assert(this->assert_status());
    //this->print();
    //std::cout << "ASSIGN_NEW_JOB: AFTER" << std::endl;


    return this->get_assigned_jobs(rank, ftemplates, fparticles, templates_idx, particles_idx, wedge_templates, wedge_particles, angle, process_pair);
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
void tom::corr::JobManagerServer<TFLOAT>::print_jobs(int rank) const {

    std::map<int, std::vector<std::size_t> >::const_iterator it = this->job_currently_assigned.find(rank);
    std::cout << "  process " << rank << " has ";
    if (it == this->job_currently_assigned.end()) {
        std::cout << "no jobs assigned.";
    } else {
        std::cout << it->second.size() << " jobs assigned: [";
        for (std::vector<std::size_t>::const_iterator it2=it->second.begin(); it2!=it->second.end(); it2++) {
            std::cout << " (";
            std::cout.width(3);
            std::cout << (*it2/this->nparticles) << ",";
            std::cout.width(3);
            std::cout << (*it2%this->nparticles) << ") ";
        }
        std::cout << "]";
    }
    std::cout << std::endl;
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
void tom::corr::JobManagerServer<TFLOAT>::print() const {

    assert(this->assert_status());

    std::cout << "Jobmanager: " << this->ntemplates << " templates and " << this->nparticles << " particles (" << (this->ntemplates*this->nparticles-this->number_assigned_jobs-this->number_finished_jobs) << " non assigned, " << this->number_assigned_jobs << " assigned and " << this->number_finished_jobs << " finished jobs)." << std::endl <<
                 "  Job assignement list: " << std::endl;
    std::size_t i=0, j=0;
    std::vector<std::pair<int, int> >::const_iterator it;
    for (it=this->job_assignment.begin(), j=0; it!=this->job_assignment.end(); it++, j++) {
        if (!i) {
            std::cout << "  ";
        }
        std::cout << "  (";
        if (it->first==TOM_MPI_MAIN__JOB_NOT_ASSIGNED) {
            std::cout << "-------";
        } else {
            std::cout << (it->first==TOM_MPI_MAIN__JOB_ASSIGNED?'A':'F') << ":";
            std::cout.width(5);
            std::cout << it->second;
        }
        std::cout.width(3);
        std::cout << this->angles[j].first <<
             " " << (void *)this->angles[j].second.get() <<
             ") ";
        if (++i >= this->nparticles) {
            std::cout << std::endl;
            i=0;
        }
    }
    for (std::map<int, std::vector<std::size_t> >::const_iterator it=this->job_currently_assigned.begin(); it!=this->job_currently_assigned.end(); it++) {
        this->print_jobs(it->first);
    }
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
bool tom::corr::JobManagerServer<TFLOAT>::get_assigned_jobs(int rank, std::vector<std::size_t> &templates_idx, std::vector<std::size_t> &particles_idx, std::vector<char> &process_pair) const {

    std::map<int, std::vector<std::size_t> >::const_iterator it_job_currently_assigned = this->job_currently_assigned.find(rank);

    if (it_job_currently_assigned == this->job_currently_assigned.end() || it_job_currently_assigned->second.empty()) {
        templates_idx.resize(0);
        particles_idx.resize(0);
        process_pair.resize(0);
        return false;
    }
    const std::vector<std::size_t> &jobs = it_job_currently_assigned->second;
    const std::size_t jobsize = jobs.size();

    assert(this->assert_status());
    #ifndef NDEBUG
    for (std::size_t ii=0; ii<jobs.size(); ii++) {
        assert(jobs.at(ii)<this->ntemplates*this->nparticles &&
               this->angles.at(jobs.at(0)).first == this->angles.at(jobs.at(ii)).first &&
               this->angles.at(jobs.at(0)).second.get() == this->angles.at(jobs.at(ii)).second.get());
    }
    #endif

    std::size_t itemplate, iparticle, i;
    std::set<std::size_t> set_templates;
    std::set<std::size_t> set_particles;
    for (i=0; i<jobsize; i++) {
        assert(jobs[i]>=0 && jobs[i]<this->ntemplates*this->nparticles &&
               this->job_assignment[(jobs[i] / this->nparticles)*this->nparticles + (jobs[i] % this->nparticles)].first==TOM_MPI_MAIN__JOB_ASSIGNED &&
               this->job_assignment[(jobs[i] / this->nparticles)*this->nparticles + (jobs[i] % this->nparticles)].second==rank);
        set_templates.insert(jobs[i] / this->nparticles);
        set_particles.insert(jobs[i] % this->nparticles);
    }
    const std::size_t ntemplates = set_templates.size();
    const std::size_t nparticles = set_particles.size();
    std::set<std::size_t>::iterator iter;


    templates_idx.resize(ntemplates);
    for (i=0, iter=set_templates.begin(); iter!=set_templates.end(); iter++, i++) {
        templates_idx[i] = *iter;
    }
    particles_idx.resize(nparticles);
    for (i=0, iter=set_particles.begin(); iter!=set_particles.end(); iter++, i++) {
        particles_idx[i] = *iter;
    }

    process_pair.resize(ntemplates*nparticles);
    for (i=0, itemplate=0; itemplate<ntemplates; itemplate++) {
        for ( iparticle=0; iparticle<nparticles; iparticle++, i++) {
            const std::pair<int, int> &ja = this->job_assignment[templates_idx[itemplate]*this->nparticles + particles_idx[iparticle]];
            process_pair[i] = ja.first==TOM_MPI_MAIN__JOB_ASSIGNED && ja.second==rank;
        }
    }

    return true;
}


/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
std::vector<std::pair<std::size_t, std::size_t> > tom::corr::JobManagerServer<TFLOAT>::get_jobs_in_state(int type) const {

    std::vector<std::pair<std::size_t, std::size_t> > res;
    res.reserve(this->ntemplates*this->nparticles);
    std::size_t itemplate, iparticle;
    const std::pair<int, int> *pjob_assignment = &this->job_assignment[0];
    for (itemplate=0; itemplate<this->ntemplates; itemplate++) {
        for (iparticle=0; iparticle<this->nparticles; iparticle++, pjob_assignment++) {
            if (pjob_assignment->first & type) {
                res.push_back(std::pair<std::size_t, std::size_t>(itemplate, iparticle));
            }
        }
    }
    return res;
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
std::size_t tom::corr::JobManagerServer<TFLOAT>::get_number_of_jobs_in_state(int type) const {

    std::size_t res = 0;
    std::size_t itemplate, iparticle;
    const std::pair<int, int> *pjob_assignment = &this->job_assignment[0];
    for (itemplate=0; itemplate<this->ntemplates; itemplate++) {
        for (iparticle=0; iparticle<this->nparticles; iparticle++, pjob_assignment++) {
            if (pjob_assignment->first & type) {
                res ++;
            }
        }
    }
    return res;
}


/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
std::pair<int,int> tom::corr::JobManagerServer<TFLOAT>::get_job_state(std::size_t itemplate, std::size_t iparticle) const {

    if (itemplate >= this->ntemplates || iparticle >= this->nparticles) {
        throw std::out_of_range("get_job_assignment tries to access ouf of bounds.");
    }
    return this->job_assignment[itemplate*this->nparticles + iparticle];
}





/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename TFLOAT>
int tom::corr::JobManagerServer<TFLOAT>::assert_status() const {


    int res = 1;
    res = res && this->ntemplates == this->ftemplates.size() &&
                 this->nparticles == this->fparticles.size() &&
                 this->ntemplates == this->wedge_templates.size() &&
                 this->nparticles == this->wedge_particles.size() &&
                 this->ntemplates*this->nparticles == this->angles.size() &&
                 this->ntemplates*this->nparticles == this->job_assignment.size();
    if (!res) {
        std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": basic members are inconsistend!" << std::endl;
        return 0;
    }

    std::size_t i;

    std::size_t n_jobs_not_assigned = 0;
    {
        // Check number_pending_jobs and job_currently_assigned and job_assignment
        std::size_t n2 = 0;
        std::map<int, std::size_t> njob_currently_assigned;
        for (i=0; i<this->ntemplates*this->nparticles; i++) {
            if (this->job_assignment[i].first==TOM_MPI_MAIN__JOB_NOT_ASSIGNED) {
                if (this->job_assignment[i].second!=TOM_MPI_MAIN__RANK_NONE) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (1)." << std::endl;
                    return 0;
                }
                n_jobs_not_assigned++;
                n2++;
            } else if (this->job_assignment[i].first==TOM_MPI_MAIN__JOB_ASSIGNED) {
                if (this->job_assignment[i].second==TOM_MPI_MAIN__RANK_NONE) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (2)." << std::endl;
                    return 0;
                }
                njob_currently_assigned[this->job_assignment[i].second] ++;
                n2++;
            } else if (this->job_assignment[i].first==TOM_MPI_MAIN__JOB_FINISHED) {
            } else {
                std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (3)." << std::endl;
                return 0;
            }
        }
        if (n2 != this->get_number_pending_jobs()) {
            std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": number_pending_jobs is inconsistend." << std::endl;
            return 0;
        }
        if (njob_currently_assigned.size() != this->job_currently_assigned.size()) {
            std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (1)." << std::endl;
            return 0;
        }
        std::map<int, std::vector<std::size_t> >::const_iterator it2;
        for (std::map<int, std::size_t>::const_iterator it=njob_currently_assigned.begin(); it!=njob_currently_assigned.end(); it++) {
            it2 = this->job_currently_assigned.find(it->first);
            if (it2==this->job_currently_assigned.end()) {
                std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (2)." << std::endl;
                return 0;
            }
            if (it2->second.size() != it->second) {
                std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (3)." << std::endl;
                return 0;
            }
        }

        for (it2=this->job_currently_assigned.begin(); it2!=this->job_currently_assigned.end(); it2++) {
            if (it2->second.empty()) {
                std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (4)." << std::endl;
                return 0;
            }
            for (std::vector<std::size_t>::const_iterator it=it2->second.begin(); it!=it2->second.end(); it++) {
                if (*it < 0 || *it>=this->ntemplates*this->nparticles) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (5)." << std::endl;
                    return 0;
                }
                if (this->job_assignment[*it].first!=TOM_MPI_MAIN__JOB_ASSIGNED || this->job_assignment[*it].second!=it2->first) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (6)." << std::endl;
                    return 0;
                }
            }
        }
    }
    {
        std::map<int, std::size_t> nangles_tag;
        for (i=0; i<this->ntemplates*this->nparticles; i++) {
            nangles_tag[this->angles.at(i).first]++;
            if (!this->angles.at(i).second.get() && (this->job_assignment.at(i).first!=TOM_MPI_MAIN__JOB_FINISHED || this->job_assignment.at(i).second!=TOM_MPI_MAIN__RANK_NONE)) {
                std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": angles is inconsistend (1)." << std::endl;
                return 0;
            }
        }


        // Check angles_not_yet_assigned
        std::size_t n = 0;
        for (std::map<int, std::vector<std::size_t> >::const_iterator it=this->angles_not_yet_assigned.begin(); it!=this->angles_not_yet_assigned.end(); it++) {
            n += it->second.size();
            for (std::vector<std::size_t>::const_iterator itv=it->second.begin(); itv!=it->second.end(); itv++) {
                if (*itv >= this->ntemplates*this->nparticles) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": angles_not_yet_assigned is inconsistend (1)." << std::endl;
                    return 0;
                }
                if (this->job_assignment.at(*itv).first!=TOM_MPI_MAIN__JOB_NOT_ASSIGNED) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": angles_not_yet_assigned is inconsistend (2)." << std::endl;
                    return 0;
                }
                if (this->angles.at(*itv).first != it->first) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": angles_not_yet_assigned is inconsistend (3)." << std::endl;
                    return 0;
                }
            }
        }
        if (n != n_jobs_not_assigned) {
            std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": angles_not_yet_assigned is inconsistend (4)." << std::endl;
            return 0;
        }
    }

    {
        for (std::map<int, std::vector<std::size_t> >::const_iterator it=this->job_currently_assigned.begin(); it!=this->job_currently_assigned.end(); it++) {
            const std::vector<std::size_t> &vv = it->second;
            for (std::size_t i=0; i<vv.size(); i++) {
                if (vv.at(i) > this->nparticles*this->ntemplates ||
                    this->job_assignment.at(vv.at(i)).first!=TOM_MPI_MAIN__JOB_ASSIGNED ||
                    this->job_assignment.at(vv.at(i)).second!=it->first ||
                    this->angles.at(vv.at(0)).first != this->angles.at(vv.at(i)).first ||
                    this->angles.at(vv.at(0)).second.get() != this->angles.at(vv.at(i)).second.get()) {
                    std::cerr << "assert_status: " << __FILE__ << ":" << __LINE__ << ": job_currently_assigned is inconsistend (1)." << std::endl;
                    return 0;
                }
            }
        }
    }

    return 1;
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
void tom::corr::ConfigDataServer::write_to_log(std::ostream &stream) const {

    stream <<   "# \n"
                "# ConfigDataServer: \n";
    if (!assert_status()) {
        stream << "# ConfigDataServer is in an inconsistent state!!!!\n";
    }

    stream <<   "volsize="                      << volsize                      << "\n"
                "binning="                      << binning                      << "\n"
                "sphere_mask_inner_radius="     << sphere_mask_inner_radius     << "\n"
                "sphere_mask_sigma="            << sphere_mask_sigma            << "\n"
                "sphere_mask_cutoff_radius="    << sphere_mask_cutoff_radius    << "\n"
                "cc_mask_radius="               << cc_mask_radius               << "\n"
                "fftw_flag="                    << fftw_flag                    << "\n"
                "reapply_mask_after_wedge="     << reapply_mask_after_wedge     << "\n"
                "fftw_wisdom_dir="              << fftw_wisdom_dir              << "\n"
                "ntemplates_amount="            << ntemplates_amount            << "\n"
                "nparticles_amount="            << nparticles_amount            << "\n"
                "peakfilename="                 << peakfilename                 << "\n"
                "outputdir="                    << outputdir                    << "\n"
                "saveccvols="                   << saveccvols                   << "\n"
                "resume="                       << resume                       << "\n"
                "force_files_exist="            << force_files_exist            << "\n"
                "logfile="                      << logfile                      << "\n"
                "nice="                         << nice                         << "\n"
                "#"                                                             << std::endl;

}





// TEMPLATE INSTANTIATIONS

template class tom::corr::JobManagerServer<float >;
template class tom::corr::JobManagerServer<double>;




