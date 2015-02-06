/****************************************************************************//**
 * \file jobmanager_client.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    13.01.2007
 *******************************************************************************/
#include "tom/corr/jobmanager_client.hpp"


#include <set>
#include <assert.h>
#include <iostream>


#include <tom/core/io.h>
#include <tom/core/volume_fcn.hpp>

#include <tom/corr/corr.hpp>
#include <tom/corr/filename_generators.hpp>

#include <helper/auto_vector.hpp>
#include <helper/snippets.hpp>




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
void tom::corr::JobManagerClient<TFLOAT, TMASKCC>::set_working_mode() {
    if (this->status==STATUS_CONFIG) {

        const std::size_t ntemplates = this->ftemplates.size();
        const std::size_t nparticles = this->fparticles.size();

        if ((ntemplates==0 && nparticles!=0) ||
            (ntemplates!=0 && nparticles==0) ||
            ntemplates != this->wedge_templates.size() ||
            nparticles != this->wedge_particles.size() ||
            ntemplates != this->templates_idx.size() ||
            nparticles != this->particles_idx.size() ||
            (ntemplates>0 && (!this->angles.get())) ||
            (this->angles.get() && (this->angles->getSizeX()!=3 || this->angles->getSizeZ()!=1)) ||
            this->process_pair.size() != ntemplates*nparticles) {
            throw std::runtime_error("set_working_mode: The job manager is in an inconsistante state (sizes). Can not switch to working mode.");
        }

        this->ntemplates = ntemplates;
        this->nparticles = nparticles;

        tom::corr::filename_generator f_gen(this->outputdir);

        if (this->saveccvols || this->resume) {
            std::size_t i, j, k;
            this->output_filenames.resize(ntemplates*nparticles);
            for (k=0, i=0; i<ntemplates; i++) {
                for (j=0; j<nparticles; j++, k++) {
                    typename tom::corr::JobManagerClient<TFLOAT, TMASKCC>::filename_triple &f = this->output_filenames[k];
                    f_gen.get_ccv(templates_idx[i], particles_idx[j]).swap(f.filename_ccvol);
                    f_gen.get_cci(templates_idx[i], particles_idx[j]).swap(f.filename_ccidx);
                    f_gen.get_ang(templates_idx[i], particles_idx[j]).swap(f.filename_angle);
                }
            }
        }

        if (this->get_peaks) {
            this->peak_list.assign(ntemplates*nparticles, tom::cc_peak<TFLOAT>());
        }

        this->status = STATUS_PROCESS;
    }

    assert(this->assert_status());
}

#include <sys/types.h>
#include <unistd.h>


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
tom::corr::JobManagerClient<TFLOAT, TMASKCC>::JobManagerClient(const tom::corr::ConfigDataClient &c)
    : outputdir(c.outputdir),
      resume(c.resume),
      saveccvols(c.saveccvols),
      get_peaks(c.return_peak_to_root),
      status(STATUS_CONFIG),
      volsize(c.volsize),
      binning(c.binning),
      fftw_flag(c.fftw_flag),
      reapply_mask_after_wedge(c.reapply_mask_after_wedge),
      force_files_exist(c.force_files_exist) {

    assert(c.assert_status());

    if (!this->binning) {
        this->binning = 1;
    }

    this->volsize_effective = c.get_volsize_effective();

    if (this->volsize_effective < 1) {
        throw std::invalid_argument("The effective volumesize after binning can not be empty.");
    }

    if ((this->resume||this->saveccvols) && this->outputdir.empty()) {
        throw std::invalid_argument("To resume it is not allowed to pass an empty output-directory");
    }

    if (!this->saveccvols && !this->get_peaks) {
        throw std::invalid_argument("You must set either saveccvols and/or get_peaks to true, otherwise no work is done.");
    }

    if (this->get_peaks && c.cc_mask_radius > 0) {
        tom::Volume<TMASKCC> *pmaskcc;
        this->maskcc.reset(pmaskcc = new tom::Volume<TMASKCC>(this->volsize_effective, this->volsize_effective, this->volsize_effective, NULL,NULL));
        tom::Volume<TMASKCC> v(this->volsize_effective, this->volsize_effective, this->volsize_effective, NULL,NULL);
        tom::init_spheremask(v, c.cc_mask_radius/this->binning, 0,0);
        tom::fftshift(v, *pmaskcc, true);
    }

    if (c.sphere_mask_inner_radius>0 || c.sphere_mask_sigma>0) {
        tom::Volume<TFLOAT> *p;
        this->mask_sphere.reset(p = new tom::Volume<TFLOAT>(this->volsize_effective, this->volsize_effective, this->volsize_effective, NULL,NULL));
        tom::init_spheremask(*p, c.sphere_mask_inner_radius/this->binning, c.sphere_mask_sigma/this->binning, c.sphere_mask_cutoff_radius/this->binning);
    }
}

/****************************************************************************//**
 * This method does some simple checks of consistency. No warranty for compleatness :)
 * Should be called by \c assert
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
int tom::corr::JobManagerClient<TFLOAT, TMASKCC>::assert_status() const {


    if (this->status==STATUS_CONFIG) {
        return true;
    }

    int res = !((this->resume||this->saveccvols) && outputdir.empty()) &&
           ((this->ntemplates==0&&this->nparticles==0) || (this->ntemplates!=0&&this->nparticles!=0)) &&
           this->ftemplates.size() == this->ntemplates &&
           this->fparticles.size() == this->nparticles &&
           this->wedge_templates.size() == this->ntemplates &&
           this->wedge_particles.size() == this->nparticles &&
           this->templates_idx.size() == this->ntemplates &&
           this->particles_idx.size() == this->nparticles &&
           (this->ntemplates==0 || this->ntemplates!=0&&this->angles.get()) &&
           (!this->angles.get() || (this->angles->getSizeX()==3 || this->angles->getSizeZ()==1)) &&
           this->process_pair.size() == this->ntemplates*this->nparticles;

    if (res) {
        res = ( (this->resume||this->saveccvols) && this->output_filenames.size()==this->nparticles*this->ntemplates) ||
              (!(this->resume||this->saveccvols) && this->output_filenames.empty());
    }

    if (res) {
        res = ( this->get_peaks && this->peak_list.size()==this->nparticles*this->ntemplates) ||
              (!this->get_peaks && this->peak_list.empty());
    }

    if (res) {
        std::set<std::size_t> s(this->templates_idx.begin(), this->templates_idx.end());
        res = s.size() == this->ntemplates;
    }
    if (res) {
        std::set<std::size_t> s(this->particles_idx.begin(), this->particles_idx.end());
        res = s.size() == this->nparticles;
    }

    if (res) {
        res = (this->get_peaks || this->saveccvols) && this->binning>0 && this->volsize>0 && this->volsize/this->binning==this->volsize_effective &&
              !(!this->get_peaks && this->maskcc.get()) && (!this->maskcc.get() || this->maskcc.get()&&this->maskcc->is_equal_size(this->volsize_effective, this->volsize_effective, this->volsize_effective));
    }

    if (res) {
        res = !this->mask_sphere.get() || (this->mask_sphere.get()&&this->mask_sphere->is_equal_size(this->volsize_effective, this->volsize_effective, this->volsize_effective));
    }


    return res;

}





/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
const std::vector<tom::cc_peak<TFLOAT> > &tom::corr::JobManagerClient<TFLOAT, TMASKCC>::get_peak_list() const {
    if (!this->get_peaks) {
        throw std::runtime_error("The jobmanager does not compute the peaks.");
    }
    if (this->status != STATUS_FINISHED) {
        throw std::runtime_error("The jobmanager has not yet processed its work.");
    }
    return this->peak_list;
}


/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
void tom::corr::JobManagerClient<TFLOAT, TMASKCC>::init_effective_jobs() {

    assert(this->assert_status() && this->status == STATUS_PROCESS);

    std::size_t i, j, itemplate, iparticle;

    if (this->resume) {
        // initialise resume_found (and peak_list, iff needed) by looking if there is a file on disk.

        tom::cc_peak<TFLOAT> ipeak;
        tom::cc_peak<TFLOAT> *p_ipeak = (this->get_peaks ? &ipeak : NULL);

        std::vector<char> new_process_pair(ntemplates*nparticles, false);
        std::vector<std::size_t> template_found(ntemplates, 0);
        std::vector<std::size_t> particle_found(nparticles, 0);

        // See for each pair iff it is already at file (and extract the peak, iff needed).
        for (j=0, itemplate=0; itemplate<this->ntemplates; itemplate++) {
            for (iparticle=0; iparticle<this->nparticles; iparticle++, j++) {
                //std::cout << "       check: " << itemplate << "," << iparticle << ": " << this->process_pair[j];
                if (this->process_pair[j]) {
                    if (this->check_existing(itemplate, iparticle, p_ipeak)) {
                        //std::cout << " found";
                        template_found[itemplate]++;
                        particle_found[iparticle]++;
                        if (this->get_peaks) {
                            peak_list[j] = ipeak;
                        }
                    } else {
                        new_process_pair[j] = true;
                    }
                } else {
                    template_found[itemplate]++;
                    particle_found[iparticle]++;
                }
                //std::cout << std::endl;
            }
        }

        this->job_ntemplates = 0;
        this->job_ftemplates.clear();
        this->job_wedge_templates.clear();
        this->job_templates_idx.clear();
        for (itemplate=0; itemplate<ntemplates; itemplate++) {
            if (template_found[itemplate] < nparticles) {
                this->job_templates_idx.push_back(itemplate);
                this->job_ftemplates.push_back(this->ftemplates[itemplate]);
                this->job_wedge_templates.push_back(this->wedge_templates[itemplate].get());
                this->job_ntemplates++;
            }
        }

        this->job_nparticles = 0;
        this->job_fparticles.clear();
        this->job_wedge_particles.clear();
        this->job_particles_idx.clear();
        for (iparticle=0; iparticle<nparticles; iparticle++) {
            if (particle_found[iparticle] < ntemplates) {
                this->job_particles_idx.push_back(iparticle);
                this->job_fparticles.push_back(this->fparticles[iparticle]);
                this->job_wedge_particles.push_back(this->wedge_particles[iparticle].get());
                this->job_nparticles++;
            }
        }

        this->job_process_pair.resize(this->job_ntemplates*this->job_nparticles);
        for (j=0, itemplate=0; itemplate<this->job_ntemplates; itemplate++) {
            for (iparticle=0; iparticle<this->job_nparticles; iparticle++, j++) {
                this->job_process_pair[j] = new_process_pair[this->job_templates_idx[itemplate]*this->nparticles + this->job_particles_idx[iparticle]];
            }
        }

    } else {
        // Redo everything: copy jobs.
        this->job_ntemplates = this->ntemplates;
        this->job_nparticles = this->nparticles;
        this->job_ftemplates = this->ftemplates;
        this->job_fparticles = this->fparticles;
        this->job_wedge_templates.resize(this->ntemplates); for (i=0; i<this->ntemplates; i++) { this->job_wedge_templates[i] = this->wedge_templates[i].get(); }
        this->job_wedge_particles.resize(this->nparticles); for (i=0; i<this->nparticles; i++) { this->job_wedge_particles[i] = this->wedge_particles[i].get(); }
        this->job_templates_idx.resize(this->ntemplates); for (i=0; i<this->ntemplates; i++) { this->job_templates_idx[i] = i; }
        this->job_particles_idx.resize(this->nparticles); for (i=0; i<this->nparticles; i++) { this->job_particles_idx[i] = i; }
        this->job_process_pair = this->process_pair;
    }
}




/****************************************************************************//**
 * \brief Looks on the filesystem if there are already the files of the result.
 *
 * \param[in] angles The list of angles. The content of \a filename_angle must
 *    equal this parameter.
 * \param[in] volsize The size of the correlation volume saved to disk.
 * \param[in] filename_ccvol Filename of the correlation volume.
 * \param[in] filename_ccidx Filename of the index volume
 * \param[in] filename_angle Filename of the angles list.
 * \param[out] ipeak If NULL, only the existence and the dimension of the
 *    volume files are checkt. If non NULL, the volumes are searched and
 *    the peak found is returned.
 * \param[in] maskcc Only relevant if \a ipeak not equal NULL. Then this
 *    is the mask (of size volsize) which is used for finding the peak
 *    (see tom::peak).
 * \returns True, if the files where found, the angles file contains the same
 *    thing as \a angles and the dimensions and datatypes of the other files
 *    match. Otherwise false.
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
bool tom::corr::JobManagerClient<TFLOAT, TMASKCC>::check_existing(std::size_t itemplate, std::size_t iparticle,
                    tom::cc_peak<TFLOAT> *ipeak) const {

    assert( ((!ipeak && !this->get_peaks) || (ipeak && this->get_peaks)) && this->assert_status() && this->status == STATUS_PROCESS && itemplate<this->ntemplates && iparticle<this->nparticles && this->process_pair.at(itemplate*this->nparticles+iparticle));

    const typename tom::corr::JobManagerClient<TFLOAT, TMASKCC>::filename_triple &f = this->output_filenames[itemplate*this->nparticles + iparticle];

    {
        tom::Volume<double> *pangles;
        std::auto_ptr<tom::Volume<double> > angles2;
        // Check existance and content of angles.
        try {
            tom::read_from_em<double>(pangles, f.filename_angle, NULL,NULL,NULL, NULL, NULL,NULL);
            angles2.reset(pangles);
        } catch (...) {
            return false;
        }
        if (!(*pangles == *this->angles)) {
            return false;
        }
    }

    tom_io_em_header header;

    //Check the content of the index-file for the angles.
    std::auto_ptr<tom::Volume<int32_t> > vccidx;
    {
        try {
            tom::Volume<int32_t> *p;
            tom::read_from_em<int32_t>(p, f.filename_ccidx, NULL,NULL,NULL, &header, NULL,NULL);
            vccidx.reset(p);
        } catch (...) {
            return false;
        }
        if (!vccidx->is_equal_size(this->volsize_effective, this->volsize_effective, this->volsize_effective) || tom_io_em_get_iotype(&header) != tom::get_tom_io_type<int32_t>()) {
            return false;
        }
        typename tom::Volume<int32_t>::element_type vmin, vmax;
        vccidx->minmax(vmin, vmax);
        if (vmin<0 || static_cast<uint32_t>(vmax)>=this->angles->getSizeY()) {
            return false;
        }
    }


    if (ipeak) {

        // Find the peak...

        std::auto_ptr<tom::Volume<TFLOAT> > vccvol;
        {
            try {
                tom::Volume<TFLOAT> *p;
                tom::read_from_em<TFLOAT>(p, f.filename_ccvol, NULL,NULL,NULL, &header, NULL,NULL);
                vccvol.reset(p);
            } catch (...) {
                return false;
            }
            if (!vccvol->is_equal_size(this->volsize_effective, this->volsize_effective, this->volsize_effective) || tom_io_em_get_iotype(&header) != tom::get_tom_io_type<TFLOAT>()) {
                return false;
            }
        }

        std::vector<tom::st_idx> peaks;
        if (this->maskcc.get()) {
            peaks = tom::peak(*vccvol, *this->maskcc);
        } else {
            peaks = tom::peak(*vccvol);
        }

        if (peaks.size() >= 1) {
            const std::size_t &px = peaks[0].x;
            const std::size_t &py = peaks[0].y;
            const std::size_t &pz = peaks[0].z;
            ipeak->x = px;
            ipeak->y = py;
            ipeak->z = pz;
            ipeak->val = vccvol->get(px, py, pz);
            ipeak->angle_idx = vccidx->get(px, py, pz);
        } else {
            *ipeak = tom::cc_peak<TFLOAT>();
        }
    } else {
        // Only check the size of the volume.
        if (tom_io_em_read_header(f.filename_ccvol.c_str(), "rb", &header, NULL) != TOM_ERR_OK ||
            tom_io_em_get_iotype(&header) != tom::get_tom_io_type<TFLOAT>() ||
            header.dims[0] != this->volsize_effective ||
            header.dims[1] != this->volsize_effective ||
            header.dims[2] != this->volsize_effective) {
            return false;
        }
    }
    return true;
}










/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
void tom::corr::JobManagerClient<TFLOAT, TMASKCC>::process() {

    if (this->status == STATUS_CONFIG) {
        throw std::runtime_error("JobManagerClient: can not process() in configuration mode");
    } else if (this->status == STATUS_FINISHED) {
        throw std::runtime_error("JobManagerClient: process already done.");
    }

    #if 0
    std::cout << "process job: " << this->ntemplates << "*" << this->nparticles << "*" << this->angles->getSizeY() << std::endl;
    for (std::size_t i=0; i<this->ntemplates; i++) {
        std::cout << "  " << i << ": " << this->templates_idx.at(i) << ": " << this->ftemplates.at(i) << std::endl;
    }
    for (std::size_t i=0; i<this->nparticles; i++) {
        std::cout << "  " << i << ": " << this->particles_idx.at(i) << ": " << this->fparticles.at(i) << std::endl;
    }
    for (std::size_t i=0; i<this->ntemplates; i++) {
        for (std::size_t j=0; j<this->nparticles; j++) {
            std::cout << "  " << (int)this->process_pair.at(i*this->nparticles+j) << "  ";
        }
        std::cout << std::endl;
    }
    #endif

    this->init_effective_jobs();

    #if 0
    for (std::size_t i=0; i<this->job_ntemplates; i++) {
        std::cout << "    " << i << ": " << this->job_templates_idx.at(i) << ": " << this->job_ftemplates.at(i) << std::endl;
    }
    for (std::size_t i=0; i<this->job_nparticles; i++) {
        std::cout << "    " << i << ": " << this->job_particles_idx.at(i) << ": " << this->job_fparticles.at(i) << std::endl;
    }
    for (std::size_t i=0; i<this->job_ntemplates; i++) {
        for (std::size_t j=0; j<this->job_nparticles; j++) {
            std::cout << "  " << (int)this->job_process_pair.at(i*this->job_nparticles+j) << "  ";
        }
        std::cout << std::endl;
    }

    std::cout << "in process: " << this->job_ftemplates.size() << "," << this->job_fparticles.size() << " - " << this->job_wedge_templates.size() << "," << this->job_wedge_particles.size() << std::endl;


    std::cout << (this->job_ntemplates&&this->job_nparticles || !this->job_ntemplates&&!this->job_nparticles) <<
           (this->job_ntemplates==this->job_ftemplates.size()) <<
           (this->job_ntemplates==this->job_templates_idx.size()) <<
           (this->job_ntemplates==this->job_wedge_templates.size()) <<
           (this->job_nparticles==this->job_fparticles.size()) <<
           (this->job_nparticles==this->job_particles_idx.size()) <<
           (this->job_nparticles==this->job_wedge_particles.size()) <<
           (this->job_ntemplates*this->job_nparticles==this->job_process_pair.size()) << " " <<
           this->job_ntemplates << this->job_nparticles << this->job_process_pair.size() << std::endl;
    #endif


    assert((this->job_ntemplates&&this->job_nparticles || !this->job_ntemplates&&!this->job_nparticles) &&
           this->job_ntemplates==this->job_ftemplates.size() &&
           this->job_ntemplates==this->job_templates_idx.size() &&
           this->job_ntemplates==this->job_wedge_templates.size() &&
           this->job_nparticles==this->job_fparticles.size() &&
           this->job_nparticles==this->job_particles_idx.size() &&
           this->job_nparticles==this->job_wedge_particles.size() &&
           this->job_ntemplates*this->job_nparticles==this->job_process_pair.size());

    if (this->job_ntemplates) {

        std::auto_ptr<tom::CorrelationHandler<TFLOAT> > corrhdl;
        if (this->saveccvols) {
            corrhdl.reset(new tom::CorrelationHandlerPeakVolume<TFLOAT>(job_ntemplates, job_nparticles, volsize_effective, volsize_effective, volsize_effective));
        } else {
            assert(this->get_peaks);
            if (maskcc.get()) {
                corrhdl.reset(new tom::CorrelationHandlerPeaks<TFLOAT>(job_ntemplates, job_nparticles, *maskcc, false));
            } else {
                corrhdl.reset(new tom::CorrelationHandlerPeaks<TFLOAT>(job_ntemplates, job_nparticles));
            }
        }

        this->do_corr3d(corrhdl.get());

        // Extract the results from corrhdl...

        std::size_t itemplate, iparticle, i, j;
        if (this->saveccvols) {
            assert(typeid(*corrhdl.get()) == typeid(tom::CorrelationHandlerPeakVolume<TFLOAT>) && this->output_filenames.size()==this->ntemplates*this->nparticles);

            const std::vector<std::pair<boost::shared_ptr<tom::Volume<TFLOAT> >, boost::shared_ptr<tom::Volume<int32_t> > > > &ccvols = dynamic_cast<tom::CorrelationHandlerPeakVolume<TFLOAT> &>(*corrhdl).getVols();
            const tom::Volume<TMASKCC> *maskcc = this->maskcc.get();
            std::vector<tom::st_idx> peaks;

            for (j=0, itemplate=0; itemplate<this->job_ntemplates; itemplate++) {
                for (iparticle=0; iparticle<this->job_nparticles; iparticle++, j++) {
                    assert(this->job_templates_idx[itemplate] < this->ntemplates && this->job_particles_idx[iparticle]<this->nparticles);
                    if (this->job_process_pair[j]) {
                        const std::pair<boost::shared_ptr<tom::Volume<TFLOAT> >, boost::shared_ptr<tom::Volume<int32_t> > > &ccvol = ccvols[j];
                        assert(ccvol.first.get() && ccvol.second.get() || !ccvol.first.get() && !ccvol.second.get());
                        if (ccvol.first.get() && ccvol.second.get()) {
                            const tom::Volume<TFLOAT> &vccvol = *ccvol.first;
                            const tom::Volume<int32_t> &vccidx = *ccvol.second;
                            i = this->job_templates_idx[itemplate]*this->nparticles + this->job_particles_idx[iparticle];
                            const typename tom::corr::JobManagerClient<TFLOAT, TMASKCC>::filename_triple &f = this->output_filenames[i];

                            if (this->get_peaks) {
                                assert(this->peak_list.size()==this->ntemplates*this->nparticles);
                                if (maskcc) {
                                    peaks = tom::peak(vccvol, *maskcc);
                                } else {
                                    peaks = tom::peak(vccvol);
                                }
                                if (!peaks.empty()) {
                                    tom::cc_peak<TFLOAT> &peak = peak_list[i];
                                    const std::size_t &px = peaks[0].x;
                                    const std::size_t &py = peaks[0].y;
                                    const std::size_t &pz = peaks[0].z;
                                    peak.x = px;
                                    peak.y = py;
                                    peak.z = pz;
                                    peak.val = vccvol.get(px, py, pz);
                                    peak.angle_idx = vccidx.get(px, py, pz);
                                }
                            }
                            vccvol.write_to_em(f.filename_ccvol, NULL);
                            vccidx.write_to_em(f.filename_ccidx, NULL);
                            this->angles->write_to_em(f.filename_angle, NULL);
                        }
                    }
                }
            }
        } else {
            assert(typeid(*corrhdl) == typeid(tom::CorrelationHandlerPeaks<TFLOAT>));
            const std::vector<std::vector<tom::cc_peak<TFLOAT> > > job_peaklist = dynamic_cast<tom::CorrelationHandlerPeaks<TFLOAT> *>(corrhdl.get())->getPeakList();
            for (j=0, itemplate=0; itemplate<this->job_ntemplates; itemplate++) {
                for (iparticle=0; iparticle<this->job_nparticles; iparticle++, j++) {
                    assert(this->job_templates_idx[itemplate] < this->ntemplates && this->job_particles_idx[iparticle]<this->nparticles);
                    if (job_process_pair[j]) {
                        i = this->job_templates_idx[itemplate] * this->nparticles + this->job_particles_idx[iparticle];
                        peak_list[i] = job_peaklist[itemplate][iparticle];
                    }
                }
            }
        }
    }

    this->status = STATUS_FINISHED;
}




/****************************************************************************//**
 * \brief
 *
 *******************************************************************************/
template<typename TFLOAT, typename TMASKCC>
void tom::corr::JobManagerClient<TFLOAT, TMASKCC>::do_corr3d(tom::CorrelationHandler<TFLOAT> *corrhdl) const {

    assert(this->assert_status() &&
           this->status == STATUS_PROCESS &&
           corrhdl &&
           ((!this->job_ntemplates&&!this->job_nparticles)||(this->job_ntemplates&&this->job_nparticles)) &&
           this->job_ntemplates == this->job_ftemplates.size() &&
           this->job_ntemplates == this->job_templates_idx.size() &&
           this->job_ntemplates == this->job_wedge_templates.size() &&
           this->job_nparticles == this->job_fparticles.size() &&
           this->job_nparticles == this->job_particles_idx.size() &&
           this->job_nparticles == this->job_wedge_particles.size() &&
           this->job_ntemplates*this->job_nparticles == this->job_process_pair.size());

    if (!this->job_ntemplates) {
        return;
    }


    std::size_t i;

    auto_vector<tom::Volume<TFLOAT> > vtemplates(this->job_ntemplates);
    auto_vector<tom::Volume<std::complex<TFLOAT> > > vparticles(this->job_nparticles);

    std::vector<const tom::Volume<TFLOAT> *> templates(this->job_ntemplates, NULL);
    std::vector<const tom::Volume<std::complex<TFLOAT> > *> particles(this->job_nparticles, NULL);


    const uint32_t binningv[3] = { this->binning, this->binning, this->binning };

    tom::Volume<TFLOAT> voltmp_t(this->volsize_effective, this->volsize_effective, this->volsize_effective, NULL,NULL);
    tom::Volume<std::complex<TFLOAT> > voltmp_f(this->volsize_effective, this->volsize_effective, this->volsize_effective/2+1, NULL,NULL);
    tom::fftw::Plan<TFLOAT> plan(voltmp_t, voltmp_f, FFTW_DESTROY_INPUT | this->fftw_flag);

    std::vector<double> variance_particles(this->job_nparticles, 0);

    tom::Volume<TFLOAT> *pv;
    tom::Volume<std::complex<TFLOAT> > *pv_f;
    std::auto_ptr<tom::Volume<TFLOAT> > pv_auto;
    tom_io_em_header header;
    std::ostringstream ss;
    double mean, variance;

    #ifdef NDEBUG
    #   define CORR3D_VERBOSE_DEBUG(x) static_cast<void>(0)
    #else
    #   define CORR3D_VERBOSE_DEBUG(x) x
    #endif

    #undef CORR3D_VERBOSE_DEBUG
    #define CORR3D_VERBOSE_DEBUG(x) static_cast<void>(0)
    CORR3D_VERBOSE_DEBUG( std::cerr << "DEBUG-OUTPUT: << " << std::endl; );

    CORR3D_VERBOSE_DEBUG( static int CNT = 0; );


    // Read the templates from disk.
    for (i=0; i<this->job_ntemplates; i++) {
        try {
            CORR3D_VERBOSE_DEBUG( std::cerr << "template " << i << ": " << this->job_ftemplates[i] << std::endl; );
            tom::read_from_em<TFLOAT>(pv, this->job_ftemplates[i], NULL,NULL,binningv, &header, NULL,NULL);
            pv_auto.reset(pv);
            //CORR3D_VERBOSE_DEBUG( pv->write_to_em(this->outputdir+"/templ"+helper::stringify(i)+"_orig.em", NULL); );
        } catch (...) {
            if (this->force_files_exist) {
                ss << "Error reading file \"" << this->job_ftemplates[i] << "\" (force_files_exist)";
                throw std::invalid_argument(ss.str());
            }
            continue;
        }
        if (!pv->is_equal_size(this->volsize_effective)) {
            if (this->force_files_exist) {
                ss << "File \"" << this->job_ftemplates[i] << "\" has wrong size. (force_files_exist)";
                throw std::invalid_argument(ss.str());
            }
            continue;
        }

        // Normalise to mean 0 and std 1.
        pv->stat(mean, variance, false);
        if (!variance) {
            if (this->force_files_exist) {
                ss << "File \"" << this->job_ftemplates[i] << "\" has variance zero. (force_files_exist)";
                throw std::invalid_argument(ss.str());
            }
            continue;
        }
        pv->template shift_scale<TFLOAT>(-mean, 1);


        CORR3D_VERBOSE_DEBUG( pv->write_to_em(this->outputdir+"/templ"+helper::stringify(CNT)+"_orig_norm.em", NULL); );


        vtemplates.assign_direct(i, pv_auto.release());
        templates[i] = pv;
    }

    const tom::Volume<TFLOAT> *mask_sphere = this->mask_sphere.get();

    CORR3D_VERBOSE_DEBUG( std::cerr << "spheremask: " << (void *)mask_sphere << std::endl; );
    //CORR3D_VERBOSE_DEBUG( if (mask_sphere) { mask_sphere->write_to_em(this->outputdir+"/mask_sphere.em", NULL); } );


    // Read the particles from disk.
    for (i=0; i<this->job_nparticles; i++) {
        try {
            CORR3D_VERBOSE_DEBUG( std::cerr << "particle " << i << ": " << this->job_fparticles[i] << std::endl; );
            tom::read_from_em<TFLOAT>(pv, this->job_fparticles[i], NULL,NULL,binningv, &header, NULL,NULL);
            pv_auto.reset(pv);
            //CORR3D_VERBOSE_DEBUG( pv->write_to_em(this->outputdir+"/parti"+helper::stringify(i)+"_orig.em", NULL); );
        } catch (...) {
            if (this->force_files_exist) {
                ss << "Error reading file \"" << this->job_fparticles[i] << "\" (force_files_exist)";
                throw std::invalid_argument(ss.str());
            }
            continue;
        }
        if (!pv->is_equal_size(this->volsize_effective)) {
            if (this->force_files_exist) {
                ss << "File \"" << this->job_fparticles[i] << "\" has wrong size. (force_files_exist)";
                throw std::invalid_argument(ss.str());
            }
            continue;
        }
        voltmp_t.setValues(*pv);


        // Applay mask to particle and/or normalise
        if (mask_sphere) {
            tom::norm_mask<TFLOAT, TFLOAT, double>(voltmp_t, *mask_sphere, tom::norm::NORM_NO_NORM, &variance_particles[i], false);
        } else {
            variance_particles[i] = voltmp_t.variance(false);
        }
        if (!variance_particles[i]) {
            if (this->force_files_exist) {
                ss << "File \"" << this->job_fparticles[i] << "\" has variance zero. (force_files_exist)";
                throw std::invalid_argument(ss.str());
            }
            continue;
        }

        //CORR3D_VERBOSE_DEBUG( voltmp_t.write_to_em(this->outputdir+"/parti"+helper::stringify(i)+"_masked.em", NULL); );

        plan.execute(voltmp_t, voltmp_f);

        // Set the mean to zero (in fourier space).
        voltmp_f.get(0,0,0) = std::complex<TFLOAT>(0,0);

        //CORR3D_VERBOSE_DEBUG( tom::Volume<TFLOAT>(voltmp_f, true ).write_to_em(this->outputdir+"/parti"+helper::stringify(CNT)+"_fd_r.em", NULL); );
        //CORR3D_VERBOSE_DEBUG( tom::Volume<TFLOAT>(voltmp_f, false).write_to_em(this->outputdir+"/parti"+helper::stringify(CNT)+"_fd_c.em", NULL); );


        vparticles.assign_direct(i, pv_f = new tom::Volume<std::complex<TFLOAT> >(voltmp_f));
        particles[i] = pv_f;
    }

    //std::cout << "before corr3d_single: " << templates.size() << "," << particles.size() << " - " << this->job_wedge_templates.size() << "," << this->job_wedge_particles.size() << std::endl;

    tom::corr3d_single(templates, particles, variance_particles, *this->angles, 0, 0, &this->job_process_pair[0],
                    mask_sphere, &this->job_wedge_templates, &this->job_wedge_particles, corrhdl, this->fftw_flag, this->reapply_mask_after_wedge);


    #undef CORR3D_VERBOSE_DEBUG
}





/****************************************************************************//**
 * Writes the current configuration to an output stream.
 *******************************************************************************/
void tom::corr::ConfigDataClient::write_to_log(std::ostream &stream) const {

    stream <<   "# \n"
                "# ConfigDataClient: \n";
    if (!assert_status()) {
        stream << "# ConfigDataClient is in an inconsistent state!!!!\n";
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
                "outputdir="                    << outputdir                    << "\n"
                "return_peak_to_root="          << return_peak_to_root          << "\n"
                "saveccvols="                   << saveccvols                   << "\n"
                "resume="                       << resume                       << "\n"
                "force_files_exist="            << force_files_exist            << "\n"
                "nice="                         << nice                         << "\n"
                "#"                                                             << std::endl;

}




// TEMPLATE INSTANTIATIONS

template class tom::corr::JobManagerClient<float , char>;
template class tom::corr::JobManagerClient<double, char>;
