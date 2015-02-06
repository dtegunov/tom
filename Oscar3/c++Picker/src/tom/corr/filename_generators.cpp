/****************************************************************************//**
 * \file filename_generators.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    07.01.2008
 *******************************************************************************/
#include <tom/corr/filename_generators.hpp>



#include <sstream>
#include <iomanip>

#include <boost/md5.hpp>


#include <helper/filesystem.hpp>

#include <tom/core/volume.hpp>




/****************************************************************************//**
 *
 *******************************************************************************/
tom::av4::filename_generator::filename_generator(const std::string &outputdir_, const std::string &prefix_)
    : outputdir(outputdir_),
      prefix(prefix_) {
    if (outputdir.empty()) {
        outputdir = ".";
    }
    if (!prefix.empty() && prefix[prefix.size()-1]!='_') {
        prefix.push_back('_');
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_peakfilename(std::size_t iit) const {

    std::ostringstream ss;
    ss << prefix << std::setw(3) << std::setfill('0') << iit << "_peakfile";

    helper::fs::path p(outputdir);
    p /= ss.str();

    return p.string();

}


/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_outputdir(std::size_t iit) const {

    std::ostringstream ss;
    ss << prefix << std::setw(3) << std::setfill('0') << iit;

    helper::fs::path p(outputdir);
    p /= ss.str();

    return p.string();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_templates_list(std::size_t iit) const {

    std::ostringstream ss;
    //ss << prefix << std::setw(3) << std::setfill('0') << iit << "_templates_list.txt";
    ss << "templates_list.txt";

    helper::fs::path p(get_outputdir(iit));
    p /= ss.str();

    return p.string();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::av4::filename_generator::get_angles_file(std::size_t iit, const tom::Volume<double> &v) const {

    std::auto_ptr<tom::Volume<double> > v_contigous;
    const void *data;
    if (v.isContiguous()) {
        data = &v.get();
    } else {
        v_contigous.reset(new tom::Volume<double>(v));
        assert(v_contigous->isContiguous() && *v_contigous==v);
        data = &v_contigous->get();
    }


    boost::md5 csum(data, v.numel()*sizeof(double));

    std::ostringstream ss;
    //ss << prefix << std::setw(3) << std::setfill('0') << iit << "_templates_list.txt";
    ss << "angles_" << csum.digest().hex_str_value() << ".em";

    helper::fs::path p(get_outputdir(iit));
    p /= ss.str();

    return p.string();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_angles_list(std::size_t iit) const {

    std::ostringstream ss;
    //ss << prefix << std::setw(3) << std::setfill('0') << iit << "_angles_list.txt";
    ss << "angles.txt";

    helper::fs::path p(get_outputdir(iit));
    p /= ss.str();

    return p.string();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_config_correlation(std::size_t iit) const {

    std::ostringstream ss;
    ss << "config.txt";

    helper::fs::path p(get_outputdir(iit));
    p /= ss.str();

    return p.string();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_particles_list_average(std::size_t iit) const {

    std::ostringstream ss;
    //ss << prefix << std::setw(3) << std::setfill('0') << iit << "_particles_list.txt";
    ss << "particles_list_average.txt";

    helper::fs::path p(get_outputdir(iit));
    p /= ss.str();

    return p.string();
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_particles_list_correlation(std::size_t iit) const {

    std::ostringstream ss;
    //ss << prefix << std::setw(3) << std::setfill('0') << iit << "_particles_list.txt";
    ss << "particles_list_correlation.txt";

    helper::fs::path p(get_outputdir(iit));
    p /= ss.str();

    return p.string();
}




/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_logfilename_correlation(std::size_t iit) const {

    std::ostringstream ss;
    ss << prefix << std::setw(3) << std::setfill('0') << iit << "_log_correlation";

    helper::fs::path p(outputdir);
    p /= ss.str();

    return p.string();
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::string  tom::av4::filename_generator::get_logfilename_average(std::size_t iit) const {

    std::ostringstream ss;
    ss << prefix << std::setw(3) << std::setfill('0') << iit << "_log_average";

    helper::fs::path p(outputdir);
    p /= ss.str();

    return p.string();
}













/****************************************************************************//**
 *
 *******************************************************************************/
tom::avg::filename_generator::filename_generator(const std::string &outputdir_)
    : outputdir(outputdir_) {
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avg_s(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avg_s_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avg_1(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avg_1_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avg_sqrtn(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avg_sqrtn_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avg_n(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avg_n_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}

/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avgwedge_s(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avgwedge_s_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avgwedge_1(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avgwedge_1_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avgwedge_sqrtn(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avgwedge_sqrtn_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::avg::filename_generator::get_avgwedge_n(std::size_t itemplate) const {
    std::ostringstream ss;
    ss << "avgwedge_n_" << std::setfill('0') << std::setw(3) << itemplate << ".em";
    helper::fs::path p(outputdir);
    p /= ss.str();
    return p.string();
}



/****************************************************************************//**
 *
 *******************************************************************************/
tom::corr::filename_generator::filename_generator(const std::string &outputdir)
    : p(outputdir), outputdir(outputdir) {
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::corr::filename_generator::get_ccv(std::size_t itemplate, std::size_t iparticle) const {
    std::ostringstream ss;
    ss.fill('_');
    ss << "ccv_" << std::setw(4) << itemplate << std::setw(5) << iparticle << ".em";
    return (p / ss.str()).string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::corr::filename_generator::get_cci(std::size_t itemplate, std::size_t iparticle) const {
    std::ostringstream ss;
    ss.fill('_');
    ss << "cci_" << std::setw(4) << itemplate << std::setw(5) << iparticle << ".em";
    return (p / ss.str()).string();
}
/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::corr::filename_generator::get_ang(std::size_t itemplate, std::size_t iparticle) const {
    std::ostringstream ss;
    ss.fill('_');
    ss << "ang_" << std::setw(4) << itemplate << std::setw(5) << iparticle << ".em";
    return (p / ss.str()).string();
}






