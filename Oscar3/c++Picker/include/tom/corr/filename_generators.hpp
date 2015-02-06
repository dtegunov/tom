/****************************************************************************//**
 * \file filename_generators.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    19.02.2008
 *******************************************************************************/
#ifndef ___INCLUDE_FILENAME_GENERATORS_HPP__
#define ___INCLUDE_FILENAME_GENERATORS_HPP__



#include <string>
#include <helper/filesystem.hpp>



namespace tom {

// forward declaration.
template<typename T> class Volume;


namespace av4 {

class filename_generator {
public:
    filename_generator(const std::string &outputdir, const std::string &prefix);
    std::string get_peakfilename(std::size_t iit) const;
    std::string get_outputdir(std::size_t iit) const;
    std::string get_config_correlation(std::size_t iit) const;
    std::string get_logfilename_correlation(std::size_t iit) const;
    std::string get_logfilename_average(std::size_t iit) const;
    std::string get_templates_list(std::size_t iit) const;
    std::string get_particles_list_correlation(std::size_t iit) const;
    std::string get_particles_list_average(std::size_t iit) const;
    std::string get_angles_list(std::size_t iit) const;
    std::string get_angles_file(std::size_t iit, const tom::Volume<double> &v) const;
private:
    std::string outputdir;
    std::string prefix;
}; // class filename_generator

} // namespace av4






namespace avg {
class filename_generator {
public:
    filename_generator(const std::string &outputdir);
    std::string get_avg_s(std::size_t itemplate) const;
    std::string get_avgwedge_s(std::size_t itemplate) const;
    std::string get_avg_1(std::size_t itemplate) const;
    std::string get_avgwedge_1(std::size_t itemplate) const;
    std::string get_avg_sqrtn(std::size_t itemplate) const;
    std::string get_avgwedge_sqrtn(std::size_t itemplate) const;
    std::string get_avg_n(std::size_t itemplate) const;
    std::string get_avgwedge_n(std::size_t itemplate) const;
private:
    std::string outputdir;
}; // class filename_generator
} // namespace avg



namespace corr {
class filename_generator {
public:
    filename_generator(const std::string &outputdir);
    std::string get_ccv(std::size_t itemplate, std::size_t iparticle) const;
    std::string get_cci(std::size_t itemplate, std::size_t iparticle) const;
    std::string get_ang(std::size_t itemplate, std::size_t iparticle) const;

private:
    helper::fs::path p;
    std::string outputdir;
}; // class filename_generator
} // namespace corr


} // namespace tom



#endif


