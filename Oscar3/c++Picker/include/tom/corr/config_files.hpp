/****************************************************************************//**
 * \file config_files.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    06.02.2008
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__CONFIG_FILES_HPP__
#define ___INCLUDE_CORR__CONFIG_FILES_HPP__


#include <vector>
#include <boost/shared_ptr.hpp>

#include "tom/core/wedge.hpp"
#include "tom/core/cc_peak.hpp"

#include "tom/corr/jobmanager_server.hpp"

#include "helper/triple.hpp"


namespace tom {


namespace av4 {
class filename_generator;
}



namespace avg {

struct st_average_result {
public:
    std::vector<std::size_t> use_idx;
    boost::shared_ptr<double> ccval_s;
    boost::shared_ptr<double> ccval_1;
    boost::shared_ptr<double> ccval_sqrtn;
    boost::shared_ptr<double> ccval_n;
    std::vector<std::pair<double, double> > fsc;
    st_average_result()
        : use_idx(),
          ccval_s(),
          ccval_1(),
          ccval_sqrtn(),
          ccval_n(),
          fsc() {
    }
}; // struct st_average_result


void write_average_to_log(  std::ostream &clog,
                            const std::vector<std::string> &ftemplates,
                            const std::vector<std::string> &fparticles,
                            const std::vector<tom::avg::st_average_result> &av,
                            const std::string &outputdir);

void parse_average(         const std::string &filename,
                            std::vector<tom::avg::st_average_result> &av);


} // namespace tom::avg


namespace corr {




void parse_angleslist(const std::string &filename, std::size_t ntemplates, std::size_t nparticles, std::vector<boost::shared_ptr<tom::Volume<double> > > &angles);
void write_angleslist(const std::string &filename, std::size_t ntemplates, std::size_t nparticles, const std::vector<boost::shared_ptr<tom::Volume<double> > > &angles, const tom::av4::filename_generator &f_gen, unsigned long iit);


template<typename T>
void parse_filelist(const std::string &filename, std::vector<std::string> &ffilenames, std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > > &wedge_desc, std::size_t binning);



template<typename T>
void write_filelist(const std::string &filename, const std::vector<std::string> &ffilenames, const std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > > &wedge_desc);



template<typename T>
void parse_peaklist(const std::string &peakfilename, std::vector<tom::cc_peak<T> > &peak_list, std::vector<helper::triple<double, double, double> > &angles, std::size_t ntemplates, std::size_t nparticles);




}
}




#endif

