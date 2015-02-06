/****************************************************************************//**
 * \file av4.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    14.02.2008
 *******************************************************************************/

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "tom/corr/av4manager.hpp"


#include "helper/snippets.hpp"
#include "helper/ostream_swapper.hpp"









/****************************************************************************//**
 *
 *******************************************************************************/
int main(int argc, char *argv[]) {

    bool show_help = false;
    for (int i=1; i<argc; i++) {
        if (!strcmp(argv[i], "-h") ||
            !strcmp(argv[i], "--help")) {
            show_help = true;
            break;
        }
    }

    if (show_help || argc<2 || argc>3) {
        std::map<std::string, std::string> parameter_description;
        std::map<std::string, std::string>::const_iterator it;
        const std::size_t nlength = 80;

        std::cerr << "Usage: " << (argc?argv[0]:"av4") << " CONFIG_FILE [LOG_FILE]"        "\n"
                     "compiled at " __DATE__ " " __TIME__ " from " __FILE__ << std::endl;

        if (show_help) {
            std::cerr << std::endl << "PARAMETERS: " << std::endl;
            tom::av4::manager::get_parameter_description().swap(parameter_description);
            for (it=parameter_description.begin(); it!=parameter_description.end(); it++) {
                helper::print_parameters(std::cerr, it->first, it->second, "  ", nlength);
            }

            std::cerr << std::endl << "PARAMETERS STATIC: " << std::endl;
            tom::av4::manager_static::get_parameter_description().swap(parameter_description);
            for (it=parameter_description.begin(); it!=parameter_description.end(); it++) {
                helper::print_parameters(std::cerr, it->first, it->second, "  ", nlength);
            }

            std::cerr << std::endl << "PARAMETERS ANGULAR_REFINEMENT: " << std::endl;
            tom::av4::manager_angular_refinement::get_parameter_description().swap(parameter_description);
            for (it=parameter_description.begin(); it!=parameter_description.end(); it++) {
                helper::print_parameters(std::cerr, it->first, it->second, "  ", nlength);
            }
        } else {
            std::cerr << "Try `" << (argc>0 ? argv[0] : "av4") << " --help' for more information." << std::endl;
        }
        return -1;
    }

    std::stringstream clog;
    std::auto_ptr<std::ofstream> flog;
    std::auto_ptr<helper::ostream_swapper> stream_swapper;
    if (argc<3 || std::string(argv[2])=="-") {
        stream_swapper.reset(new helper::ostream_swapper(clog, std::cout));
    } else {
        flog.reset(new std::ofstream(argv[2]));
        if (!flog->good()) {
            std::cerr << "could not open log file \"" << argv[2] << "\" for writing." << std::endl;
            return -2;
        }
        stream_swapper.reset(new helper::ostream_swapper(clog, *flog));
    }
    clog.precision(20);
    clog << "# start " << (argc>0?argv[0]:"av4") << " at " << helper::now2str() << " (compiled at " __DATE__ " " __TIME__ ")." << std::endl;
    {
        std::vector<char> cbuff(HOST_NAME_MAX+2, 0);
        int i = gethostname(&cbuff[0], cbuff.size()-1);
        clog << "# process id " << getpid() << " on " << (i==0 ? &cbuff[0] : "UNKNOWN") << std::endl;
    }


    const std::string config_file = argv[1];
    std::auto_ptr<tom::av4::manager> c;
    clog << "# Parsing config file \"" << config_file << "\"." << std::endl;
    try {
        c.reset(tom::av4::manager::create(config_file, clog));
    } catch (std::exception &e) {
        std::cerr << "Error parsing the config file \"" << config_file << "\": " << e.what() << std::endl;
        clog << "Error parsing the config file \"" << config_file << "\": " << e.what() << std::endl;
        return -2;
    }

    clog << "# CONFIGURATION: \n\n" << c->getConfigStr() << "\n";

    try {
        c->process();
    } catch (std::exception &e) {
        clog << "# PROCESSING ABORTED DUE TO EXCEPTION: " << e.what() << "." << std::endl;
        throw;
    } catch (...) {
        clog << "# PROCESSING ABORTED DUE TO UNKNOWN EXCEPTION." << std::endl;
        throw;
    }

    clog << "# finished at " << helper::now2str() << std::endl;
    return 0;
}



















