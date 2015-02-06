/****************************************************************************//**
 * \file c_exec.hpp
 * \brief Class for a triple similar to std::pair.
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    30.01.2008
 *******************************************************************************/
#ifndef ___INCLUDE_C_EXEC_HPP__
#define ___INCLUDE_C_EXEC_HPP__



namespace tom {
namespace mpi {




/****************************************************************************//**
 *
 *******************************************************************************/
class exec_params {
public:
    virtual ~exec_params() { }
    virtual void execvp(std::string &file, std::vector<std::string> &argv) const = 0;
}; // class exec_params



/****************************************************************************//**
 *
 *******************************************************************************/
class exec_params_mpirun: public exec_params {
public:
    exec_params_mpirun(int np): np(np) { if (np < 1) { throw std::runtime_error("number of processes is less then 1."); } };
    virtual ~exec_params_mpirun() { }
    virtual void execvp(std::string &file, std::vector<std::string> &argv) const;

private:
    int np;
}; // class exec_params_mpirun



/****************************************************************************//**
 *
 *******************************************************************************/
class exec_params_poe: public exec_params {
public:
    exec_params_poe(int np): np(np) { if (np < 1) { throw std::runtime_error("number of processes is less then 1."); } };
    virtual ~exec_params_poe() { }
    virtual void execvp(std::string &file, std::vector<std::string> &argv) const;

private:
    int np;
}; // class exec_params_poe


std::auto_ptr<tom::mpi::exec_params> create_from_config_parameter(const std::string &s, int np);


} // namespace mpi
} // namespace tom





// Inline functions.
#include <unistd.h>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <algorithm>
#include <ctype.h>


/****************************************************************************//**
 *
 *******************************************************************************/
void tom::mpi::exec_params_mpirun::execvp(std::string &file, std::vector<std::string> &argv) const {

    std::vector<std::string> argv_new;


    argv_new.push_back("mpirun");
    argv_new.push_back("-np");
    argv_new.push_back(boost::lexical_cast<std::string>(np));

    argv_new.push_back(file);

    // Skip the first argument because its supposed to be the name of the programm itself.
    std::size_t i = 1;
    while (i < argv.size()) {
        argv_new.push_back(argv[i]);
        i++;
    }

    argv.swap(argv_new);
    file = "mpirun";
}




/****************************************************************************//**
 *
 *******************************************************************************/
void tom::mpi::exec_params_poe::execvp(std::string &file, std::vector<std::string> &argv) const {

    std::vector<std::string> argv_new;

    argv_new.push_back("poe");
    argv_new.push_back(file);
    // Skip the first argument because its supposed to be the name of the programm itself.
    std::size_t i = 1;
    while (i < argv.size()) {
        argv_new.push_back(argv[i]);
        i++;
    }
    argv_new.push_back("-procs");
    argv_new.push_back(boost::lexical_cast<std::string>(np));

    argv.swap(argv_new);
    file = "poe";
}


/****************************************************************************//**
 *
 *******************************************************************************/
std::auto_ptr<tom::mpi::exec_params> create_from_config_parameter(const std::string &s_, int np) {
    std::auto_ptr<tom::mpi::exec_params> res;

    std::string s(s_);
    std::transform(s.begin(), s.end(), s.begin(), &tolower);

    if (s == "poe") {
        res.reset(new tom::mpi::exec_params_poe(np));
    } else if (s == "mpirun") {
        res.reset(new tom::mpi::exec_params_mpirun(np));
    } else { // default
        throw std::invalid_argument("The Parallel Operating Environment '" + s_ + "' is not recognized. Try one of 'poe', 'mpirun'.");
    }

    return res;
}


#endif




