/****************************************************************************//**
 * \file mpi_fcn.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    08.01.2007
 *******************************************************************************/
#ifndef ___INCLUDE_MPI__MPI_FCN_HPP__
#define ___INCLUDE_MPI__MPI_FCN_HPP__


#include <mpi.h>
#include <vector>
#include <boost/shared_ptr.hpp>


#include "tom/core/wedge.hpp"
#include "tom/core/cc_peak.hpp"
#include "tom/os3/os3_structures.hpp"




namespace tom {
namespace mpi {





void startup(std::vector<std::string> &vargv, int argc, char *argv[]);
int load_fftw_wisdom(const std::string &fftw_wisdom_dir, MPI_Comm comm_world);
int save_fftw_wisdom(MPI_Comm comm_world);


void bcast_cmd_line(std::vector<std::string> &cmd, int argc, char *argv[]);
int bcast_string(int root, MPI_Comm comm, std::string &s);
template<typename T> int bcast_volume(int root, MPI_Comm comm, tom::Volume<T> *&v);

int send_vstring(int dst, int tag, MPI_Comm comm, const std::vector<std::string> &v);
int recv_vstring(int src, int tag, MPI_Comm comm, std::vector<std::string> &v);

int os3_send_job_request(int dst, int tag, MPI_Comm comm,int senderRank);
int os3_recieve_job_request(MPI_Comm comm,int& senderRank);
int send_os3_job(int dst, int tag, MPI_Comm comm, const tom::os3_job& job);
int recv_os3_job(int src, int tag, MPI_Comm comm, tom::os3_job& job);

template<typename T> int send_wedges(int dst, int tag, MPI_Comm comm, const std::vector<tom::WedgeDescriptor<T> *> &wedges);
template<typename T> int recv_wedges(int src, int tag, MPI_Comm comm, std::vector<boost::shared_ptr<tom::Wedge<T> > > &wedges);


template<typename T> int send_volume(int dst, int tag, MPI_Comm comm, const tom::Volume<T> &v);
template<typename T> int recv_volume(int src, int tag, MPI_Comm comm, std::auto_ptr<tom::Volume<T> > &v);


/// Gets the MPI_Datatype from an template-function.
template<typename T> inline MPI_Datatype get_MPI_Type();

template<typename TFLOAT> MPI_Datatype create_cc_peak_type();
template<typename TFLOAT> void send_peaklist(const std::vector<tom::cc_peak<TFLOAT> > &peak_list, int dst, int msg_tag, MPI_Comm comm);
template<typename TFLOAT> void recv_peaklist(std::vector<tom::cc_peak<TFLOAT> > &peak_list, int src, int msg_tag, MPI_Comm comm);


int gather_hostname(MPI_Comm comm, std::vector<std::pair<pid_t, std::string> > &res, int root=-1);

} // namespace tom::mpi
} // namespace tom





// INLINE functions.

namespace tom { namespace mpi {
template<          > inline MPI_Datatype get_MPI_Type<float        >() { return MPI_FLOAT; }
template<          > inline MPI_Datatype get_MPI_Type<double       >() { return MPI_DOUBLE; }
template<          > inline MPI_Datatype get_MPI_Type<int          >() { return MPI_INT; }
template<          > inline MPI_Datatype get_MPI_Type<char         >() { return MPI_CHAR; }
template<          > inline MPI_Datatype get_MPI_Type<unsigned long>() { return MPI_UNSIGNED_LONG; }
template<          > inline MPI_Datatype get_MPI_Type<unsigned     >() { return MPI_UNSIGNED; }
template<          > inline MPI_Datatype get_MPI_Type<long         >() { return MPI_LONG; }
template<          > inline MPI_Datatype get_MPI_Type<unsigned char>() { return MPI_UNSIGNED_CHAR; }
#ifdef __USE_MPI_BOOL
// MPI::BOOL (and C++ language binding) is MPI2 only. (thus problems compiling with pure MPI1.1
template<          > inline MPI_Datatype get_MPI_Type<bool>() { return MPI::BOOL; }
#endif
}}






#endif
















