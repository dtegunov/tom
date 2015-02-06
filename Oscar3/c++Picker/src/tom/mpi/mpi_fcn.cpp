/****************************************************************************//**
 * \file mpi_fcn.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    08.01.2008
 *******************************************************************************/
#include "tom/mpi/mpi_fcn.hpp"


#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <assert.h>
#include <typeinfo>


#include <boost/static_assert.hpp>


#include "tom/core/fftw_plan.hpp"


#ifndef THREAD_SAFE
#  error Define THREAD_SAFE.
#endif




#if THREAD_SAFE
#  define __MAKE_STATIC__
#else
#  define __MAKE_STATIC__ static
#endif


namespace {

/****************************************************************************//**
 * This class creates an MPI_Datatype for tom::cc_peak.
 *
 * The type is only created upon the first access to the data-type, not uppon
 * creation. This is done, to allow instanciations of the class at file scope.
 * In that case the constructor it is executed at program start, but before
 * initialising the type, MPI_Init must be called.
 *******************************************************************************/
template<typename TFLOAT>
class MPI_Type_cc_peak {
public:
    MPI_Type_cc_peak(): type(MPI_DATATYPE_NULL) {
        static tom::cc_peak<TFLOAT> peak0[2];
        {
            // Check that MPI_Init were called.
            int flag;
            MPI_Initialized(&flag);
            if (!flag) {
                throw std::runtime_error("MPI_Init not called before accessing MPI_Type_cc_peak.");
            }
        }
        MPI_Datatype type[6] = {    tom::mpi::get_MPI_Type<typename tom::cc_peak<TFLOAT>::idx_type>(),
                                    tom::mpi::get_MPI_Type<typename tom::cc_peak<TFLOAT>::idx_type>(),
                                    tom::mpi::get_MPI_Type<typename tom::cc_peak<TFLOAT>::idx_type>(),
                                    tom::mpi::get_MPI_Type<typename tom::cc_peak<TFLOAT>::angle_idx_type>(),
                                    tom::mpi::get_MPI_Type<typename tom::cc_peak<TFLOAT>::val_type>(),
                                    MPI_UB };
        int blocklen[6] = { 1, 1, 1, 1, 1, 1 };
        MPI_Aint disp[6];
        MPI_Aint base;
        int i;
        MPI_Address(&peak0[0].x,         &disp[0]);
        MPI_Address(&peak0[0].y,         &disp[1]);
        MPI_Address(&peak0[0].z,         &disp[2]);
        MPI_Address(&peak0[0].angle_idx, &disp[3]);
        MPI_Address(&peak0[0].val,       &disp[4]);
        MPI_Address(&peak0[1].x,         &disp[5]);

        base = disp[0];
        for (i=0; i<6; i++) { disp[i] -= base; }

        MPI_Type_struct(6, blocklen, disp, type, &this->type);
        MPI_Type_commit(&this->type);
    }
    ~MPI_Type_cc_peak() {
        MPI_Type_free(&this->type);
    }

    /// Gets a copy of the pointer to the type in the class. Never change or free it.
    MPI_Datatype get_type() {
        return this->type;
    };

    /// The result of \c clone_type must be freed using \c MPI_Type_free.
    MPI_Datatype clone_type() {
        MPI_Datatype t2;
        MPI_Type_contiguous(1, this->type, &t2);
        MPI_Type_commit(&t2);
        return t2;
    };


private:
    MPI_Type_cc_peak(const MPI_Type_cc_peak &c): type(c.clone_type()) { } // hide copyctor
    MPI_Datatype type;
};
}


/****************************************************************************//**
 * \brief create MPI_Datatype for an array of type tom::cc_peak.
 *
 * The result of the function must be freed calling MPI_Type_free. It is already
 * commited.
 *******************************************************************************/
template<typename TFLOAT>
MPI_Datatype tom::mpi::create_cc_peak_type() {
    return ::MPI_Type_cc_peak<TFLOAT>().clone_type();
}

/****************************************************************************//**
 * \brief Send a vector of peaks through MPI.
 *
 * \param[in] peak_list A vector containing the tom::cc_peak<>s to send.
 * \param[in] dst The rank of the destination process in the MPI-comminicator.
 * \param[in] msg_tag The tag used for the send operation.
 * \param[in] comm The communicater used for sending.
 *
 * This function is the sender side operation of tom::mpi::recv_peaklist.
 *******************************************************************************/
template<typename TFLOAT>
void tom::mpi::send_peaklist(const std::vector<tom::cc_peak<TFLOAT> > &peak_list, int dst, int msg_tag, MPI_Comm comm) {
    // The count (i.e. peak_list.size()) may also be zero.

    #if 0
    std::cerr << "send peaklist: " << peak_list.size() << std::endl;
    int i = 0;
    typename std::vector<tom::cc_peak<TFLOAT> >::const_iterator it;
    for (it=peak_list.begin(); it!=peak_list.end(); it++,i++) {
        std::cerr << i << ": [" << it->x << "," << it->y << "," << it->z << "] " << std::endl;
    }
    #endif

    //std::cout << "HELLO WORLD: " << __FILE__ << ":" << __LINE__ << std::endl;

    MPI_Send(   const_cast<tom::cc_peak<TFLOAT> *>(&peak_list[0]),
                peak_list.size(),
                ::MPI_Type_cc_peak<TFLOAT>().get_type(),
                dst, msg_tag, comm);

    // std::cout << "HELLO WORLD: " << __FILE__ << ":" << __LINE__ << std::endl;
}

/****************************************************************************//**
 * \brief Receives a vector of peaks through MPI.
 *
 * \param[out] peak_list A vector containing the tom::cc_peak<>s after receiving.
 * \param[in] src The rank of the destination process in the MPI-comminicator.
 * \param[in] msg_tag The tag used for the send operation.
 * \param[in] comm The communicater used for sending.
 *
 * This function is the receiver side operation of tom::mpi::send_peaklist.
 *******************************************************************************/
template<typename TFLOAT>
void tom::mpi::recv_peaklist(std::vector<tom::cc_peak<TFLOAT> > &peak_list, int src, int msg_tag, MPI_Comm comm) {

    // std::cout << "HELLO WORLD: " << __FILE__ << ":" << __LINE__ << std::endl;

    MPI_Type_cc_peak<TFLOAT> c;
    MPI_Datatype peak_type = c.get_type();

    MPI_Status status;
    int count;

    MPI_Probe(src, msg_tag, comm, &status);
    MPI_Get_count(&status, peak_type, &count);

    peak_list.resize(count);
    MPI_Recv(&peak_list[0], count, peak_type, src, msg_tag, comm, &status);

    // std::cout << "HELLO WORLD: " << __FILE__ << ":" << __LINE__ << std::endl;
    #if 0
    std::cerr << "recieved peaklist: " << peak_list.size() << std::endl;
    int i = 0;
    typename std::vector<tom::cc_peak<TFLOAT> >::iterator it;
    for (it=peak_list.begin(); it!=peak_list.end(); it++,i++) {
        std::cerr << i << ": [" << it->x << "," << it->y << "," << it->z << "] " << std::endl;
    }
    #endif


}



/****************************************************************************//**
 * \brief Checks for the MPI_IO attribute and broadcast the input arguments.
 *
 * Should be called soon after MPI_Init
 *******************************************************************************/
void tom::mpi::startup(std::vector<std::string> &vargv, int argc, char *argv[]) {

    #if 0
    /* Check MPI_IO for io-capabilities. */
    int *mpi_io_ptr;
    int flag;
    MPI_Attr_get(MPI_COMM_WORLD, MPI_IO, &mpi_io_ptr, &flag);

    if (!flag) {
        std::cout << "Attribute MPI_IO: not defined (not MPI-compilant)!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    } else if (*mpi_io_ptr == MPI_PROC_NULL) {
        std::cout << "Attribute MPI_IO: MPI_PROC_NULL: program needs file io. EXIT" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -2);
    } else if (*mpi_io_ptr == MPI_ANY_SOURCE) {
    } else if (0) {
        std::cout << "Attribute MPI_IO: returned " << *mpi_io_ptr << ": program needs file-io for MPI_ANY_SOURCE. EXIT" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -3);
    }
    #endif

    /* Check that command line options are the same on every process. */
    tom::mpi::bcast_cmd_line(vargv, argc, argv);
}



/****************************************************************************//**
 * \brief Broadcast the command line arguments to all processes.
 *
 * \param[out] cmd Vector containing the command line arguments.
 * \param[in] argc The number of arguments.
 * \param[in] argv The command line arguments.
 *
 * For implementations of MPI where not each process gets a copy of the
 * command line.
 * It uses MPI_COMM_WORLD as communicater and does no error checking of the
 * MPI-Functions.
 *******************************************************************************/
void tom::mpi::bcast_cmd_line(std::vector<std::string> &cmd, int argc, char *argv[]) {

    struct {
        int argc;
        int rank;
    } sendbuf, recvbuf;

    int my_rank_world;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_world);

    sendbuf.argc = argc;
    sendbuf.rank = my_rank_world;
    MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

    int all_the_same;
    int has_max = argc==recvbuf.argc;
    MPI_Allreduce(&has_max, &all_the_same, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    cmd.resize(recvbuf.argc);
    if (all_the_same) {
        for (int i=0; i<recvbuf.argc; i++) {
            cmd.at(i) = argv[i];
        }
    } else {
        int slen;
        if (my_rank_world == recvbuf.rank) {
            for (int i=0; i<recvbuf.argc; i++) {
                slen = strlen(argv[i]) + 1;
                MPI_Bcast(&slen, 1, MPI_INT, recvbuf.rank, MPI_COMM_WORLD);
                MPI_Bcast(argv[i], slen, MPI_CHAR, recvbuf.rank, MPI_COMM_WORLD);
                cmd.at(i) = argv[i];
            }
        } else {
            std::vector<char> buff;
            for (int i=0; i<recvbuf.argc; i++) {
                MPI_Bcast(&slen, 1, MPI_INT, recvbuf.rank, MPI_COMM_WORLD);
                if (buff.size() < static_cast<std::size_t>(slen)) {
                    buff.resize(slen);
                }
                MPI_Bcast(&buff[0], slen, MPI_CHAR, recvbuf.rank, MPI_COMM_WORLD);
                cmd.at(i) = &buff[0];
            }
        }
    }
}





namespace {
struct st_process_info {
    pid_t id;
    char name[MPI_MAX_PROCESSOR_NAME+1];
    st_process_info() {
        memset(name, 0, sizeof(name));
    }
};
}
/****************************************************************************//**
 * \param[in] comm The comunicater used where to fetch the hostname. All
 *    Processes in the communication must call the function together.
 * \param[out] res Returns a vector with the process id and
 *    \c MPI_Get_processor_name for each process in \a comm.
 * \param[in] root The rank of the process which recieves the gathered hostnames.
 *    If \a root is less then 0, each process gets the hostnames (Allgather).
 *    The root must be equal on all processes in comm.
 *******************************************************************************/
int tom::mpi::gather_hostname(MPI_Comm comm, std::vector<std::pair<pid_t, std::string> > &res, int root) {


    int rank;
    int numprocs;
    int i;
    int resultlen;
    int res_status;

    MPI_Comm_rank(comm, &rank);
    if (1 || root < 0 || root==rank) {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &numprocs);
    } else {
        rank = 0;
        numprocs = 1;
    }


    std::vector<st_process_info> info(numprocs>2?numprocs:2);
    info[rank].id = getpid();
    MPI_Get_processor_name(info[rank].name, &resultlen);

    MPI_Datatype _type[3] = { tom::mpi::get_MPI_Type<pid_t>(), MPI_CHAR, MPI_UB };

    int blocklen[3] = { 1, MPI_MAX_PROCESSOR_NAME+1, 1 };
    MPI_Aint disp[3];
    MPI_Aint base;
    MPI_Address(&info[0],            &base);
    MPI_Address(&info[0].id,         &disp[0]);
    MPI_Address(&info[0].name[0],    &disp[1]);
    MPI_Address(&info[1].id,         &disp[2]);
    for (i=0; i<3; i++) { disp[i] -= base; }

    MPI_Datatype type;

    MPI_Type_struct(3, blocklen, disp, _type, &type);
    MPI_Type_commit(&type);

    try {
        if (root < 0) {
            res_status = MPI_Allgather(&info[rank], 1, type, &info[0], 1, type, comm);
        } else {
            res_status = MPI_Gather(&info[rank], 1, type, &info[0], 1, type, root, comm);
        }
    } catch (...) {
        MPI_Type_free(&type);
        throw;
    }
    MPI_Type_free(&type);

    if (res_status == MPI_SUCCESS && (root<0||root==rank)) {
        res.resize(numprocs);
        std::vector<st_process_info               >::iterator sit = info.begin();
        std::vector<std::pair<pid_t, std::string> >::iterator dit = res.begin();
        for (; dit!=res.end(); sit++, dit++) {
            sit->name[MPI_MAX_PROCESSOR_NAME] = 0;
            dit->first = sit->id;
            dit->second = sit->name;
        }
    } else {
        res.resize(0);
    }

    return res_status;
}




/****************************************************************************//**
 * \brief Loads the fftw wisdom from file.
 *
 * Only one process per MPI_Get_processor_name reads acctually the file.
 * The other ones receive it from that process. This is a collective operation
 * for all processes in comm.
 *******************************************************************************/
int tom::mpi::load_fftw_wisdom(const std::string &fftw_wisdom_dir, MPI_Comm comm_world) {

    int res = 0;

    if (fftw_wisdom_dir.empty()) {
        res = tom::fftw::set_wisdom_name(NULL, 1);
    } else {

        const std::size_t l = MPI_MAX_PROCESSOR_NAME+1;

        int my_rank_world;
        MPI_Comm_rank(comm_world, &my_rank_world);
        int numprocs_world;
        MPI_Comm_size(comm_world, &numprocs_world);
        int i;

        std::string fftw_wisdom_filename;
        std::vector<char> all_proc_name(numprocs_world*l);


        char proc_name[l];
        {
            // The the processor name.
            int resultlen;
            MPI_Get_processor_name(proc_name, &resultlen);

            // Create the fftw_wisdom_filename.
            std::ostringstream ss;
            ss << fftw_wisdom_dir;
            if (!fftw_wisdom_dir.empty()) {
                if (fftw_wisdom_dir[fftw_wisdom_dir.size()-1] != '/') {
                    ss << '/';
                }
            }
            ss << "fftw_wisdom_" << proc_name << ".txt";
            fftw_wisdom_filename = ss.str();

            // Send the procname to the root.
            MPI_Gather(proc_name, l, MPI_CHAR, &all_proc_name[0], l, MPI_CHAR, 0, comm_world);
        }


        std::vector<int> color(numprocs_world);
        if (my_rank_world == 0) {
            // Group the processes by their processesor name.
            std::map<std::string, std::vector<int> > pn_grouped;
            for (i=0; i<numprocs_world; i++) {
                pn_grouped[&all_proc_name[i*l]].push_back(i);
            }
            i = 0;
            for (std::map<std::string, std::vector<int> >::const_iterator it=pn_grouped.begin(); it!=pn_grouped.end(); it++, i++) {
                for (std::vector<int>::const_iterator itv=it->second.begin(); itv!=it->second.end(); itv++) {
                    color[*itv] = i;
                }
            }
        }
        MPI_Scatter(&color[0], 1, MPI_INT, &i, 1, MPI_INT, 0, comm_world);

        MPI_Comm comm;
        int my_rank_comm;
        MPI_Comm_split(comm_world, i, my_rank_world, &comm);
        MPI_Comm_rank(comm, &my_rank_comm);

        int wisdom_length;
        std::vector<char> wisdom;
        if (my_rank_comm==0) {
            res = tom::fftw::set_wisdom_name(fftw_wisdom_filename.c_str(), 1);
            if (res) {
                char *cwisdom = fftw_export_wisdom_to_string();
                if (cwisdom) {
                    wisdom.resize(strlen(cwisdom)+1);
                    strcpy(&wisdom[0], cwisdom);
                    fftw_free(cwisdom);
                }
            }
            wisdom_length = wisdom.size();
            MPI_Bcast(&wisdom_length, 1, MPI_INT, 0, comm);
            if (wisdom_length) {
                MPI_Bcast(&wisdom[0], wisdom_length, MPI_CHAR, 0, comm);
            }
        } else {
            res = tom::fftw::set_wisdom_name(fftw_wisdom_filename.c_str(), 0);
            MPI_Bcast(&wisdom_length, 1, MPI_INT, 0, comm);
            if (wisdom_length) {
                wisdom.resize(wisdom_length);
                MPI_Bcast(&wisdom[0], wisdom_length, MPI_CHAR, 0, comm);
                res = fftw_import_wisdom_from_string(&wisdom[0]);
            }
        }
        MPI_Comm_free(&comm);

        //std::cerr << my_rank_world << ":" << proc_name << ": [" << fftw_export_wisdom_to_string() << "]" << std::endl;
    }



    return res;
}



/****************************************************************************//**
 * \brief saves the fftw wisdom back to file.
 *
 * Only one process per MPI_Get_processor_name acctually writes the file.
 * The other ones receive it from that process. This is a collective operation
 * for all processes in comm.
 *
 * The filename is as loaded before by load_fftw_wisdom.
 *******************************************************************************/
int tom::mpi::save_fftw_wisdom(MPI_Comm comm_world) {

    int res = 0;

    const std::size_t l = MPI_MAX_PROCESSOR_NAME+1;
    int my_rank_world;
    int numprocs_world;
    MPI_Comm_rank(comm_world, &my_rank_world);
    MPI_Comm_size(comm_world, &numprocs_world);

    int i;
    std::string fftw_wisdom_filename = tom::fftw::get_wisdom_name();
    std::vector<char> all_proc_name(numprocs_world*l);


    char proc_name[l];
    {
        // The the processor name.
        int resultlen;
        MPI_Get_processor_name(proc_name, &resultlen);

        // Send the procname to the root.
        MPI_Gather(proc_name, l, MPI_CHAR, &all_proc_name[0], l, MPI_CHAR, 0, comm_world);
    }


    std::vector<int> color(numprocs_world);
    if (my_rank_world == 0) {
        // Group the processes by their processesor name.
        std::map<std::string, std::vector<int> > pn_grouped;
        for (i=0; i<numprocs_world; i++) {
            pn_grouped[&all_proc_name[i*l]].push_back(i);
        }
        i = 0;
        for (std::map<std::string, std::vector<int> >::const_iterator it=pn_grouped.begin(); it!=pn_grouped.end(); it++, i++) {
            for (std::vector<int>::const_iterator itv=it->second.begin(); itv!=it->second.end(); itv++) {
                color[*itv] = i;
            }
        }
    }
    MPI_Scatter(&color[0], 1, MPI_INT, &i, 1, MPI_INT, 0, comm_world);

    MPI_Comm comm;
    MPI_Comm_split(comm_world, i, my_rank_world, &comm);
    int my_rank_comm;
    int numprocs_comm;
    MPI_Comm_rank(comm, &my_rank_comm);
    MPI_Comm_size(comm, &numprocs_comm);


    int wisdom_length;
    std::vector<char> wisdom;
    if (my_rank_comm==0) {
        MPI_Status status;
        int src;
        for (i=1; i<numprocs_comm; i++) {
            MPI_Probe(MPI_ANY_SOURCE, 0, comm, &status);
            src = status.MPI_SOURCE;
            MPI_Get_count(&status, MPI_CHAR, &wisdom_length);
            wisdom.resize(wisdom_length+1);
            MPI_Recv(&wisdom[0], wisdom_length, MPI_CHAR, src, 0, comm, &status);
            fftw_import_wisdom_from_string(&wisdom[0]);
        }
        res = tom::fftw::set_wisdom_name(fftw_wisdom_filename.c_str(), 0);
        tom::fftw::save_wisdom();
    } else {
        char *cwisdom = fftw_export_wisdom_to_string();
        if (cwisdom) {
            wisdom.resize(strlen(cwisdom)+1);
            strcpy(&wisdom[0], cwisdom);
            fftw_free(cwisdom);
        }
        MPI_Send(&wisdom[0], wisdom.size(), MPI_CHAR, 0, 0, comm);
    }

    MPI_Comm_free(&comm);


    return res;
}




/****************************************************************************//**
 * \brief Sends a std::string to all processes in comm using a broadcast.
 *
 * This is a collective operation.
 *******************************************************************************/
int tom::mpi::bcast_string(int root, MPI_Comm comm, std::string &s) {

    int res;

    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    int i;
    __MAKE_STATIC__ std::vector<char> cbuff;
    if (my_rank == root) {
        i = s.size()+1;
        if ((res = MPI_Bcast(&i, 1, MPI_INT, root, comm)) != MPI_SUCCESS) {
            return res;
        }
        if (cbuff.size() < static_cast<std::size_t>(i)) { cbuff.resize(i); }
        strcpy(&cbuff[0], s.c_str());
        res = MPI_Bcast(&s[0], i, MPI_CHAR, root, comm);
    } else {
        if ((res = MPI_Bcast(&i, 1, MPI_INT, root, comm)) != MPI_SUCCESS) {
            return res;
        }
        if (cbuff.size() < static_cast<std::size_t>(i)) { cbuff.resize(i); }
        if ((res = MPI_Bcast(&cbuff[0], i, MPI_CHAR, root, comm)) == MPI_SUCCESS) {
            s = &cbuff[0];
        }
    }
    return res;
}





/****************************************************************************//**
 * \brief Broadcasts a tom::Volume.
 *
 * This is a collective operation.
 *******************************************************************************/
template<typename T>
int tom::mpi::bcast_volume(int root, MPI_Comm comm, tom::Volume<T> *&v) {

    int res;

    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    unsigned long dims[3] = { 0,0,0 };
    std::auto_ptr<tom::Volume<T> > av;
    tom::Volume<T> *pv = v;

    if (my_rank == root) {
        if (v) {
            dims[0] = v->getSizeX();
            dims[1] = v->getSizeY();
            dims[2] = v->getSizeZ();
        }
        if ((res=MPI_Bcast(dims, 3, MPI_UNSIGNED_LONG, root, comm)) != MPI_SUCCESS) {
            return res;
        }
        if (dims[0]) {
            if (!v->isContiguous()) {
                av.reset(pv = new tom::Volume<T>(*v));
            }
            res = MPI_Bcast(&pv->get(), dims[0]*dims[1]*dims[2], tom::mpi::get_MPI_Type<T>(), root, comm);
        }
    } else {
        if ((res=MPI_Bcast(dims, 3, MPI_UNSIGNED_LONG, root, comm)) != MPI_SUCCESS) {
            return res;
        }
        if (dims[0] == 0) {
            v = NULL;
            return res;
        }
        av.reset(pv = new tom::Volume<T>(dims[0], dims[1], dims[2], NULL,NULL));

        if ((res = MPI_Bcast(&pv->get(), dims[0]*dims[1]*dims[2], tom::mpi::get_MPI_Type<T>(), root, comm)) == MPI_SUCCESS) {
            v = av.release();
        }
    }

    return res;
}




/****************************************************************************//**
 * \brief Sends a vector of strings to a process.
 *
 * \param[in] dst Rank of the reciever.
 * \param[in] tag Tag of the send.
 * \param[in] comm The intra-communicater where sender and reciever are.
 * \param[in] v A std::vector containing strings.
 * \return The error code of the called mpi-functions.
 *
 * This is the sender side part to tom::mpi::recv_vstring.
 *******************************************************************************/
int tom::mpi::send_vstring(int dst, int tag, MPI_Comm comm, const std::vector<std::string> &v) {

    int res;

    const std::size_t nsize = v.size();

    if (nsize < 1) {
        unsigned long bb = 0;
        return MPI_Send(&bb, 1, MPI_UNSIGNED_LONG, dst, tag, comm);
    }

    std::size_t i;
    std::vector<unsigned long> ulbuff(nsize);
    std::size_t total_length = 0;
    for (i=0; i<nsize; i++) {
        ulbuff[i] = (total_length += v[i].size() + 1);
    }

    if ((res=MPI_Send(&ulbuff[0], nsize, MPI_UNSIGNED_LONG, dst, tag, comm)) == MPI_SUCCESS) {
        __MAKE_STATIC__ std::vector<char> cbuff;
        cbuff.resize(total_length);

        std::size_t last_start_idx = 0;
        for (i=0; i<nsize; i++) {
            strcpy(&cbuff[last_start_idx], v[i].c_str());
            last_start_idx = ulbuff[i];
        }
        res = MPI_Ssend(&cbuff[0], total_length, MPI_CHAR, dst, tag, comm);
    }

    return res;
}



/****************************************************************************//**
 * \brief Recieves a vector of strings to a process.
 *
 * \param[in] src Rank of the sender.
 * \param[in] tag Tag of the send.
 * \param[in] comm The intra-communicater where sender and reciever are.
 * \param[out] v A std::vector containing strings.
 * \return The error code of the called mpi-functions.
 *
 * This is the sender side part to tom::mpi::recv_vstring.
 *******************************************************************************/
int tom::mpi::recv_vstring(int src, int tag, MPI_Comm comm, std::vector<std::string> &v) {

    int res;
    int nsize;
    MPI_Status status;

    MPI_Probe(src, tag, comm, &status);
    MPI_Get_count(&status, MPI_UNSIGNED_LONG, &nsize);

    std::vector<unsigned long> ulbuff(nsize);
    if ((res=MPI_Recv(&ulbuff[0], nsize, MPI_UNSIGNED_LONG, src, tag, comm, &status)) == MPI_SUCCESS) {
        if (nsize==1 && ulbuff[0]==0) {
            v.resize(0);
        } else {
            const std::size_t total_length = ulbuff.back();

            __MAKE_STATIC__ std::vector<char> cbuff;
            cbuff.resize(total_length);

            if ((res=MPI_Recv(&cbuff[0], total_length, MPI_CHAR, src, tag, comm, &status)) != MPI_SUCCESS) {
                return res;
            }
            std::size_t last_start_idx = 0;
            v.resize(nsize);
            for (int i=0; i<nsize; i++) {
                v[i] = &cbuff[last_start_idx];
                last_start_idx = ulbuff[i];
            }
        }
    }
    return res;
}





/****************************************************************************//**
 * \brief Sends wedges to a process.
 *
 * \param[in] dst Rank of the reciever.
 * \param[in] tag Used tag for the communication.
 * \param[in] comm Communicator.
 * \param[in] wedges A vector with pointers to tom::WedgeDescriptor.
 *    Elements set to null arrive as tom::NoWedge at the reciever.
 *    Elements of the vector can point to the same instance of a WedgeDescriptor.
 * \return The error status of the last MPI-command that failed (or MPI_SUCCESS).
 *
 * This is the counterpart to recv_wedges. The differene is that here a
 * vector of tom::WedgeDescriptor is given and the reciever gets the instanciated
 * wedges.\n
 * As these are template functions, the counterpart must have the same template type.
 *******************************************************************************/
template<typename T>
int tom::mpi::send_wedges(int dst, int tag, MPI_Comm comm, const std::vector<tom::WedgeDescriptor<T> *> &wedges) {

    int res;

    std::map<int, std::vector<int> > idxlist;
    std::map<int, std::vector<int> >::iterator it;

    std::size_t nsize = wedges.size();
    std::size_t i;

    if (nsize < 1) {
        int ii = 0;
        return MPI_Send(&ii, 1, MPI_INT, dst, tag, comm);
    }

    {
        // Group the indices of the wedges by their type.
        const tom::WedgeDescriptor<T> *pwd;
        for (i=0; i<nsize; i++) {
            if ((pwd = wedges[i])) {
                idxlist[pwd->getTypeID()].push_back(i);
            } else {
                idxlist[tom::WedgeDescriptor_NoWedge<T>::getClassTypeID()].push_back(i);
            }
        }
    }

    {
        // Pack the grouped list of indices to send them to the reciever.
        // The packing works like this:
        // | ntypes | id1 | id1_n | id1_1 | id1_2 ... | id2 | id2_n | id2_1 ...
        // The first element is the total number of different wedge types (ntypes).
        // Then for each of the ntypes its id (id1, id2, ...), the number of
        // indices with the same id (id1_n, id2_n, ...) and the actual indices
        // itself.
        int total_size = 1+2*idxlist.size();
        std::vector<int> ii(total_size+nsize);
        ii[0] = idxlist.size();
        i = 0;
        int j, k;
        for (it=idxlist.begin(); it!=idxlist.end(); it++) {
            std::vector<int> &v = (*it).second;
            ii[++i] = (*it).first;
            ii[++i] = total_size;
            k = v.size();
            for (j=0; j<k; j++) {
                ii[j+total_size] = v[j];
            }
            total_size += k;
        }

        res = MPI_Send(&ii[0], ii.size(), MPI_INT, dst, tag, comm);
    }

    if (res == MPI_SUCCESS) {
        std::size_t i, vsize;
        for (it=idxlist.begin(); it != idxlist.end(); it++) {
            const int &ctypeid = (*it).first;
            const std::vector<int> &v = (*it).second;
            vsize = v.size();

            if (ctypeid == tom::WedgeDescriptor_EMWedge<T>::getClassTypeID()) {
                std::vector<std::size_t> vunsigned(vsize);
                std::vector<std::string> vstring(vsize);
                for (i=0; i<vsize; i++) {
                    assert(typeid(*wedges.at(v.at(i))) == typeid(tom::WedgeDescriptor_EMWedge<T>));
                    const tom::WedgeDescriptor_EMWedge<T> &pp = dynamic_cast<const tom::WedgeDescriptor_EMWedge<T> &>(*wedges[v[i]]);
                    vstring[i] = pp.getFilename();
                    vunsigned[i] = pp.getBinning();
                }
                if ((res=tom::mpi::send_vstring(dst, tag, comm, vstring)             ) != MPI_SUCCESS) {
                    return res;
                }
                if ((res=MPI_Send(&vunsigned[0], vsize, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS) {
                    return res;
                }
            } else if (ctypeid == tom::WedgeDescriptor_SimpleWedge<T>::getClassTypeID()) {
                std::vector<double> vdouble(2*vsize);
                for (i=0; i<vsize; i++) {
                    assert(typeid(*wedges.at(v.at(i))) == typeid(const tom::WedgeDescriptor_SimpleWedge<T>));
                    const tom::WedgeDescriptor_SimpleWedge<T> &pp = dynamic_cast<const tom::WedgeDescriptor_SimpleWedge<T> &>(*wedges[v[i]]);
                    vdouble[i*2+0] =  pp.getAngle();
                    vdouble[i*2+1] =  pp.getCutoff_radius();
                }
                if ((res=MPI_Ssend(&vdouble[0], 2*vsize, MPI_DOUBLE, dst, tag, comm)) != MPI_SUCCESS) {
                    return res;
                }
            } else {
                assert(ctypeid == tom::WedgeDescriptor_NoWedge<T>::getClassTypeID() );
                // Here should only be NoWedge, or there is an not yet known subclass of WedgeDescriptor.
            }
        }
    }

    return res;
}








/****************************************************************************//**
 * \brief Recieves wedges from a process.
 *
 * \param[in] src Rank of the reciever.
 * \param[in] tag Used tag for the communication.
 * \param[in] comm Communicator.
 * \param[out] wedges A vector with pointers to tom::Wedge.
 *    Elements set to null at the sender arrive as tom::NoWedge at the reciever.
 * \return The error status of the last MPI-command that failed (or MPI_SUCCESS).
 *
 * This is the counterpart to send_wedges.
 * If one or more wedges are basically the same, also the content of \a wedges
 * point to the same instance of tom::Wedge. This is desired, as the implementation
 * of wedges caches the wedge under its last rotation. So using all the wedges under
 * the same rotation saves the work to actually rotate the volumes.\n
 * As these are template functions, the counterpart must have the same template type.
 *******************************************************************************/
template<typename T>
int tom::mpi::recv_wedges(int src, int tag, MPI_Comm comm,
                    std::vector<boost::shared_ptr<tom::Wedge<T> > > &wedges) {

    int res;

    MPI_Status status;

    MPI_Probe(src, tag, comm, &status);

    int count;
    MPI_Get_count(&status, MPI_INT, &count);
    std::vector<int> ii(count);
    if ((res=MPI_Recv(&ii[0], count, MPI_INT, src, tag, comm, &status)) != MPI_SUCCESS) {
        return res;
    }
    if (count==0 || ii[0]==0) {
        wedges.resize(0);
    } else {
        int i, j;
        std::map<int, std::vector<int> > idxlist;
        std::map<int, std::vector<int> >::iterator it;

        const int nitems = ii.at(0);
        {
            int end_idx = count;
            int start_idx;
            for (i=nitems-1; i>=0; i--) {
                std::vector<int> &v = idxlist[ii.at(i*2+1)];
                start_idx = ii.at(i*2+2);
                for (j=start_idx; j<end_idx; j++) {
                    v.push_back(ii.at(j));
                }
            }
        }

        std::vector<boost::shared_ptr<tom::Wedge<T> > > wedges_local(count - 1 - 2*nitems);

        {
            std::size_t i, vsize;
            int ctypeid;
            for (it=idxlist.begin(); it != idxlist.end(); it++) {
                ctypeid = (*it).first;
                const std::vector<int> &v = (*it).second;
                vsize = v.size();

                if (ctypeid == tom::WedgeDescriptor_EMWedge<T>::getClassTypeID()) {
                    std::vector<std::string> vstring;
                    if ((res=tom::mpi::recv_vstring(src, tag, comm, vstring)) != MPI_SUCCESS) {
                        return res;
                    }
                    assert(vstring.size() == vsize);

                    std::vector<std::size_t> vunsigned(vsize);
                    #ifndef NDEBUG
                    {
                        MPI_Probe(src, tag, comm, &status);
                        int count;
                        MPI_Get_count(&status, tom::mpi::get_MPI_Type<std::size_t>(), &count);
                        assert((std::size_t)count == vsize);
                    }
                    #endif
                    if ((res=MPI_Recv(&vunsigned[0], vsize, tom::mpi::get_MPI_Type<std::size_t>(), src, tag, comm, &status)) != MPI_SUCCESS) {
                        return res;
                    }

                    // Create the map pw to create only one simple-wedge object per parameter combination.
                    std::map<std::pair<std::size_t, std::string>, boost::shared_ptr<tom::Wedge<T> > > pw;
                    typename std::map<std::pair<std::size_t, std::string>, boost::shared_ptr<tom::Wedge<T> > >::const_iterator pw_it;
                    for (i=0; i<vsize; i++) {
                        const std::string &filename = vstring[i];
                        const std::size_t &binning = vunsigned[i];
                        std::pair<std::size_t, std::string> iwedge(binning, filename);
                        pw_it = pw.find(iwedge);
                        if (pw_it == pw.end()) {
                            pw_it = pw.insert(std::make_pair(iwedge, boost::shared_ptr<tom::Wedge<T> >(new tom::VolumeEMWedge<T>(filename, binning)))).first;
                        }
                        wedges_local.at(v[i]) = (*pw_it).second;
                    }
                } else if (ctypeid == tom::WedgeDescriptor_SimpleWedge<T>::getClassTypeID()) {
                    std::vector<double> dbuff(2*vsize);
                    if ((res=MPI_Recv(&dbuff[0], 2*vsize, MPI_DOUBLE, src, tag, comm, &status)) != MPI_SUCCESS) {
                        return res;
                    }

                    // Create the map pw to create only one simple-wedge object per parameter combination.
                    std::map<std::pair<double, double>, boost::shared_ptr<tom::Wedge<T> > > pw;
                    typename std::map<std::pair<double, double>, boost::shared_ptr<tom::Wedge<T> > >::iterator pw_it;
                    std::pair<double, double> key;
                    for (i=0; i<vsize; i++) {
                        key.first = dbuff[i*2+0];
                        key.second = dbuff[i*2+1];
                        pw_it = pw.find(key);
                        if (pw_it == pw.end()) {
                            // Allocate a new simplewedge.
                            boost::shared_ptr<tom::Wedge<T> >  cntptr_w(new tom::SimpleWedge<T>(key.first, key.second));
                            std::pair<std::pair<double,double>, boost::shared_ptr<tom::Wedge<T> > > new_val(key, cntptr_w);
                            pw_it = pw.insert(new_val).first;
                        }
                        wedges_local.at(v[i]) = (*pw_it).second;
                    }
                } else {
                    assert(ctypeid == tom::WedgeDescriptor_NoWedge<T>::getClassTypeID());
                    boost::shared_ptr<tom::Wedge<T> > pw(new tom::NoWedge<T>());
                    for (i=0; i<vsize; i++) {
                        wedges_local.at(v[i]) = pw;
                    }
                }
            }
        }
        wedges.swap(wedges_local);
    }
    return res;
}







/****************************************************************************//**
 * \brief Send a tom::Volume to a process.
 *
 * \param[in] dst Rank of the reciever.
 * \param[in] tag Tag of the direct communication.
 * \param[in] comm Communicator
 * \param[in] v Volume to be transmitted.
 * \returns Error-code from the last MPI-Command.
 *
 * This is the counterpart to tom::mpi::recv_volume.
 * As these are template functions, the counterpart must have the same template type.
 *******************************************************************************/
template<typename T>
int tom::mpi::send_volume(int dst, int tag, MPI_Comm comm, const tom::Volume<T> &v) {

    int res;

    long unsigned dims[3] = { v.getSizeX(), v.getSizeY(), v.getSizeZ() };
    if ((res=MPI_Send(dims, 3, MPI_UNSIGNED_LONG, dst, tag, comm)) != MPI_SUCCESS) {
        return res;
    }

    const std::size_t stridex = v.getStrideX();
    const std::size_t stridey = v.getStrideY();
    const std::size_t stridez = v.getStrideZ();

    MPI_Datatype vtype;

    std::auto_ptr<tom::Volume<T> > vtmp;
    T *pv = const_cast<T *>(&v.get());

    if (v.isContiguous()) {
        MPI_Type_contiguous(dims[0]*dims[1]*dims[2], get_MPI_Type<T>(), &vtype);
    } else {
        if (stridex%sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
            vtmp.reset(new tom::Volume<T>(v));
            pv = &vtmp->get();
        } else if (stridex>sizeof(T) && stridey==dims[0]*stridex && stridez==dims[1]*stridey) {
            MPI_Type_vector(dims[0]*dims[1]*dims[2], 1, stridex/sizeof(T), get_MPI_Type<T>(), &vtype);
        } else if (stridex==sizeof(T) && stridey>dims[0]*stridex && stridez==dims[1]*stridey) {
            MPI_Type_vector(dims[1]*dims[2], dims[0], stridey/sizeof(T), get_MPI_Type<T>(), &vtype);
        } else if (stridex==sizeof(T) && stridey==dims[0]*stridex && stridez>dims[1]*stridey) {
            MPI_Type_vector(dims[2], dims[0]*dims[1], stridez/sizeof(T), get_MPI_Type<T>(), &vtype);
        } else {
            vtmp.reset(new tom::Volume<T>(v));
            pv = &vtmp->get();
        }
    }
    MPI_Type_commit(&vtype);
    try {
        res = MPI_Ssend(pv, 1, vtype, dst, tag, comm);
    } catch (...) {
        MPI_Type_free(&vtype);
        throw;
    }
    MPI_Type_free(&vtype);

    return res;
}


/****************************************************************************//**
 * \brief Recieve a tom::Volume from a process.
 *
 * \param[in] src Rank of the sender.
 * \param[in] tag Tag of the direct communication.
 * \param[in] comm Communicator
 * \param[out] v Pointer to the recieved volume.
 * \returns Error-code from the last MPI-Command. In case of error, \a v remains unchanged.
 *
 * This is the counterpart to tom::mpi::send_volume.
 * As these are template functions, the counterpart must have the same template type.
 *******************************************************************************/
template<typename T>
int tom::mpi::recv_volume(int src, int tag, MPI_Comm comm, std::auto_ptr<tom::Volume<T> > &v) {

    int res;

    MPI_Status status;
    long unsigned dims[3];

    if ((res = MPI_Recv(dims, 3, MPI_UNSIGNED_LONG, src, tag, comm, &status)) == MPI_SUCCESS) {
        std::auto_ptr<tom::Volume<T> > vtmp(new tom::Volume<T>(dims[0], dims[1], dims[2], NULL,NULL));

        if ((res=MPI_Recv(&vtmp->get(), dims[0]*dims[1]*dims[2], get_MPI_Type<T>(), src, tag, comm, &status)) == MPI_SUCCESS) {
            v = vtmp;
        }
    }
    return res;
}




/****************************************************************************//**
 * \brief Request a job from the server.
 *
 * \param[in] dest Rank of the destination.
 * \param[in] tag Tag of the direct communication.
 *
 *
 *******************************************************************************/
int tom::mpi::os3_send_job_request(int dst, int tag, MPI_Comm comm,int senderRank){
	int res;

	if (( res=MPI_Send(&senderRank, 1, MPI_INT, 0, tag, comm)) != MPI_SUCCESS) {
        return res;
    }
	return MPI_SUCCESS;
}


/****************************************************************************//**
 * \brief Recieves a job request.
 *
 * \param[in] comm MPI_Comm
 * \param[in] senderRank Will contain request id if recieve has been successful.
 *
 *
 *******************************************************************************/
int tom::mpi::os3_recieve_job_request(MPI_Comm comm,int& senderRank){
	int res;

	MPI_Status status;

	try{
		res = MPI_Recv(&senderRank, 1, MPI_INT, senderRank, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		throw;
	}
	return MPI_SUCCESS;
}

/****************************************************************************//**
 * \brief Send a job to a worker
 *
 * \param[in] dest Rank of the destination.
 * \param[in] tag Tag of the direct communication.
 * \param[in] job The job to be done
 *
 *
 *******************************************************************************/
int tom::mpi::send_os3_job(int dst, int tag, MPI_Comm comm, const tom::os3_job& job) {
	int res;



	std::size_t stringSize;

	if((res=MPI_Send(const_cast<int*>(&job.jobType), 1,MPI_INT, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	stringSize = job.volumePath.size()+1;
	//send size of the string first, and then the string afterwards
	if ( (res=MPI_Send(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
		return res;
	if((res=MPI_Send(const_cast<char *>(job.volumePath.c_str()), stringSize, MPI_CHAR, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	stringSize = job.patternPath.size()+1;
	if ( (res=MPI_Send(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
		return res;
	if((res=MPI_Send(const_cast<char *>(job.patternPath.c_str()), stringSize, MPI_CHAR, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	stringSize = job.maskPath.size()+1;
	if ( (res=MPI_Send(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
		return res;
	if((res=MPI_Send(const_cast<char *>(job.maskPath.c_str()), stringSize, MPI_CHAR, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	stringSize = job.psfPath.size()+1;
	if ( (res=MPI_Send(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
		return res;
	if((res=MPI_Send(const_cast<char *>(job.psfPath.c_str()), stringSize, MPI_CHAR, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	stringSize = job.resultPath.size()+1;
	if ( (res=MPI_Send(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
		return res;
	if((res=MPI_Send(const_cast<char *>(job.resultPath.c_str()), stringSize, MPI_CHAR, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<uint*>(job.subregion), 6, MPI_UNSIGNED, dst, tag, comm)) != MPI_SUCCESS)
        return res;
	if((res=MPI_Send(const_cast<uint*>(job.writeBackSubregion), 6, MPI_UNSIGNED, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<uint*>(job.binning), 3, MPI_UNSIGNED, dst, tag, comm)) != MPI_SUCCESS)
        return res;
	if((res=MPI_Send(const_cast<uint*>(job.sampling), 3, MPI_UNSIGNED, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<std::size_t*>(&job.id), 1,tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<std::size_t*>(&job.numberJobs), 1,tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<double*>(job.weight), 1,MPI_DOUBLE, dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<std::size_t*>(&job.numberParticles), 1,tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
        return res;

    #ifdef __USE_MPI_BOOL
	if((res=MPI_Send(const_cast<bool*>(&job.filterSave), 1,tom::mpi::get_MPI_Type<bool>(), dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<bool*>(&job.listSave), 1,tom::mpi::get_MPI_Type<bool>(), dst, tag, comm)) != MPI_SUCCESS)
        return res;

	if((res=MPI_Send(const_cast<bool*>(&job.stackSave), 1,tom::mpi::get_MPI_Type<bool>(), dst, tag, comm)) != MPI_SUCCESS)
        return res;
    #else
    {
        char c;
        c = job.filterSave == true;
        if((res=MPI_Send(&c, 1, tom::mpi::get_MPI_Type<char>(), dst, tag, comm)) != MPI_SUCCESS) {
            return res;
        }

        c = job.listSave == true;
        if((res=MPI_Send(&c, 1,tom::mpi::get_MPI_Type<char>(), dst, tag, comm)) != MPI_SUCCESS)
            return res;

        c = job.stackSave == true;
        if((res=MPI_Send(&c, 1,tom::mpi::get_MPI_Type<char>(), dst, tag, comm)) != MPI_SUCCESS)
            return res;
    }
    #endif

	std::size_t sizeAngleList = job.angleList.size();
	if ( (res=MPI_Send(&sizeAngleList, 1, tom::mpi::get_MPI_Type<std::size_t>(), dst, tag, comm)) != MPI_SUCCESS)
		return res;

	for(std::size_t i=0;i<sizeAngleList;i++){
		double angles[3] = {0,0,0};
		angles[0] = job.angleList[i].phi;
		angles[1] = job.angleList[i].psi;
		angles[2] = job.angleList[i].theta;
		//send each element of list
		if((res=MPI_Send(const_cast<double*>(angles), 3,MPI_DOUBLE, dst, tag, comm)) != MPI_SUCCESS)
        	return res;
	}

	return res;

}





/****************************************************************************//**
 * \brief Receive a os3_job
 *
 * \param[in] src .
 * \param[in] tag Tag of the direct communication.
 * \param[in] comm
 * \param[in] job The job to be done
 *******************************************************************************/
int tom::mpi::recv_os3_job(int src, int tag, MPI_Comm comm, tom::os3_job& job){

    MPI_Status status;
    int res = MPI_SUCCESS;
	try{
		res = MPI_Recv(&job.jobType, 1, MPI_INT, src, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving job.type!\n";
		return res;
	}
	std::size_t stringSize;
	try{
		res = MPI_Recv(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);
		job.volumePath = std::string(stringSize, '0');
		res = MPI_Recv(const_cast<char *>(job.volumePath.c_str()), stringSize, MPI_CHAR, src, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving job.volumePath!\n";
		return res;
	}

	try{
		res = MPI_Recv(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);
		job.patternPath= std::string(stringSize, '0');
		res = MPI_Recv(const_cast<char *>(job.patternPath.c_str()), stringSize, MPI_CHAR, src, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving job.patternPath!\n";
		return res;
	}
	try{
		res = MPI_Recv(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);
		job.maskPath= std::string(stringSize, '0');
		res = MPI_Recv(const_cast<char *>(job.maskPath.c_str()), stringSize, MPI_CHAR, src, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving job.maskPath!\n";
		return res;
	}

	try{
		res = MPI_Recv(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);
		job.psfPath= std::string(stringSize, '0');
		res = MPI_Recv(const_cast<char *>(job.psfPath.c_str()), stringSize, MPI_CHAR, src, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving job.psfPath!\n";
		return res;
	}

	try{
		res = MPI_Recv(&stringSize, 1, tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);
		job.resultPath= std::string(stringSize, '0');
		res = MPI_Recv(const_cast<char *>(job.resultPath.c_str()), stringSize, MPI_CHAR, src, MPI_ANY_TAG, comm, &status);
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving job.resultPath!\n";
		return res;
	}

	try{
		res = MPI_Recv(job.subregion, 6, MPI_UNSIGNED, src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(job.writeBackSubregion, 6, MPI_UNSIGNED, src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(job.binning, 3, MPI_UNSIGNED, src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(job.sampling, 3, MPI_UNSIGNED, src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(&job.id, 1,tom::mpi::get_MPI_Type<size_t>(), src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(&job.numberJobs, 1,tom::mpi::get_MPI_Type<size_t>(), src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(job.weight, 3,MPI_DOUBLE, src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(&job.numberParticles, 1,tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);
        #ifdef __USE_MPI_BOOL
		res = MPI_Recv(&job.filterSave, 1,tom::mpi::get_MPI_Type<bool>(), src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(&job.listSave, 1,tom::mpi::get_MPI_Type<bool>(), src, MPI_ANY_TAG, comm, &status);
		res = MPI_Recv(&job.stackSave, 1,tom::mpi::get_MPI_Type<bool>(), src, MPI_ANY_TAG, comm, &status);
        #else
        {
            char c;
            res = MPI_Recv(&c, 1,tom::mpi::get_MPI_Type<char>(), src, MPI_ANY_TAG, comm, &status);
            job.filterSave = c != 0;
            res = MPI_Recv(&c, 1,tom::mpi::get_MPI_Type<char>(), src, MPI_ANY_TAG, comm, &status);
            job.listSave = c != 0;
            res = MPI_Recv(&c, 1,tom::mpi::get_MPI_Type<char>(), src, MPI_ANY_TAG, comm, &status);
            job.stackSave = c != 0;
        }
        #endif
	}
	catch(...){
		std::cerr << "tom::mpi::recv_os3_job : error receiving numeric values (subregion, binning, ...) !\n";
		return res;
	}

	std::size_t sizeAngleList;
	res = MPI_Recv(&sizeAngleList, 1, tom::mpi::get_MPI_Type<std::size_t>(), src, MPI_ANY_TAG, comm, &status);

	tom::angleTupel tupel;
	for(std::size_t i=0;i<sizeAngleList;i++){
		double angles[3] = {0,0,0};
		try{
			res = MPI_Recv(angles, 3,MPI_DOUBLE, src, MPI_ANY_TAG, comm, &status);
		}
		catch(...){
			std::cerr << "tom::mpi::recv_os3_job : error receiving angle tuple!\n";
			return res;
		}
		tupel.phi = angles[0];
		tupel.psi = angles[1];
		tupel.theta = angles[2];

		job.angleList.push_back(tupel);
	}

	return MPI_SUCCESS;
}





// template instantiations.
template int tom::mpi::bcast_volume<char  >(int root, MPI_Comm comm, tom::Volume<char  > *&v);
template int tom::mpi::bcast_volume<int   >(int root, MPI_Comm comm, tom::Volume<int   > *&v);
template int tom::mpi::bcast_volume<float >(int root, MPI_Comm comm, tom::Volume<float > *&v);


template int tom::mpi::send_wedges<float >(int dst, int tag, MPI_Comm comm, const std::vector<tom::WedgeDescriptor<float > *> &wedges);
template int tom::mpi::send_wedges<double>(int dst, int tag, MPI_Comm comm, const std::vector<tom::WedgeDescriptor<double> *> &wedges);
template int tom::mpi::recv_wedges<float >(int src, int tag, MPI_Comm comm, std::vector<boost::shared_ptr<tom::Wedge<float > > > &wedges);
template int tom::mpi::recv_wedges<double>(int src, int tag, MPI_Comm comm, std::vector<boost::shared_ptr<tom::Wedge<double> > > &wedges);


template int tom::mpi::send_volume<float >(int dst, int tag, MPI_Comm comm, const tom::Volume<float > &v);
template int tom::mpi::send_volume<double>(int dst, int tag, MPI_Comm comm, const tom::Volume<double> &v);
template int tom::mpi::recv_volume<float >(int src, int tag, MPI_Comm comm, std::auto_ptr<tom::Volume<float > > &v);
template int tom::mpi::recv_volume<double>(int src, int tag, MPI_Comm comm, std::auto_ptr<tom::Volume<double> > &v);


template void tom::mpi::send_peaklist<float >(const std::vector<tom::cc_peak<float > > &peak_list, int dst, int msg_tag, MPI_Comm comm);
template void tom::mpi::send_peaklist<double>(const std::vector<tom::cc_peak<double> > &peak_list, int dst, int msg_tag, MPI_Comm comm);
template void tom::mpi::recv_peaklist<float >(std::vector<tom::cc_peak<float > > &peak_list, int src, int msg_tag, MPI_Comm comm);
template void tom::mpi::recv_peaklist<double>(std::vector<tom::cc_peak<double> > &peak_list, int src, int msg_tag, MPI_Comm comm);


template MPI_Datatype tom::mpi::create_cc_peak_type<float >();
template MPI_Datatype tom::mpi::create_cc_peak_type<double>();



