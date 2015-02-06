/****************************************************************************//**
 * \file exec.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    18.02.2007
 *******************************************************************************/
#ifndef ___INCLUDE_HELPER_EXEC_HPP__
#define ___INCLUDE_HELPER_EXEC_HPP__



#include <errno.h>
#include <vector>
#include <stdexcept>
#include <boost/shared_array.hpp>




namespace helper {



template<typename T> void exec(T e, std::vector<char> &stdout_, std::vector<char> &stderr_, int *status_);





/****************************************************************************//**
 * class for handling a pipe.
 *
 * Currently not thread save because of its use of strerror.
 *******************************************************************************/
class pipe {
public:
    pipe() {
        if (::pipe(fd) != 0) {
            const int errno_ = errno;
            throw std::runtime_error(std::string("error creating pipe: ") + strerror(errno_));
        }
        fd_open[0] = fd_open[1] = true;
    }
    ~pipe() {
        close();
    }
    void close() {
        close_fd(0);
        close_fd(1);
    }
    void close_nothrow() {
        close_nothrow_fd(0);
        close_nothrow_fd(1);
    }
    void close_fd(int i) {
        if (fd_open[i=i?1:0]) {
            fd_open[i] = false;
            if (::close(fd[i]) != 0) {
                const int errno_ = errno;
                throw std::runtime_error(std::string("error closing pipe: ") + strerror(errno_));
            }
        }
    }
    void close_nothrow_fd(int i) {
        if (fd_open[i = i?1:0]) {
            fd_open[i] = false;
            ::close(fd[i]);
        }
    }
    int operator[](int i) const {
        i = i ? 1 : 0;
        if (!fd_open[i]) { throw std::runtime_error("file descriptor already closed."); }
        return fd[i];
    }
    bool is_open_fd(int i) const {
        return fd_open[i ? 1 : 0];
    }
private:
    int fd[2];
    bool fd_open[2];
}; // class helper::pipe



/****************************************************************************//**
 *
 *******************************************************************************/
class c_execv {
public:
    c_execv(const std::string &file, const std::vector<std::string> &argv, bool use_path);
    inline int operator()() const {
        int (*e)(const char *file, char *const argv[]) = use_path ? &::execvp : &::execv;
        #if 0
        std::cerr << ">>>>>>> " << &file_[0] << std::endl;
        const char *const *p = &argv_p[0];
        int i=0;
        while (*p) {
            std::cerr << ">>>>>>> " << i << ": _" << *p << "_" << std::endl;
            p++;
            i++;
        }
        #endif
        errno = 0;
        e(&file_[0], &argv_p[0]);
        return errno;
    }
private:
    c_execv(const c_execv &c) { }
    c_execv &operator=(const c_execv &c) { return *this; }

    bool use_path;
    std::vector<char> file_;
    std::vector<boost::shared_array<char> > argv_v;
    std::vector<char *> argv_p;
}; // class helper::c_execvp


} // namespace helper






// inline functions.



/****************************************************************************//**
 *
 *******************************************************************************/
inline helper::c_execv::c_execv(const std::string &file, const std::vector<std::string> &argv, bool use_path)
    : use_path(use_path),
      file_(),
      argv_v(),
      argv_p() {
    if (file.empty()) {
        throw std::invalid_argument("the file parameter can not be empty.");
    }
    file_.reserve(file.size()+1);
    file_.assign(file.begin(), file.end());
    file_.push_back(0);

    std::vector<std::string> argv2;
    if (argv.empty()) {
        argv2.push_back(file);
    }
    const std::vector<std::string> &a = argv.empty() ? argv2 : argv;

    argv_v.reserve(a.size());
    argv_p.reserve(a.size()+1);
    for (std::vector<std::string>::const_iterator it=a.begin(); it!=a.end(); it++) {
        char *c = new char[it->size()+1];
        argv_v.push_back(boost::shared_array<char>(c));
        argv_p.push_back(c);
        std::strcpy(c, it->c_str());
    }
    argv_p.push_back(NULL);
}


/****************************************************************************//**
 * \brief runs an other program by calling exec and fork.
 *
 * \param[in] e Template argument. Must be callable as function without parameters
 *    and return an errorstatus in case of failure (preferably \c errno).
 *    If no error occures in \a e, the function should not return as the process
 *    is replaced by a new process.
 *******************************************************************************/
template<typename T>
inline void helper::exec(T e, std::vector<char> &stdout_, std::vector<char> &stderr_, int *status_) {

    stdout_.resize(0);
    stderr_.resize(0);
    int status__ = 0;
    int &status = status_ ? *status_ : status__;
    status = 0;
    int errno_;
    int ret;

    helper::pipe fd_error, fd_stdout, fd_stderr;

    std::size_t stdout_size = 0;
    std::size_t stderr_size = 0;

    pid_t pid = fork();
    if (pid < 0) {
        errno_ = errno;
        fd_error.close_nothrow();
        fd_stdout.close_nothrow();
        fd_stderr.close_nothrow();
        throw std::runtime_error(std::string("fork failed: ") + strerror(errno_));
    } else if (0 == pid) {
        // child process.
        std::string error_msg = "unknown error.";
        try {

            // set the close on exec flag for the error pipe.
            {
                // get the flags of the file descriptor.
                int fd_ = fd_error[1];
                if ((ret = fcntl(fd_, F_GETFD)) < 0) {
                    errno_ = errno;
                    throw std::runtime_error(std::string("fcntl(F_GETFD) of error pipe: ") + strerror(errno_));
                }

                // make the write end close-on-exec
                int fd_flags = ret | FD_CLOEXEC;
                if (fcntl(fd_, F_SETFD, fd_flags) == -1) {
                    errno_ = errno;
                    throw std::runtime_error(std::string("fcntl(F_SETFD) of error pipe: ") + strerror(errno_));
                }
            }
            fd_error.close_fd(0);

            // Close the STDIN of the cild.
            if (close(STDIN_FILENO) != 0) {
                errno_ = errno;
                throw std::runtime_error(std::string("could not close STDIN: ") + strerror(errno_));
            }

            // connect the STDOUT of the child with the parent.
            if (dup2(fd_stdout[1], STDOUT_FILENO) == -1) {
                errno_ = errno;
                throw std::runtime_error(std::string("error redirect STDOUT to parent: ") + strerror(errno_));
            }
            fd_stdout.close();

            // connect the STDERR of the child with the parent.
            if (dup2(fd_stderr[1], STDERR_FILENO) == -1) {
                errno_ = errno;
                throw std::runtime_error(std::string("error redirect STDERR to parent: ") + strerror(errno_));
            }
            fd_stderr.close();

            errno = 0;
            try {
                errno_ = e();
                error_msg = std::string("exec returned with errno (") + strerror(errno_) + ").";
            } catch (std::exception &e) {
                errno_ = errno;
                error_msg = std::string("exec throws exception (") + e.what() + ").";
            } catch (...) {
                errno_ = errno;
                error_msg = std::string("exec throws unknown exception: (") + strerror(errno_) + ").";
            }
        } catch (std::exception &e) {
            // try to write the error message.
            error_msg = e.what();
        } catch (...) {
        }
        error_msg = "child process: " + error_msg;
        // try to write the error message to the parent.
        try {
            std::cout << __FILE__ << ": " << __LINE__ << " " << error_msg << std::endl;
            std::cout <<

            write(fd_error[1], error_msg.c_str(), error_msg.size()+1)

            << std::endl;
        } catch (...) {
            std::cout << __FILE__ << ": " << __LINE__ << std::endl;
        }
        std::cerr << __FILE__ << ": " << __LINE__ << std::endl;
        // che child exits with error code -1.
        exit(-1);
    } else {
        // parent process.
        pid_t wpid;

        const std::size_t BUFFSIZE = 512;
        ssize_t cread = 0;
        try {
            try {
                // Close the write end of the pipes.
                fd_error.close_fd(1);
                fd_stdout.close_fd(1);
                fd_stderr.close_fd(1);

                // Read the error message from the error-pipe.
                int fd = fd_error[0];
                std::size_t ccbuff = 0;
                std::vector<char> cbuff;
                cread = 0;
                do {
                    ccbuff += cread;
                    cbuff.resize(ccbuff+BUFFSIZE, 0);
                    do {
                        cread = read(fd, &cbuff[ccbuff], BUFFSIZE);
                    } while (cread==-1 && errno==EINTR); // In case of returned by a signal: repeat reading.
                } while (cread > 0);
                if (cread == -1) {
                    errno_ = errno;
                    throw std::runtime_error(std::string("error reading the error message fron the child: ") + strerror(errno_));
                }
                // Just to be sure, that the buffer is null teminated.
                cbuff.push_back(0); cbuff[ccbuff+=cread] = 0;
                if (ccbuff > 0) {
                    throw std::runtime_error(&cbuff[0]);
                }

                fd_error.close();
            } catch (...) {
                fd_error.close_nothrow();
                fd_stdout.close_nothrow();
                fd_stderr.close_nothrow();
                throw;
            }

            // At this point, the child process successfully did exec.
            // Now read the STDOUT and STDERR...
            fd_set rfds;
            bool fd_stdout_open = true;
            bool fd_stderr_open = true;
            int fd_stdout_errno = 0;
            int fd_stderr_errno = 0;
            int fd_stderr_ = fd_stderr[0];
            int fd_stdout_ = fd_stdout[0];
            int nfds = (fd_stderr_>fd_stdout_ ? fd_stderr_ : fd_stdout_) + 1;
            int i;
            // Now read from the pipes into the vectors.
            // Use select to read only when there is really data available.
            // Otherwise it likely will deadlock.
            try {
                // Repeat as long there are open streams (pipes).
                do {
                    // Init the bit-mask.
                    FD_ZERO(&rfds);
                    if (fd_stdout_open) { FD_SET(fd_stdout_, &rfds); }
                    if (fd_stderr_open) { FD_SET(fd_stderr_, &rfds); }
                    i = select(nfds, &rfds, NULL, NULL, NULL);
                    if (i < 0) {
                        // error in select.
                        errno_ = errno;
                        if (i==-1 && errno_ == EINTR) { continue; }
                        throw std::runtime_error(std::string("error waiting for reading streams (select): ") + strerror(errno_));
                    }
                    if (FD_ISSET(fd_stdout_, &rfds)) {
                        // Read from stream.
                        stdout_.resize(stdout_size+BUFFSIZE, 0);
                        cread = read(fd_stdout_, &stdout_[stdout_size], BUFFSIZE);
                        if (cread == -1) {
                            if ((errno_=errno) != EINTR) {
                                fd_stdout_open = false;
                                fd_stdout_errno = errno_;
                                fd_stdout.close_nothrow();
                            }
                        } else if (cread == 0) {
                            // EOF
                            fd_stdout_open = false;
                            fd_stdout.close_nothrow();
                        } else {
                            stdout_size += cread;
                        }
                    }
                    if (FD_ISSET(fd_stderr_, &rfds)) {
                        // Read from stream.
                        stderr_.resize(stderr_size+BUFFSIZE, 0);
                        cread = read(fd_stderr_, &stderr_[stderr_size], BUFFSIZE);
                        if (cread == -1) {
                            if ((errno_=errno) != EINTR) {
                                fd_stderr_open = false;
                                fd_stderr_errno = errno_;
                                fd_stderr.close_nothrow();
                            }
                        } else if (cread == 0) {
                            // EOF
                            fd_stderr_open = false;
                            fd_stderr.close_nothrow();
                        } else {
                            stderr_size += cread;
                        }
                    }
                } while ((fd_stdout_open||fd_stderr_open));
            } catch (...) {
                fd_stdout.close_nothrow();
                fd_stderr.close_nothrow();
                throw;
            }
            if (fd_stdout_errno!=0 && fd_stderr_errno!=0) {
                throw std::runtime_error(std::string("Error reading STDOUT and STDERR from child process: OUT: ") + strerror(fd_stdout_errno) + "; ERR: " + strerror(fd_stderr_errno) + ".");
            } else if (fd_stdout_errno != 0) {
                throw std::runtime_error(std::string("Error reading STDOUT from child process: ") + strerror(fd_stdout_errno) + ".");
            } else if (fd_stderr_errno != 0) {
                throw std::runtime_error(std::string("Error reading STDERR from child process: ") + strerror(fd_stderr_errno) + ".");
            }
        } catch (...) {
            do {
                wpid = waitpid(pid, &status, 0);
            } while (wpid==-1 && errno==EINTR); // if a signal woke the parent, sleep again....
            throw;
        }

        // wait for child process.
        do {
            wpid = waitpid(pid, &status, 0);
        } while (wpid==-1 && errno==EINTR); // if a signal woke the parent, sleep again....
        if (wpid == -1) {
            errno_ = errno;
            throw std::runtime_error(std::string("Error waiting for child: ") + strerror(errno) + ".");
        }
    }

    stdout_.resize(stdout_size);
    stderr_.resize(stderr_size);
}




#endif










