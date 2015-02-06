/****************************************************************************//**
 * \file ostream_swapper.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    24.01.2008
 * Taken from http://groups.google.com/group/borland.public.cppbuilder.language/browse_thread/thread/52e6bab7756110b0/dfc46e2bfaaa6cbc?lnk=st&#dfc46e2bfaaa6cbc
 *******************************************************************************/
#ifndef ___INCLUDE_OSTREAM_SWAPPER_HPP__
#define ___INCLUDE_OSTREAM_SWAPPER_HPP__



#include <ostream>


namespace helper {

class ostream_swapper {

public:
    ostream_swapper(std::ostream &stream, std::ostream &sreplace)
        : swapped_back(false), rdbuf_(stream.rdbuf()), stream_(stream) {
        stream.rdbuf(sreplace.rdbuf());
    }
    ~ostream_swapper() {
        this->swap_back();
    }
    void swap_back() {
        if (!swapped_back) {
            stream_.rdbuf(rdbuf_);
            swapped_back = true;
        }
    }
private:
    bool swapped_back;
    std::streambuf *rdbuf_;
    std::ostream & stream_;
};


}


#endif

