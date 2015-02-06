/****************************************************************************//**
 * \file filesystem.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    05.03.2007
 *******************************************************************************/
#ifndef ___INCLUDE_HELPER_FILESYSTEM_HPP__
#define ___INCLUDE_HELPER_FILESYSTEM_HPP__



#ifndef NO_BOOST_FILESYSTEM

#include <boost/filesystem.hpp>


namespace helper {
    namespace fs = boost::filesystem;
} // namespace helper


#else

#error not_yet_implemented

// create a vary basic implementation of the boost/filesystem api.

namespace helper {
namespace fs {


    class path: public std::string {
    };

    inline bool is_regular( const path & ph ) {
        return false;
    }

    inline bool is_directory( const path & ph ) {
        return false;
    }

    inline bool create_directory( const path & dir_ph ) {
        return false;
    }

    inline bool remove( const path & ph ) {
        return false;
    }

    inline bool exists( const path & ph ) {
        return false;
    }

} // namespace fs
} // namespace helper

#endif






#endif


