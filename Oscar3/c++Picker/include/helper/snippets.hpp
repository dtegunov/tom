/****************************************************************************//**
 * \file snippets.hpp
 * \brief Small functions for general purpose.
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    25.01.2008
 *******************************************************************************/
#ifndef ___INCLUDE_SNIPPETS_HPP__
#define ___INCLUDE_SNIPPETS_HPP__


#include <vector>
#include <string>
#include <iosfwd>


#include <boost/shared_ptr.hpp>





namespace helper {

template<typename T>
bool str2range(const std::string &instr, std::vector<T> &range, char delimiter, char delimiter_range);




std::string now2str();


template<typename T, typename _Compare>
void unify_shared_vector(std::vector<boost::shared_ptr<T> > &v, const _Compare &c = _Compare());



template<typename T, typename TIDX>
void subset(const std::vector<T> &v, const std::vector<TIDX> &vidx, std::vector<T> &vsub, bool unique, bool sort);


inline std::vector<const char *> svector2c(std::vector<std::string> &v) {
    std::vector<const char *> r(v.size()+1);
    std::vector<const char *>::iterator rit = r.begin();
    std::vector<std::string>::iterator vit = v.begin();
    for (; vit!=v.end(); vit++, rit++) {
        *rit = vit->c_str();
    }
    *rit = NULL;
    return r;
}



std::string split_line(std::string &s, std::size_t max_line_length);
void print_parameters(std::ostream &clog, const std::string &name, const std::string &text, const std::string &prefix, std::size_t max_line_length);



}




#endif


