/***********************************************************************//**
 * \file snippets.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    05.02.2007
 **************************************************************************/
#include "helper/snippets.hpp"


#include <set>
#include <map>
#include <ostream>


#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>







/****************************************************************************//**
 *
 *
 *******************************************************************************/
std::string helper::now2str() {
    time_t now = time(NULL);
    std::string s(asctime(localtime(&now)));
    s.erase(s.end()-1);
    return s;
}






/****************************************************************************//**
 *
 *******************************************************************************/
std::string helper::split_line(std::string &s, std::size_t max_line_length) {
    #define DELIMITER " \n\t"
    std::string s0;
    if (max_line_length < 1) {

    } else if (max_line_length >= s.size()) {
        s0.swap(s);
        s.clear();
    } else {
        std::string s_;
        s.substr(0, max_line_length+1).swap(s_);
        std::string::size_type idx = s_.find_last_of(DELIMITER);
        if (idx == std::string::npos) {
            idx = s.find_first_of(DELIMITER, max_line_length-1);
        }

        if (idx != std::string::npos) {
            s.substr(0, idx).swap(s0);
            s.erase(0, idx+1);
        } else {
            s0.swap(s);
            s.clear();
        }
    }

    return s0;
    #undef DELIMITER
}

/****************************************************************************//**
 *
 *******************************************************************************/
void helper::print_parameters(std::ostream &clog, const std::string &name, const std::string &text, const std::string &prefix, std::size_t max_line_length) {

    std::string s(text);
    std::size_t i;
    if (!name.empty()) {
        i = name.size() + 2;
        clog << name << ": " << helper::split_line(s, i>max_line_length ? 0 : max_line_length-i) << "\n";
    }

    i = prefix.size();
    if (max_line_length < i) {
        max_line_length = 1;
    } else {
        max_line_length -= i;
    }

    while (!s.empty()) {
        clog << prefix << helper::split_line(s, max_line_length) << "\n";
    }

}




/****************************************************************************//**
 * \param[in] instr The input string to be parsed to an integer range.
 *   The single elements are separated by ','. '-' specifies a range (with both
 *   numbers included. Negative numbers are not allowed. Spaces are not allowed.
 * \param[out] range The found numbers. Sorted in accending order.
 * \return True if the string could be successfully parsed. Otherwise false.
 *
 * Each value comes only one time in the result. Empty ranges occure never.
 * The output \a range will be sorted.
 *******************************************************************************/
template<typename T>
bool helper::str2range(const std::string &instr, std::vector<T> &range, char delimiter, char delimiter_range) {

    range.resize(0);
    if (instr.empty()) {
        return false;
    }

    std::set<T> srange;
    std::string temp;
    std::string str = instr;

    std::string::size_type pos;
    T val1, val2, i;
    try {
        do {
            if ((pos=str.find(delimiter, 0)) == std::string::npos) {
                pos = str.size();
            }
            temp = str.substr(0, pos);
            str.erase(0, pos + 1);
            pos = temp.find(delimiter_range, 0);
            if (pos == std::string::npos) {
                srange.insert(boost::lexical_cast<T>(temp));
            } else if (temp.find(delimiter_range, pos+1) == std::string::npos) {
                val1 = boost::lexical_cast<T>(temp.substr(0, pos));
                temp.erase(0, pos+1);
                val2 = boost::lexical_cast<T>(temp);
                for (i=val1; i<val2 || (!(i<val2)&&!(val2<i)); i++) { // Use only the < operator of the template-type T.
                    srange.insert(i);
                }
            } else {
                return false;
            }
        } while (!str.empty());
    } catch (boost::bad_lexical_cast &e) {
        return false;
    }
    range.reserve(srange.size());
    for (std::set<std::size_t>::const_iterator it=srange.begin(); it!=srange.end(); it++) {
        range.push_back(*it);
    }
    return true;
}




/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename T, typename TIDX>
void helper::subset(const std::vector<T> &v, const std::vector<TIDX> &vidx, std::vector<T> &vsub, bool unique, bool sort) {

    std::vector<TIDX> vsorted;
    const std::vector<TIDX> *pvidx = &vidx;
    if (unique) {
        std::set<TIDX> vidx_sorted;
        if (sort) {
            vidx_sorted.insert(vidx.begin(), vidx.end());
            vsorted.assign(vidx_sorted.begin(), vidx_sorted.end());
        } else {
            vsorted.reserve(vidx.size());
            for (typename std::vector<TIDX>::const_iterator it=vidx.begin(); it!=vidx.end(); it++) {
                if (vidx_sorted.insert(*it).second) {
                    vsorted.push_back(*it);
                }
            }
        }
        pvidx = &vsorted;
    } else if (sort) {
        std::multiset<TIDX> vidx_sorted(vidx.begin(), vidx.end());
        vsorted.assign(vidx_sorted.begin(), vidx_sorted.end());
        pvidx = &vsorted;
    }

    vsub.resize(0);
    vsub.reserve(pvidx->size());
    for (typename std::vector<TIDX>::const_iterator it=pvidx->begin(); it!=pvidx->end(); it++) {
        vsub.push_back(v[*it]);
    }
}



/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename T, typename _Compare>
void helper::unify_shared_vector(std::vector<boost::shared_ptr<T> > &v, const _Compare &c) {

    std::map<T *, std::vector<std::size_t>, _Compare> m(c);

    typename std::vector<boost::shared_ptr<T> >::const_iterator vit;
    typename std::vector<std::size_t>::const_iterator vit2;
    typename std::map<T *, std::vector<std::size_t>, _Compare>::const_iterator mit;
    std::size_t i;

    for (i=0, vit=v.begin(); vit!=v.end(); vit++, i++) {
        std::vector<std::size_t> &vv = m[vit->get()];
        vv.push_back(i);
    }

    for (mit=m.begin(); mit!=m.end(); mit++) {
        const std::vector<std::size_t> &vv = mit->second;
        boost::shared_ptr<T> &p = v[vv[0]];
        for (vit2=vv.begin()+1; vit2!=vv.end(); vit2++) {
            v[*vit2] = p;
        }
    }
}




// template instantiations.
template bool helper::str2range(const std::string &instr, std::vector<std::size_t> &range, char delimiter, char delimiter_range);



#include "tom/core/volume.hpp"
template void helper::unify_shared_vector<tom::Volume<double>, tom::volume_less<double> >(std::vector<boost::shared_ptr<tom::Volume<double> > > &v, const tom::volume_less<double> &c);



template void helper::subset<std::string, std::size_t>(const std::vector<std::string> &v, const std::vector<std::size_t> &vidx, std::vector<std::string> &vsub, bool unique, bool sort);

