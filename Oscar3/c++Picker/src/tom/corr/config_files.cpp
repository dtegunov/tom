/****************************************************************************//**
 * \file config_files.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    06.02.2008
 *******************************************************************************/
#include "tom/corr/config_files.hpp"




#include <sstream>
#include <map>
#include <set>
#include <string>
#include <fstream>
#include <iomanip>

#include <boost/lexical_cast.hpp>
#include "boost/tuple/tuple.hpp"

#include "helper/snippets.hpp"
#include "helper/GetPot"


#include "tom/corr/filename_generators.hpp"



namespace {
/****************************************************************************//**
 *
 *******************************************************************************/
void str_trim_spaces(std::string &str) {

    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if((std::string::npos == startpos) || (std::string::npos == endpos)) {
        str = "";
    } else {
        str = str.substr(startpos, endpos-startpos+1);
    }
}

}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::corr::write_filelist(const std::string &filename, const std::vector<std::string> &ffilenames, const std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > > &wedge_desc) {

    std::ofstream f(filename.c_str(), std::ios_base::trunc | std::ios_base::out);

    if (!f) {
        std::ostringstream ss; ss << "ERROR: writing filelist: Could not open \"" << filename << "\"."; throw std::runtime_error(ss.str());
    }

    if (ffilenames.size() != wedge_desc.size()) {
        std::ostringstream ss; ss << "ERROR: writing filelist: length of filenames and wedge is not equal."; throw std::runtime_error(ss.str());
    }

    std::vector<std::string>::const_iterator itf = ffilenames.begin();
    typename std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > >::const_iterator itw = wedge_desc.begin();

    for (; itf!=ffilenames.end(); itf++, itw++) {
        f << *itf << "\n" <<
             (*itw)->getConfigEntry() << std::endl;
    }
}








/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::corr::parse_filelist(const std::string &filename, std::vector<std::string> &ffilenames, std::vector<boost::shared_ptr<tom::WedgeDescriptor<T> > > &wedge_desc, std::size_t binning) {

    ffilenames.resize(0);
    wedge_desc.resize(0);

    std::ifstream f(filename.c_str());
    if (!f.good()) {
        std::ostringstream ss; ss << "Parsing error: Could not open \"" << filename << "\".";
        throw std::runtime_error(ss.str());
    }

    std::size_t linecnt=0;
    std::string s1, s2, s21, s22, s23;
    std::istringstream iss;

    if (!binning) {
        binning = 1;
    }

    do {
        std::getline(f, s1);

        linecnt++;

        //std::cout << std::setw(3) << linecnt << ": \"" << s1 << "\"" << std::endl;

        if (f.good()) {

            if (s1.empty()) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": " << linecnt << ": Empty filename.";
                throw std::runtime_error(ss.str());
            }


            std::getline(f, s2);
            linecnt++;

            //std::cout << std::setw(3) << linecnt << ": \"" << s2 << "\"" << std::endl;

            ffilenames.push_back(s1);

            try {
                wedge_desc.push_back(boost::shared_ptr<tom::WedgeDescriptor<T> >(tom::WedgeDescriptor<T>::fromConfigEntry(s2, binning)));
            } catch (std::runtime_error &e) {
                std::ostringstream ss;
                ss << "Parsing error: \"" << filename << "\": " << linecnt << ": " << e.what();
                throw std::runtime_error(ss.str());
            }
        } else if (!s1.empty()) {
            std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": " << linecnt << ": file does not contain an even number of lines (with filename, wedge alternating).";
            throw std::invalid_argument(ss.str());
        }
    } while (f.good());
    //assert(ffilenames.size() == wedge_desc.size());

    if (ffilenames.empty()) {
        throw std::runtime_error("The filelist " + filename + " is empty.");
    }
}





/****************************************************************************//**
 *
 *******************************************************************************/
void tom::corr::write_angleslist(const std::string &filename, std::size_t ntemplates, std::size_t nparticles, const std::vector<boost::shared_ptr<tom::Volume<double> > > &angles, const tom::av4::filename_generator &f_gen, unsigned long iit) {


    std::ofstream f(filename.c_str(), std::ios_base::trunc | std::ios_base::out);

    if (!f) {
        std::ostringstream ss; ss << "ERROR: writing angles list: Could not open \"" << filename << "\"."; throw std::runtime_error(ss.str());
    }

    if (ntemplates*nparticles != angles.size()) {
        throw std::runtime_error("ERROR: angles list has the wrong size.");
    }

    std::size_t itemplate, iparticle, i;
    std::map<const tom::Volume<double> *, std::vector<std::size_t>, tom::volume_less<double> > vmap;
    for (i=0; i<ntemplates*nparticles; i++) {
        vmap[angles[i].get()].push_back(i);
    }

    std::map<const tom::Volume<double> *, std::vector<std::size_t>, tom::volume_less<double> >::const_iterator it;
    std::vector<std::size_t>::const_iterator vit;
    std::string f_angles;
    for (it=vmap.begin(); it!=vmap.end(); it++) {
        const tom::Volume<double> *const &vol = it->first;
        if (vol) {
            f_gen.get_angles_file(iit, *vol).swap(f_angles);
            const std::vector<std::size_t> &v = it->second;
            vol->write_to_em(f_angles, NULL);
            for (vit=v.begin(); vit!=v.end(); vit++) {
                itemplate = *vit / nparticles;
                iparticle = *vit % nparticles;
                f << std::setw(3) << std::setfill(' ') << itemplate << "  " <<
                     std::setw(5) << std::setfill(' ') << iparticle << "  " <<
                     f_angles << std::endl;
            }
        }
    }
}




/****************************************************************************//**
 *
 *******************************************************************************/
void tom::corr::parse_angleslist(const std::string &filename, std::size_t ntemplates, std::size_t nparticles, std::vector<boost::shared_ptr<tom::Volume<double> > > &angles) {


    std::ifstream f(filename.c_str());
    if (!f.good()) {
        std::ostringstream ss; ss << "Parsing error: Could not open \"" << filename << "\".";
        throw std::runtime_error(ss.str());
    }

    std::string s, s1, s2, s3;
    std::size_t linecnt = 0;
    std::istringstream iss;
    std::vector<std::size_t> ranget, rangep;
    std::map<std::size_t, std::pair<std::size_t, std::string> > anglesm;
    std::size_t i;
    do {
        std::getline(f, s);
        linecnt++;
        if (!s.empty() && s[0]!='#') {
            iss.clear();
            iss.str(s);
            s1.clear(); s2.clear();
            iss >> s1 >> s2;
            std::getline(iss, s3);

            if (!helper::str2range(s1, ranget, ',', '-') || ranget.empty() || !helper::str2range(s2, rangep, ',', '-') || rangep.empty()) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": " << linecnt << ": the line must begin with two integer values referencing the index of template and particle respectively (zero based).";
                throw std::runtime_error(ss.str());
            }
            if (ranget[ranget.size()-1]>=ntemplates || rangep[rangep.size()-1]>=nparticles) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": " << linecnt << ": There are only " << ntemplates << "x" << nparticles << " volumes. The index is out of bound.";
                throw std::runtime_error(ss.str());
            }

            for (std::vector<std::size_t>::const_iterator it=ranget.begin(); it!=ranget.end(); it++) {
                for (std::vector<std::size_t>::const_iterator ip=rangep.begin(); ip!=rangep.end(); ip++) {
                    i = *it * nparticles + *ip;
                    if (anglesm.find(i) == anglesm.end()) {
                        std::pair<std::size_t, std::string> &ai = anglesm[i];
                        ai.first = linecnt;
                        ai.second = s3;
                    } else {
                        std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": " << linecnt << ": Angle defined for (" << *it << "," << *ip << ") but it appeared before in line " << anglesm[i].first << ".";
                        throw std::runtime_error(ss.str());
                    }
                }
            }

        }
    } while (f.good());

    angles.assign(ntemplates*nparticles, boost::shared_ptr<tom::Volume<double> >());
    tom::Volume<double> *pangles;
    std::auto_ptr<tom::Volume<double> > apangles;
    std::map<std::string, boost::shared_ptr<tom::Volume<double> > > volm;
    for (std::map<std::size_t, std::pair<std::size_t, std::string> >::const_iterator it=anglesm.begin(); it!=anglesm.end(); it++) {
        s1 = it->second.second;
        str_trim_spaces(s1);
        if (s1.size() >= 2 && s1[0]=='"' && s1[s1.size()-1]=='"') {
            s1 = s1.substr(1, s1.size()-2);
        }
        //std::cout << " angle: " << std::setw(3) << (it->first/nparticles) << "x" << std::setw(3) << (it->first%nparticles) << ": line " << std::setw(3) << it->second.first << " = \"" << s1 << "\"" << std::endl;

        if (volm.find(s1) == volm.end()) {
            try {
                tom::read_from_em<double>(pangles, s1, NULL,NULL,NULL, NULL, NULL,NULL);
                apangles.reset(pangles);
            } catch (std::exception &e) {
                std::ostringstream ss; ss << "Could not load the angle list \"" << s1 << "\" specified in " << filename << ":" << it->second.first << ": " << e.what() << ".";
                throw std::runtime_error(ss.str());
            }
            if (apangles->getSizeX()!=3 || apangles->getSizeZ()!=1) {
                std::ostringstream ss; ss << "The angle list in em-file \"" << s1 << "\" is not of size 3xNx1 as expected (phi,psi,theta in deg). Specified in " << filename << ":" << it->second.first << ".";
                throw std::runtime_error(ss.str());
            }
            volm[s1] = boost::shared_ptr<tom::Volume<double> >(apangles.release());
        }
        assert(it->first < ntemplates*nparticles);
        angles[it->first] = volm[s1];
    }

}





/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::corr::parse_peaklist(const std::string &filename, std::vector<tom::cc_peak<T> > &peak_list, std::vector<helper::triple<double, double, double> > &anglesv, std::size_t ntemplates, std::size_t nparticles) {

    std::ifstream f(filename.c_str());
    if (!f) {
        std::ostringstream ss; ss << "Parsing error: Could not open \"" << filename << "\".";
        throw std::runtime_error(ss.str());
    }

    std::string s;
    std::vector<std::string> line_tok(10);
    std::istringstream iss;
    std::size_t linecnt=0, i;

    std::map<std::pair<std::size_t, std::size_t>, boost::tuple<std::size_t, tom::cc_peak<T>, helper::triple<double,double,double> > > map;

    helper::triple<double,double,double> angle;
    tom::cc_peak<T> peak;
    std::pair<std::size_t, std::size_t> ipair;
    std::pair<typename std::map<std::pair<std::size_t, std::size_t>, boost::tuple<std::size_t, tom::cc_peak<T>, helper::triple<double,double,double> > >::iterator, bool> rinsert;

    do {
        linecnt++;
        std::getline(f, s);
        if (s.empty() || s[0]=='#') {
        } else {
            iss.clear();
            iss.str(s);
            i = 0;
            while (i<line_tok.size() && iss >> line_tok[i]) { i++; }
            if (i != line_tok.size()) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": to few tokens in line " << linecnt << "."; throw std::runtime_error(ss.str());
            }
            iss >> std::ws;
            if (!iss.eof()) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": to many tokens in line " << linecnt << "."; throw std::runtime_error(ss.str());
            }

            try {
                i = 0;
                ipair.first = boost::lexical_cast<std::size_t>(line_tok[i++]);
                ipair.second = boost::lexical_cast<std::size_t>(line_tok[i++]);
                peak.val = boost::lexical_cast<typename tom::cc_peak<T>::val_type>(line_tok[i++]);
                peak.x = boost::lexical_cast<typename tom::cc_peak<T>::idx_type>(line_tok[i++]);
                peak.y = boost::lexical_cast<typename tom::cc_peak<T>::idx_type>(line_tok[i++]);
                peak.z = boost::lexical_cast<typename tom::cc_peak<T>::idx_type>(line_tok[i++]);
                peak.angle_idx = boost::lexical_cast<typename tom::cc_peak<T>::angle_idx_type>(line_tok[i++]);
                angle.x = boost::lexical_cast<double>(line_tok[i++]);
                angle.y = boost::lexical_cast<double>(line_tok[i++]);
                angle.z = boost::lexical_cast<double>(line_tok[i++]);
            } catch (boost::bad_lexical_cast &e) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": parsing line " << linecnt << "."; throw std::runtime_error(ss.str());
            }

            if (ipair.first  >= ntemplates) { std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": template index to large (line " << linecnt << ")."; throw std::runtime_error(ss.str()); }
            if (ipair.second >= nparticles) { std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": particle index to large (line " << linecnt << ")."; throw std::runtime_error(ss.str()); }

            rinsert = map.insert(std::make_pair(ipair, boost::tuples::make_tuple(linecnt, peak, angle)));
            if (!rinsert.second) {
                std::ostringstream ss; ss << "Parsing error: \"" << filename << "\": The peak for the combination " << ipair.first << "x" << ipair.second << " was first defined in line " << boost::tuples::get<0>(rinsert.first->second) << " and again in line " << linecnt << "."; throw std::runtime_error(ss.str());
            }
        }
    } while (f);

    typename std::map<std::pair<std::size_t, std::size_t>, boost::tuple<std::size_t, tom::cc_peak<T>, helper::triple<double,double,double> > >::const_iterator mapit;

    peak_list.assign(ntemplates*nparticles, tom::cc_peak<T>());
    anglesv.assign(ntemplates*nparticles, helper::triple<double, double, double>(0,0,0));
    for (mapit=map.begin(); mapit!=map.end(); mapit++) {
        assert(mapit->first.first<ntemplates && mapit->first.second<nparticles);

        i = mapit->first.first*nparticles + mapit->first.second;
        peak_list[i] = boost::tuples::get<1>(mapit->second);
        anglesv[i] = boost::tuples::get<2>(mapit->second);
    }
}






/****************************************************************************//**
 *
 *******************************************************************************/
void tom::avg::parse_average(   const std::string &filename,
                                std::vector<tom::avg::st_average_result> &av) {

    av.resize(0);



    std::ifstream f(filename.c_str());
    if (!f) { throw std::runtime_error("ERROR parsing average log: Could not open \"" + filename + "\"."); }

    std::string s;
    std::istringstream iss;
    std::size_t linecnt=0, i, j, k, itemplate;
    double x, y;
    std::string s1, s2;
    long signed int ntemplates = -1;

    std::vector<std::size_t> av_use_idx;
    std::ostringstream ss;

    std::map<std::size_t, std::size_t> iparticle_line;

    do {
        linecnt++;
        std::getline(f, s);
        if (s.empty() || s[0]=='#') {
        } else {
            iss.clear();
            iss.str(s);
            s1.clear(); iss >> std::ws >> s1;
            if (s1 == "NTEMPLATES") {
                if (ntemplates != -1) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": option NTEMPLATES occures several times.";
                    throw std::runtime_error(ss.str());
                }
                s1.clear(); iss >> std::ws >> s1;
                try {
                    ntemplates = boost::lexical_cast<long signed int>(s1);
                } catch (boost::bad_lexical_cast &e) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": option NTEMPLATES expects the number of classes as parameter.";
                    throw std::runtime_error(ss.str());
                }
                if (ntemplates < 0) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": option NTEMPLATES expects the number of classes as parameter (must be positive).";
                    throw std::runtime_error(ss.str());
                } else if (static_cast<std::size_t>(ntemplates) > av.size()) {
                    av.resize(ntemplates);
                } else if (static_cast<std::size_t>(ntemplates) < av.size()) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": option NTEMPLATES specifies " << ntemplates << " classes, but there are some fields with higher index.";
                    throw std::runtime_error(ss.str());
                }
                s1.clear(); iss >> std::ws >> s1;
                if (!s1.empty()) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": To many tokens. Format is: \"NTEMPLATES n\".";
                    throw std::runtime_error(ss.str());
                }
            } else if (s1 == "NUM") {
                s1.clear(); iss >> std::ws >> s1;
                try {
                    itemplate = boost::lexical_cast<std::size_t>(s1);
                } catch (boost::bad_lexical_cast &e) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": expects the index of the template as first parameter.";
                    throw std::runtime_error(ss.str());
                }
                if (ntemplates == -1) {
                    if (itemplate >= av.size()) {
                        av.resize(itemplate+1);
                    }
                } else if (itemplate >= static_cast<std::size_t>(ntemplates)) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": given template index of " << itemplate << ", but there are only " << ntemplates << " templates.";
                    throw std::runtime_error(ss.str());
                }
                s1.clear(); iss >> std::ws >> s1;
                try {
                    i = boost::lexical_cast<std::size_t>(s1);
                } catch (boost::bad_lexical_cast &e) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": NUM expects the number of particles as second parameter.";
                    throw std::runtime_error(ss.str());
                }
                av_use_idx.resize(i);
                for (j=0; j<i; j++) {
                    s1.clear(); iss >> std::ws >> s1;
                    try {
                        av_use_idx[j] = k = boost::lexical_cast<std::size_t>(s1);
                    } catch (boost::bad_lexical_cast &e) {
                        ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": NUM expects the indices of the " << i << " particles.";
                        throw std::runtime_error(ss.str());
                    }
                    std::pair<std::map<std::size_t, std::size_t>::iterator, bool> it = iparticle_line.insert(std::make_pair(k, linecnt));
                    if (!it.second) {
                        ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": index " << k << " appeaed already before (in line " << it.first->second << ").";
                        throw std::runtime_error(ss.str());
                    }
                }
                s1.clear(); iss >> std::ws >> s1;
                if (!s1.empty()) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": To many tokens. Format is: NUM ITEMPLATE NPARTICLE IPARTICLE1 IPARTICLE2 ... IPARTICLEN.";
                    throw std::runtime_error(ss.str());
                }
                av[itemplate].use_idx.swap(av_use_idx);
            } else if (s1=="CC_N" || s1=="CC_1" || s1=="CC_SQRTN" || s1=="CC_S") {
                s2.swap(s1);
                s1.clear(); iss >> std::ws >> s1;
                try {
                    itemplate = boost::lexical_cast<std::size_t>(s1);
                } catch (boost::bad_lexical_cast &e) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": expects the index of the template as first parameter.";
                    throw std::runtime_error(ss.str());
                }
                if (ntemplates == -1) {
                    if (itemplate >= av.size()) {
                        av.resize(itemplate+1);
                    }
                } else if (itemplate >= static_cast<std::size_t>(ntemplates)) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": given template index of " << itemplate << ", but there are only " << ntemplates << " templates.";
                    throw std::runtime_error(ss.str());
                }
                s1.clear(); iss >> std::ws >> s1;
                try {
                    x = boost::lexical_cast<double>(s1);
                } catch (boost::bad_lexical_cast &e) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": floating point number expected as second parameter.";
                    throw std::runtime_error(ss.str());
                }
                if (s2=="CC_S") {
                    av[itemplate].ccval_s.reset(new double(x));
                } else if (s2=="CC_1") {
                    av[itemplate].ccval_1.reset(new double(x));
                } else if (s2=="CC_SQRTN") {
                    av[itemplate].ccval_sqrtn.reset(new double(x));
                } else if (s2=="CC_N") {
                    av[itemplate].ccval_n.reset(new double(x));
                } else {
                    assert(0);
                }
            } else if (s1 == "FSC") {
                s1.clear(); iss >> std::ws >> s1;
                try {
                    itemplate = boost::lexical_cast<std::size_t>(s1);
                } catch (boost::bad_lexical_cast &e) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": expects the index of the template as first parameter.";
                    throw std::runtime_error(ss.str());
                }
                if (ntemplates == -1) {
                    if (itemplate >= av.size()) {
                        av.resize(itemplate+1);
                    }
                } else if (itemplate >= static_cast<std::size_t>(ntemplates)) {
                    ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": given template index of " << itemplate << ", but there are only " << ntemplates << " templates.";
                    throw std::runtime_error(ss.str());
                }
                std::vector<std::pair<double, double> > v;
                s1.clear(); iss >> std::ws >> s1;
                do {
                    try {
                        //std::cerr << "_1__" << s1 << "____" << std::endl;
                        x = boost::lexical_cast<double>(s1);
                        s1.clear(); iss >> std::ws >> s1;
                        //std::cerr << "_2__" << s1 << "____" << std::endl;
                        y = boost::lexical_cast<double>(s1);
                    } catch (boost::bad_lexical_cast &e) {
                        ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": floating point numbers expected.";
                        throw std::runtime_error(ss.str());
                    }
                    s1.clear(); iss >> std::ws >> s1;
                    v.push_back(std::make_pair(x, y));
                } while(!s1.empty());
                av[itemplate].fsc.swap(v);
            } else {
                std::ostringstream ss; ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": line start \"" << s1 << "\" not recognized."; throw std::runtime_error(ss.str());
            }
        }
    } while (f);
    if (ntemplates == -1) {
        ss << "ERROR parsing average log: \"" << filename << "\":" << linecnt << ": missing parameter NTEMPLATES.";
        throw std::runtime_error(ss.str());
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
void tom::avg::write_average_to_log(  std::ostream &clog,
                            const std::vector<std::string> &ftemplates,
                            const std::vector<std::string> &fparticles,
                            const std::vector<tom::avg::st_average_result> &av,
                            const std::string &outputdir) {

    assert(av.size() == ftemplates.size());
    #ifndef NDEBUG
    std::set<std::size_t> m;
    std::size_t n = 0;
    for (std::vector<tom::avg::st_average_result>::const_iterator it=av.begin(); it!=av.end(); it++) {
        n += it->use_idx.size();
        m.insert(it->use_idx.begin(), it->use_idx.end());
    }
    std::set<std::size_t>::const_iterator i=m.end();
    if (!m.empty()) { i--; }
    assert(m.size()==n && (m.empty() || (*i<fparticles.size())));
    #endif


    tom::avg::filename_generator f_gen(outputdir);

    const std::size_t ntemplates = ftemplates.size();
    std::size_t itemplate;
    std::vector<std::size_t>::const_iterator culvit; // Const Unsigned Long Vector ITerator

    clog << "NTEMPLATES " << ntemplates << std::endl;

    // Now compute the average, the FSC, the CC and save it to file.
    for (itemplate=0; itemplate<ntemplates; itemplate++) {

        const tom::avg::st_average_result &ar = av[itemplate];

        const std::vector<std::size_t> &v_ = ar.use_idx;
        clog << "\n"
                "# template " << std::setw(8) << itemplate << " (\"" << ftemplates[itemplate] << "\") has the maximum correlation in " << v_.size() << " particles." << std::endl;


        for (culvit=v_.begin(); culvit!=v_.end(); culvit++) {
            clog << "#    particle " << std::setw(5) << *culvit << " (\"" << fparticles[*culvit] << "\")." << std::endl;
        }

        clog << "# avg_s:          \"" << f_gen.get_avg_s(itemplate) << "\"" << std::endl;
        clog << "# avgwedge_s:     \"" << f_gen.get_avgwedge_s(itemplate) << "\"" << std::endl;
        clog << "# avg_1:          \"" << f_gen.get_avg_1(itemplate) << "\"" << std::endl;
        clog << "# avgwedge_1:     \"" << f_gen.get_avgwedge_1(itemplate) << "\"" << std::endl;
        clog << "# avg_sqrtn:      \"" << f_gen.get_avg_sqrtn(itemplate) << "\"" << std::endl;
        clog << "# avgwedge_sqrtn: \"" << f_gen.get_avgwedge_sqrtn(itemplate) << "\"" << std::endl;
        clog << "# avg_n:          \"" << f_gen.get_avg_n(itemplate) << "\"" << std::endl;
        clog << "# avgwedge_n:     \"" << f_gen.get_avgwedge_n(itemplate) << "\"" << std::endl;

        clog << "NUM " << std::setw(3) << itemplate << "  " << std::setw(5) << v_.size();
        if (v_.empty()) {
            clog << "\n# WARNING: no particles have peak for template #" << itemplate << ". Skip averaging and correlation." << std::endl;
            continue;
        }

        clog << "       ";
        for (culvit=v_.begin(); culvit!=v_.end(); culvit++) {
            clog << ' ' << *culvit;
        }
        clog << std::endl;
        if (!ar.fsc.empty()) {
            clog << "FSC " << std::setw(3) << itemplate << "  ";
            const std::size_t n_shells = ar.fsc.size();
            for (std::size_t i=0; i<n_shells; i++) {
                clog << ar.fsc[i].first << " " << ar.fsc[i].second << "   ";
            }
            clog << std::endl;
        }

        if (ar.ccval_s.get()) {
            clog << "CC_S " << std::setw(3) << itemplate << "  " << std::setw(20) << *ar.ccval_s << std::endl;
        }
        if (ar.ccval_1.get()) {
            clog << "CC_1  " << std::setw(3) << itemplate << "  " << std::setw(20) << *ar.ccval_1  << std::endl;
        }
        if (ar.ccval_sqrtn.get()) {
            clog << "CC_SQRTN  " << std::setw(3) << itemplate << "  " << std::setw(20) << *ar.ccval_sqrtn  << std::endl;
        }
        if (ar.ccval_n.get()) {
            clog << "CC_N " << std::setw(3) << itemplate << "  " << std::setw(20) << *ar.ccval_n << std::endl;
        }
    }
}












template void tom::corr::parse_filelist<float >(const std::string &filename, std::vector<std::string> &ffilenames, std::vector<boost::shared_ptr<tom::WedgeDescriptor<float > > > &wedge_desc, std::size_t binning);
template void tom::corr::parse_filelist<double>(const std::string &filename, std::vector<std::string> &ffilenames, std::vector<boost::shared_ptr<tom::WedgeDescriptor<double> > > &wedge_desc, std::size_t binning);

template void tom::corr::parse_peaklist<float >(const std::string &peakfilename, std::vector<tom::cc_peak<float > > &peak_list, std::vector<helper::triple<double, double, double> > &angles, std::size_t ntemplates, std::size_t nparticles);
template void tom::corr::parse_peaklist<double>(const std::string &peakfilename, std::vector<tom::cc_peak<double> > &peak_list, std::vector<helper::triple<double, double, double> > &angles, std::size_t ntemplates, std::size_t nparticles);


template void tom::corr::write_filelist<double>(const std::string &filename, const std::vector<std::string> &ffilenames, const std::vector<boost::shared_ptr<tom::WedgeDescriptor<double> > > &wedge_desc);

