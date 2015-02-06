/****************************************************************************//**
 * \file fftw_plan.hpp
 * \brief Headerfile with inline-replacements of the fftw-functions as templates.
 * \author  Thomas Haller
 * \version 0.1
 * \date    28.11.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__FFTW_PLAN_HPP__
#define ___INCLUDE_CORE__FFTW_PLAN_HPP__

#include <fftw3.h>


#include "tom/core/volume.hpp"
#include <boost/shared_ptr.hpp>

#define TOM_FFTW_PLAN_INVALID   0
#define TOM_FFTW_PLAN_3D        1
#define TOM_FFTW_PLAN_2D        2


namespace tom {

namespace fftw {



template<typename T> void setFcnMemFFTW();
template<typename T> void *(*fftw_malloc())(size_t n);
template<typename T> void (*fftw_free())(void *ptr);
const char *get_wisdom_name();
int         set_wisdom_name(const char *fname, int load);
void        clear_wisdom_name();
int         save_wisdom();




template<typename T>
class Plan {
public:
    typedef enum { FFTW_DFT_R2C, FFTW_DFT_C2R, FFTW_DFT_FORWARD, FFTW_DFT_BACKWARD } fftw_plan_type;

	Plan(tom::fftw::Plan<T> &originalPlan);
    Plan(tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst, unsigned flags);
    Plan(tom::Volume<std::complex<T> > &vsrc, tom::Volume<T> &vdst, unsigned flags);
    ~Plan();

    static int is_valid_dft_r2c(const tom::Volume<T> &vsrc, const tom::Volume<std::complex<T> > &vdst);
    static int is_valid_dft_c2r(const tom::Volume<std::complex<T> > &vsrc, const tom::Volume<T> &vdst);



    void execute(tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst) const;
    void execute(tom::Volume<std::complex<T> > &vsrc, tom::Volume<T> &vdst) const;

private:
    fftw_plan_type type;

    tom::Volume<std::complex<T> > *p_freq0, *p_freq1;
    tom::Volume<T> *p_time0, *p_time1;

	boost::shared_ptr<void> plan;
};


std::string flag2str(unsigned flag);
unsigned str2flag(const std::string &s);

} // namespace fftw
} // namespace tom





template<typename T>
inline void tom::fftw::setFcnMemFFTW() {
    tom::setFcnMem(::tom::fftw::fftw_malloc<T>(), ::tom::fftw::fftw_free<T>());
}

namespace tom {
namespace fftw {
template<> inline void *(*fftw_malloc<float >())(size_t n) { return &::fftwf_malloc; }
template<> inline void *(*fftw_malloc<double>())(size_t n) { return &::fftw_malloc ; }
template<> inline void (*fftw_free<float >())(void *ptr) { return &::fftwf_free; }
template<> inline void (*fftw_free<double>())(void *ptr) { return &::fftw_free ; }
}
}






#endif




