/****************************************************************************//**
 * \file volume_fcn.hpp
 * \brief The header file for single functions manipulating tom::Volume.
 * \author  Thomas Haller
 * \version 0.1
 * \date    10.12.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__VOLUME_FCN_HPP__
#define ___INCLUDE_CORE__VOLUME_FCN_HPP__



#include <complex>


#include "tom/core/volume.hpp"
#include "helper/triple.hpp"


namespace tom {

template<typename T> bool is_double();
template<typename T> bool is_float();
template<typename T> bool is_float_complex();
template<typename T> bool is_double_complex();


template <typename T> int get_tom_io_type();


typedef helper::triple<std::size_t, std::size_t, std::size_t> st_idx;


template<typename T> void make_binary(tom::Volume<T> &v);


template<typename T> void init_spheremask(tom::Volume<T> &mask, float radius, float sigma, float max_radius);
template<typename T> void init_spheremask(tom::Volume<T> &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);


template<typename T> void hermitian_symmetry_to_full(      tom::Volume<T> &v                      );
template<typename T> void hermitian_symmetry_to_full(const tom::Volume<T> &vsrc, tom::Volume<T> &v);

template<typename T> void fftshift(      tom::Volume<T> &v   ,                       bool is_ifftshift);
template<typename T> void fftshift(const tom::Volume<T> &vsrc, tom::Volume<T> &v   , bool is_ifftshift);

//template<typename T> void shift(tom::Volume<std::complex<T> > &v, std::size_t sizez, double shiftx, double shifty, double shiftz);

namespace norm {
    enum ntype {
        NORM_NO_NORM,
        NORM_STDDEV_1,
        NORM_STDDEV_SAMPLE_1
    };
}
template<typename T, typename T2, typename T3>
void norm_mask(tom::Volume<T> &v, const tom::Volume<T2> &mask, tom::norm::ntype stddev_type, T3 *variance, bool is_boolean_mask);

template<typename T, typename T2>
void update_correlation_volume(const Volume<T> &v, Volume<T> &cc, Volume<T2> &vindex, T2 index);


template<typename T>              std::vector<st_idx> peak(const tom::Volume<T> &v);
template<typename T, typename T2> std::vector<st_idx> peak(const tom::Volume<T> &v, const tom::Volume<T2> &mask);



template<typename T> void rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz);
template<typename T> void rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post);
template<typename T> void rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta);


template<typename T, typename T2> void element_wise_add          (tom::Volume<T> &b, const tom::Volume<T2> &a);
template<typename T, typename T2> void element_wise_sub          (tom::Volume<T> &b, const tom::Volume<T2> &a);
template<typename T, typename T2> void element_wise_multiply     (tom::Volume<T> &b, const tom::Volume<T2> &v);
template<typename T, typename T2> void element_wise_div          (tom::Volume<T> &b, const tom::Volume<T2> &a, T inf_value);
template<typename T, typename T2> void element_wise_conj_multiply(tom::Volume<T> &b, const tom::Volume<T2> &a);
template<typename T1, typename T2, typename T3> void element_wise_conj_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2, const tom::Volume<T3> &v3);

template<typename T> void element_wise_set_below_threshold(tom::Volume<T> &a, T threshold, T value);


template<typename T> void element_wise_power(tom::Volume<T> &v, T exponent);
template<typename T> void element_wise_max(tom::Volume<T> &v, T val);


template<typename T> void fourier_shell_correlation(const tom::Volume<std::complex<T> > &F1, bool is_fft_shifted1, const tom::Volume<std::complex<T> > &F2, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &r);


template<typename T, typename TOP> void element_wise_operation(tom::Volume<T> &v, TOP o);


namespace math {
template<typename T> T log (T x);
template<typename T> T exp (T x);
template<typename T> T sqrt(T x);
template<typename T> T acos(T x);
template<typename T> T abs (T x);
template<typename T> T ceil(T x);
}


}






// Inline functions...

#include <tom/core/volume_loop.hpp>


namespace tom {
namespace fcn {
template<typename T, typename TOP>
struct st_element_wise_operation {
    TOP o_;
    st_element_wise_operation(TOP o): o_(o) { }
    void operator()(T &v, std::size_t, std::size_t, std::size_t) { v = o_(v); }
};
} // namespace tom::fcn
} // namespace tom
template<typename T, typename TOP>
inline void tom::element_wise_operation(tom::Volume<T> &v, TOP o) {
    tom::loop::for_each<T, tom::Volume<T> &, tom::fcn::st_element_wise_operation<T, TOP> >(v, tom::fcn::st_element_wise_operation<T, TOP>(o));
}








namespace tom {
namespace math {
template<          > inline float  log (float  x) { return ::logf(x); }
template<          > inline double log (double x) { return std::log(x); }
template<          > inline float  exp (float  x) { return ::expf(x); }
template<          > inline double exp (double x) { return std::exp(x); }
template<          > inline float  sqrt(float  x) { return ::sqrtf(x); }
template<          > inline double sqrt(double x) { return std::sqrt(x);  }
template<          > inline float  acos(float  x) { return ::acosf(x); }
template<          > inline double acos(double x) { return std::acos(x); }
template<          > inline float  ceil(float  x) { return ::ceilf(x); }
template<          > inline double ceil(double x) { return std::ceil(x); }

template<          > inline long int abs (long int x) { return std::labs(x); }
template<          > inline int      abs (int      x) { return std::abs(x); }
template<          > inline float    abs (float    x) { return ::fabsf(x); }
template<          > inline double   abs (double   x) { return std::fabs(x); }
}
}


namespace tom {
/** \brief Returns true if the template function  is for float */
template <typename T> inline bool is_float       ()         { return false; }
template <          > inline bool is_float<float>()         { return true ; }
/** \brief Returns true if the template function  is for double */
template <typename T> inline bool is_double        ()        { return false; }
template <          > inline bool is_double<double>()        { return true ; }
/** \brief Returns true if the template function is for float complex numbers */
template <typename T> inline bool is_float_complex                      () { return false; }
template <          > inline bool is_float_complex<std::complex<float> >() { return true ; }
/** \brief Returns true if the template function  is for double complex numbers */
template <typename T> inline bool is_double_complex                       ()  { return false; }
template <          > inline bool is_double_complex<std::complex<double> >()  { return true ; }
}


template<typename T>
void tom::fftshift(tom::Volume<T> &v, bool is_ifftshift) {
    tom::Volume<T> v2(v);
    tom::fftshift(v2, v, is_ifftshift);
}



/****************************************************************************//**
 * \brief initialises volume with sphere mask.
 * \param[out] mask
 * \param[in] radius
 * \param[in] sigma
 * \param[in] max_radius
 *
 * Calls init_spheremask with the center parameters set to the center of the
 * volume. Beware, that the center is always integer and is rounded towards
 * infinity, in case of odd mask size.\n
 *******************************************************************************/
template<typename T>
inline void tom::init_spheremask(tom::Volume<T> &mask, float radius, float sigma, float max_radius) {
    ::tom::init_spheremask<T>(mask, radius, sigma, max_radius, mask.getSizeX() / 2, mask.getSizeY() / 2, mask.getSizeZ() / 2);
}


/****************************************************************************//**
 * \brief Rotates the volume around the its center.
 *
 * Beware in case of volume with even size, the center is rounded towards to
 * next integer.
 *******************************************************************************/
template <typename T>
inline void tom::rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta) {
    tom::rotate(src, v_rot, phi, psi, theta, ceil((src.getSizeX()-1)/2.), ceil((src.getSizeY()-1)/2.), ceil((src.getSizeZ()-1)/2.));
}


/****************************************************************************//**
 * \brief Rotates the volume.
 *******************************************************************************/
template<typename T>
void tom::rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz) {
    tom::rotate(src, v_rot, phi, psi, theta, centerx, centery, centerz, 0,0,0, 0,0,0 );
}





#endif

