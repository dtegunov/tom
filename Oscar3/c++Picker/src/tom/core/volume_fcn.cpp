/***********************************************************************//**
 * \file volume_fcn.cpp
 * \brief Implementation of functions manipulating the volume class.
 * \author  Thomas Haller
 * \version 0.1
 * \date    22.11.2007
 **************************************************************************/
#include "tom/core/volume_fcn.hpp"



#include <assert.h>
#include <fftw3.h>
#include <limits>

#include <tom/core/volume_loop.hpp>
#include <tom/core/io.h>
#include <tom/core/transform.h>
#include <tom/core/transform.hpp>






namespace tom {
/****************************************************************************//**
 * \brief Returns an integer specifying the data-type.
 *
 * \returns Number as specified in \link tom_defines.h \endlink for defines
 *   such as TOM_IO_TYPE_DOUBLE. Returns TOM_IO_TYPE_UNKNOWN
 *   if the type is not defined there.
 *******************************************************************************/
template<typename T> int get_tom_io_type                      () { return TOM_IO_TYPE_UNKNOWN; }
template<          > int get_tom_io_type<char                >() { return TOM_IO_TYPE_INT8; }
template<          > int get_tom_io_type<int8_t              >() { return TOM_IO_TYPE_INT8; }
template<          > int get_tom_io_type<int16_t             >() { return TOM_IO_TYPE_INT16; }
template<          > int get_tom_io_type<int32_t             >() { return TOM_IO_TYPE_INT32; }
template<          > int get_tom_io_type<float               >() { return TOM_IO_TYPE_FLOAT; }
template<          > int get_tom_io_type<fftwf_complex       >() { return TOM_IO_TYPE_COMPLEX4; }
template<          > int get_tom_io_type<std::complex<float> >() { return TOM_IO_TYPE_COMPLEX4; }
template<          > int get_tom_io_type<double              >() { return TOM_IO_TYPE_DOUBLE; }
template<          > int get_tom_io_type<int64_t             >() { return TOM_IO_TYPE_INT64; }
}







namespace {
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__get_mean {
    for_each__tom__norm_mask__get_mean(): n(0), sum(0) { }
    std::size_t n; TPRECISION sum;
    inline void operator()(const T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            this->n ++;
            this->sum += v;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__apply_mask_no_variance {
    for_each__tom__norm_mask__apply_mask_no_variance(TPRECISION mean): sum(0), mean(mean) { }
    TPRECISION sum;
    TPRECISION const mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            const TPRECISION v2 = (v - this->mean) * mask;
            v = v2;
            this->sum += v2;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__set_mean0 {
    for_each__tom__norm_mask__set_mean0(TPRECISION mean): mean(mean) { }
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            v = v - this->mean;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__set_mean0_all {
    for_each__tom__norm_mask__set_mean0_all(TPRECISION mean): mean(mean) { }
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            v = v - this->mean;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__sub_meanv {
    for_each__tom__norm_mask__sub_meanv(TPRECISION mean): sum(0), sum2(0), mean(mean) { }
    TPRECISION sum, sum2;
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            const TPRECISION v2 = (v - this->mean) * mask;
            v = v2;
            this->sum += v2;
            this->sum2 += v2*v2;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__sub_meanv_bool {
    for_each__tom__norm_mask__sub_meanv_bool(TPRECISION mean): sum(0), sum2(0), mean(mean) { }
    TPRECISION sum, sum2;
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            const TPRECISION v2 = v - this->mean;
            v = v2;
            this->sum += v2;
            this->sum2 += v2*v2;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_under_mask_5 {
    for_each__tom__norm_under_mask_5(TPRECISION mean, TPRECISION stddev): mean(mean), stddev(stddev) { assert(stddev != 0); }
    const TPRECISION mean, stddev;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            v = (v - this->mean) / this->stddev;
        }
    }
};
}



/****************************************************************************//**
 * \brief Normalises the input volume to mean 0 and standard deviation 1 under mask
 *
 * \param[in,out] v The volume to be normalised.
 * \param[in] mask A mask under which \a v is to be normalized. \a mask must have the same
 *   size as the volume \a v. The mask should contain numbers between 0 and 1.
 * \param[in] stddev_type How to finally scale the volume. It can be one of the integer
 *   constants NORM_NO_NORM, NORM_STDDEV_1 and NORM_STDDEV_SAMPLE_1.
 * \param[out] variance If not NULL, the variance of the volume is returned.
 *   If \a stddev_type == NORM_NO_NORM then it is the variance
 *   of the volume after executing the function (with divisor N, not (N-1)). If \a stddev_type is
 *   NORM_STDDEV_1 or NORM_STDDEV_SAMPLE_1 it always 1.
 *
 * All computations are done with the precision of the template type \a T3. \n
 * What this function does is the following: first the volume \a v is scaled to have
 * a mean of 0 under the \a mask (i.e. the voxels of \a v where the corresponding
 * voxel of the \a mask is != 0). Then each element of the \a v is multiplied with
 * \a mask and the volume is again shifted to have a mean of 0 under the mask.
 * If \a stddev_type == NORM_NO_NORM
 * then the variance is computed and the function exits. Otherwise a final scaling is done
 * to set the standard deviation. If NORM_STDDEV_1 it takes the standard
 * deviation with divisor N (N the number of voxels). If NORM_STDDEV_SAMPLE_1
 * it takes the sample standard deviation with divisor (N-1).
 *******************************************************************************/
template<typename T, typename TMASK, typename TPRECISION>
void tom::norm_mask(tom::Volume<T> &v, const tom::Volume<TMASK> &mask, tom::norm::ntype stddev_type, TPRECISION *variance, bool is_boolean_mask) {

    if (stddev_type != tom::norm::NORM_NO_NORM &&
        stddev_type != tom::norm::NORM_STDDEV_1 &&
        stddev_type != tom::norm::NORM_STDDEV_SAMPLE_1) {
        throw std::invalid_argument("The normalisation type of the deviation is not defined.");
    }


    std::size_t val_n;
    const std::size_t numel = v.getSizeX() * v.getSizeY() * v.getSizeZ();
    TPRECISION val_sum, val_sum2, val_mean, val_stddev, val_variance;

    val_n = 0;
    val_sum = 0.;
    val_sum2 = 0.;
    {
        // First compute the number of elements with a mask-value ~= 0 and its mean.
        for_each__tom__norm_mask__get_mean<T, TMASK, TPRECISION> s;
        tom::loop::for_each<const T, const tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__get_mean<T, TMASK, TPRECISION> &>(v, mask, s);
        val_n = s.n;
        val_sum = s.sum;
    }
    val_mean = val_sum / static_cast<TPRECISION>(val_n);

    if (stddev_type == tom::norm::NORM_NO_NORM && !variance) {
        // First compute the number of elements with a mask-value ~= 0 and its mean.
        if (is_boolean_mask) {
            if (val_mean != 0.) {
                // Finally normalise again under the mask.
                tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__set_mean0_all<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0_all<T, TMASK, TPRECISION>(val_mean));
            }
        } else {
            for_each__tom__norm_mask__apply_mask_no_variance<T, TMASK, TPRECISION> s(val_mean);
            tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__apply_mask_no_variance<T, TMASK, TPRECISION> &>(v, mask, s);
            val_sum = s.sum;
            val_mean = val_sum / static_cast<TPRECISION>(val_n);
            if (val_mean != 0.) {
                // Finally normalise again under the mask.
                tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION>(val_mean));
            }
        }
    } else {
        // Shift inside the mask to give mean==0 and multiply with the mask.
        if (is_boolean_mask) {
            // First compute the number of elements with a mask-value ~= 0 and its mean.
            for_each__tom__norm_mask__sub_meanv_bool<T, TMASK, TPRECISION> s(val_mean);
            tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__sub_meanv_bool<T, TMASK, TPRECISION> &>(v, mask, s);
            val_sum = s.sum;
            val_sum2 = s.sum2;
        } else {
            // First compute the number of elements with a mask-value ~= 0 and its mean.
            for_each__tom__norm_mask__sub_meanv<T, TMASK, TPRECISION> s(val_mean);
            tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__sub_meanv<T, TMASK, TPRECISION> &>(v, mask, s);
            val_sum = s.sum;
            val_sum2 = s.sum2;
        }

        val_mean = val_sum / static_cast<TPRECISION>(val_n);
        val_variance = (val_sum2 - static_cast<TPRECISION>(val_n)*val_mean*val_mean) / static_cast<TPRECISION>(numel - (stddev_type==tom::norm::NORM_STDDEV_SAMPLE_1 ? 1 : 0));

        if (stddev_type != tom::norm::NORM_NO_NORM) {
            val_stddev = sqrt(val_variance);
            if (val_mean!=0. || (val_stddev!=1.&&val_stddev!=0.)) {
                // Finally normalise again under the mask.
                if (val_stddev!=1. && val_stddev!=0.) {
                    tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_under_mask_5<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_under_mask_5<T, TMASK, TPRECISION>(val_mean, val_stddev));
                } else {
                    tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION>(val_mean));
                }
            }
            if (val_stddev) {
                val_variance = 1.;
            } else {
                val_variance = 0.;
            }
        } else {
            if (val_mean!=0.) {
                // Finally normalise again under the mask.
                tom::loop::for_each<T, tom::Volume<T> &, const TMASK, const tom::Volume<TMASK> &, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION>(val_mean));
            }
        }
        if (variance) {
            *variance = val_variance;
        }
    }
}





namespace {
template<typename T1, typename T2>
struct for_each__tom__element_wise_add {
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { v1 = v1 + v2; }
};
}
/****************************************************************************//**
 * \brief Adds two volumes elementwise.
 *
 * \param[in,out] v1 The left side operand of the elementwise addition and the
 *   destination volume at the same time.
 * \param[in] v2 The right side operand of the addition.
 *
 * Self assignment is not a problem. In that case the volume is added to itself.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_add(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each<T1, tom::Volume<T1> &, const T2, const tom::Volume<T2> &, ::for_each__tom__element_wise_add<T1, T2> >(v1, v2, ::for_each__tom__element_wise_add<T1, T2>());
}





namespace {
template<typename T1, typename T2>
struct for_each__tom__element_wise_sub {
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { v1 = v1 - v2; }
};
}
/****************************************************************************//**
 * \brief Subtracts two volumes elementwise.
 *
 * \param[in,out] v1 The left side operand of the elementwise subtraction and the
 *   destination volume at the same time.
 * \param[in] v2 The right side operand of the subtraction.
 *
 * Self assignment is not a problem. In that case the volume is subtracted by itself. Yields 0 volume.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_sub(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each<T1, tom::Volume<T1> &, const T2, const tom::Volume<T2> &, ::for_each__tom__element_wise_sub<T1, T2> >(v1, v2, ::for_each__tom__element_wise_sub<T1, T2>());
}




namespace {
template<typename T1, typename T2>
struct for_each__tom__element_wise_div {
    for_each__tom__element_wise_div(T1 inf_value): inf_value(inf_value) { }
    T1 inf_value;
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) {
        if (v2) {
            v1 = v1 / v2;
        } else {
            v1 = inf_value;
        }
    }
};
}
/****************************************************************************//**
 * \brief Elementwise division
 *
 * \param[in,out] b The left side operand of the elementwise division (numerator)
 *   and the destination volume at the same time.
 * \param[in] a The right side operand of the addition (denominator)
 * \param[in] inf_value Value to assign if the corresponding voxel of @a a is 0.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_div(tom::Volume<T1> &v1, const tom::Volume<T2> &v2, T1 inf_value) {
    tom::loop::for_each<T1, tom::Volume<T1> &, const T2, const tom::Volume<T2> &, ::for_each__tom__element_wise_div<T1, T2> >(v1, v2, ::for_each__tom__element_wise_div<T1, T2>(inf_value));
}




namespace {
template<typename T>
struct for_each__element_wise_set_below_threshold {
    for_each__element_wise_set_below_threshold(T threshold, T value): value(value), threshold(threshold) { }
    T value, threshold;
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) const {
        if (a < threshold) { a = value; }
    }
}; // struct for_each__element_wise_div_threshold
} // namespace
/****************************************************************************//**
 * \brief Set lower limit of the volume.
 *
 * \param[in,out] a The volume to be limited.
 * \param[in] threshold The threshold
 * \param[in] value The value.
 *
 * Every value of the volume \a a which is smaller! then threshold is set
 * to value.
 *******************************************************************************/
template<typename T>
void tom::element_wise_set_below_threshold(tom::Volume<T> &a, T threshold, T value) {
    tom::loop::for_each<T, tom::Volume<T> &, const ::for_each__element_wise_set_below_threshold<T> &>(a, ::for_each__element_wise_set_below_threshold<T>(threshold, value));
}






namespace {
template<typename T1, typename T2> inline void for_each__tom__element_wise_multiply_(T1 &v1, const T2 &v2) { v1 = v1 * v2; }
#if defined __xlC__
// in the complex header from the currently used xlC compiler the methods real and imag to not return a reference...
// strange :)
template<> inline void for_each__tom__element_wise_multiply_<std::complex<float >, double>(std::complex<float > &v1, const double &v2) { v1.real(v1.real()*v2); v1.imag(v1.imag()*v2); }
template<> inline void for_each__tom__element_wise_multiply_<std::complex<double>, float >(std::complex<double> &v1, const float  &v2) { v1.real(v1.real()*v2); v1.imag(v1.imag()*v2); }
#else
template<> inline void for_each__tom__element_wise_multiply_<std::complex<float >, double>(std::complex<float > &v1, const double &v2) { v1.real() = v1.real()*v2; v1.imag() = v1.imag()*v2; }
template<> inline void for_each__tom__element_wise_multiply_<std::complex<double>, float >(std::complex<double> &v1, const float  &v2) { v1.real() = v1.real()*v2; v1.imag() = v1.imag()*v2; }
#endif
template<typename T1, typename T2>
struct for_each__tom__element_wise_multiply {
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { for_each__tom__element_wise_multiply_<T1, T2>(v1, v2); }
};
}
/****************************************************************************//**
 * \brief Elementwise multiplication
 *
 * \param[in,out] v1 The left side operand of the elementwise multiplication
 *   and the destination volume at the same time.
 * \param[in] v2 The right side operand of the multiplication.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each<T1, tom::Volume<T1> &, const T2, const tom::Volume<T2> &, ::for_each__tom__element_wise_multiply<T1, T2> >(v1, v2, ::for_each__tom__element_wise_multiply<T1, T2>());
}





namespace {
template<typename T1, typename T2>
struct for_each__tom__element_wise_conj_multiply {
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { v1 = conj(v1) * v2; }
};
}
/****************************************************************************//**
 * \brief Elementwise multiplication
 *
 * \param[in,out] v1 The left side operand of the elementwise multiplication
 *   and the destination volume at the same time. The conjugate complex
 *   value of b is taken.
 * \param[in] v2 The right side operand of the multiplication.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_conj_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each<T1, tom::Volume<T1> &, const T2, const tom::Volume<T2> &, ::for_each__tom__element_wise_conj_multiply<T1, T2> >(v1, v2, ::for_each__tom__element_wise_conj_multiply<T1, T2>());
}



namespace {
template<typename T1, typename T2, typename T3>
struct for_each__tom__element_wise_conj_multiply3 {
    inline void operator()(T1 &v1, const T2 &v2, const T3 &v3, std::size_t, std::size_t, std::size_t) { v1 = conj(v2) * v3; }
};
}
/****************************************************************************//**
 * \brief Elementwise multiplication
 *
 * \param[in,out] v1 The left side operand of the elementwise multiplication
 *   and the destination volume at the same time. The conjugate complex
 *   value of b is taken.
 * \param[in] v2 The right side operand of the multiplication.
 *******************************************************************************/
template<typename T1, typename T2, typename T3>
void tom::element_wise_conj_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2, const tom::Volume<T3> &v3) {
    tom::loop::for_each<T1, tom::Volume<T1> &, const T2, const tom::Volume<T2> &, const T3, const tom::Volume<T3> &, ::for_each__tom__element_wise_conj_multiply3<T1, T2, T3> >(v1, v2, v3, ::for_each__tom__element_wise_conj_multiply3<T1, T2, T3>());
}




/****************************************************************************//**
 * \brief Shift zero-frequency component to center of spectrum.
 *
 * \param[in] vsrc The source volume.
 * \param[out] v The shifted volume.
 * \param[in] is_ifftshift If true it behaves as \c ifftshift from MATLAB.
 *   Otherwise as \c fftshift. In case of even volume sizes, there is not
 *   difference. Otherwise ifftshift is the inverse operation of the fftshift.
 *
 * See the documentation of fftshift/ifftshift from MATLAB :)
 * \todo The current implementation could be done better.
 *******************************************************************************/
template<typename T>
void tom::fftshift(const Volume<T> &vsrc, Volume<T> &v, bool is_ifftshift) {

    if (vsrc.getSizeX() != v.getSizeX() || vsrc.getSizeY() != v.getSizeY() || vsrc.getSizeZ() != v.getSizeZ()) {
        throw std::runtime_error("The volumes must have the same size.");
    }

    const std::size_t sizex = v.getSizeX();
    const std::size_t sizey = v.getSizeY();
    const std::size_t sizez = v.getSizeZ();
    const std::size_t sizex2 = sizex/2 + (is_ifftshift ? 0 : sizex%2);
    const std::size_t sizey2 = sizey/2 + (is_ifftshift ? 0 : sizey%2);
    const std::size_t sizez2 = sizez/2 + (is_ifftshift ? 0 : sizez%2);

    std::size_t x, y, z;
    std::size_t z_tmp, y_tmp, x_tmp;
    if (v.isContiguous() && vsrc.isContiguous()) {
        const T *pvsrc = &vsrc.get();
        T *pv = &v.get();
        std::size_t i;
        i = 0;
        for (z=0; z<sizez; z++) {
            z_tmp = ((z+sizez2)%sizez) * sizey;
            for (y=0; y<sizey; y++) {
                y_tmp = (z_tmp + (y+sizey2)%sizey) * sizex;
                for (x=0; x<sizex; x++) {
                    x_tmp = y_tmp + (x+sizex2)%sizex;
                    pv[i++] = pvsrc[x_tmp];
                }
            }
        }
    } else {
        for (z=0; z<sizez; z++) {
            z_tmp = (z+sizez2)%sizez;
            for (y=0; y<sizey; y++) {
                y_tmp = (y+sizey2)%sizey;
                for (x=0; x<sizex; x++) {
                    v.get(x,y,z) = vsrc.get((x+sizex2)%sizex, y_tmp, z_tmp);
                }
            }
        }
    }
}







namespace {
template<typename T1, typename T2>
struct for_each__tom__peakv {
    for_each__tom__peakv(T1 max): max(max) { }
    T1 max;
    std::vector<tom::st_idx> res;
    inline void operator()(const T1 &v1, const T2 &v2, const std::size_t &x, const std::size_t &y, const std::size_t &z) {
        if (v2 && !(v1<this->max)) {
            tom::st_idx idx;
            idx.x = x;
            idx.y = y;
            idx.z = z;
            if (v1 > this->max) {
                this->max = v1;
                this->res.clear();
            }
            res.push_back(idx);
        }
    }
};
}
/****************************************************************************//**
 * Gets the index of the maximum value.
 *******************************************************************************/
template<typename T1, typename T2>
std::vector<tom::st_idx> tom::peak(const tom::Volume<T1> &v, const tom::Volume<T2> &mask) {

    for_each__tom__peakv<T1, T2> s(- std::numeric_limits<T1>::max());

    tom::loop::for_each<const T1, const tom::Volume<T1> &, const T2, const tom::Volume<T2> &, ::for_each__tom__peakv<T1, T2> &>(v, mask, s);

    return s.res;
}

/****************************************************************************//**
 * \brief Expands the reduced complex volume to its full dimension.
 * \param[in] vsrc The complex volume. It can be either reduced
 *   (with dimension z only of size z/2+1) or full. However only
 *   the upper half is considered.
 * \param[out] v The complex volume. It has the dimension of the full
 *   (output) volume. The upper half (0:sizex-1, 0:sizey-1, 0:sizez/2+1)
 *   contains the complex values as returned by a R2C fftw.
 *   The second half (0:sizex-1, 0:sizey-1, sizez/2+2:sizez) i filled with
 *   the conjugate complex values to fullfill the hermitian symmetrie.
 *******************************************************************************/
template<typename T>
void tom::hermitian_symmetry_to_full(const tom::Volume<T> &vsrc, tom::Volume<T> &v) {
    if (vsrc.getSizeZ() != v.getSizeZ() && vsrc.getSizeZ() != v.getSizeZ()/2+1) {
        throw std::invalid_argument("The input volume must be either reduced or full complex.");
    }
    std::size_t sizez = v.getSizeZ()/2+1;
    if (&vsrc != &v) {
        T *data = const_cast<T *>(&vsrc.get());
        tom::Volume<T>(v, NULL, v.getSizeX(), v.getSizeY(), sizez, v.getStrideX(), v.getStrideY(), v.getStrideZ()).setValues(
            tom::Volume<T>(data, vsrc.getSizeX(), vsrc.getSizeY(), sizez, vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ(), false, NULL));
    }
    tom::hermitian_symmetry_to_full(v);
}






namespace {
template<typename T>
inline T local_fcn_conjugate_complex(const T &a) {
    return conj(a);
}
template<> inline float  local_fcn_conjugate_complex(const float  &a) { return a; }
template<> inline double local_fcn_conjugate_complex(const double &a) { return a; }
}
/****************************************************************************//**
 * \brief Expands the reduced volume to its full dimension.
 *
 * \param[in,out] v The volume. It has the dimension of the full
 *   (output) volume. The upper half (0:sizex-1, 0:sizey-1, 0:sizez/2+1)
 *   contains the already initialized values.
 *   This function copies that half into the second half accoding to
 *   the properties of the hermitian symetrie.\n
 *   For example the fftw real to complex transformations return a the half
 *   volume of complex numbers.
 *
 * Can also be a pure real volume (T==float or double)
 *******************************************************************************/
template<typename T>
void tom::hermitian_symmetry_to_full(tom::Volume<T> &v) {
    typedef T local_TYPE;
    local_TYPE *local_A = &v.get();
    local_TYPE *v_data = &v.get();
    std::size_t local_sizex = v.getSizeX();
    std::size_t local_sizey = v.getSizeY();
    std::size_t local_sizez = v.getSizeZ();
    std::size_t local_stridex = v.getStrideX();
    std::size_t local_stridey = v.getStrideY();
    std::size_t local_stridez = v.getStrideZ();

    {
        std::size_t x, y, z;
        std::size_t x2, y2, z2;
        std::size_t y2_cnt, x2_cnt;

        tom::loop::ptr_add_byte_offset(local_A, (local_sizez/2+1)*local_stridez);
        if (!(local_stridex%sizeof(local_TYPE)) &&
            !(local_stridey%sizeof(local_TYPE)) &&
            !(local_stridez%sizeof(local_TYPE))) {

            local_stridex /= sizeof(local_TYPE);
            local_stridey /= sizeof(local_TYPE);
            local_stridez /= sizeof(local_TYPE);

            if (local_stridex == 1) {
                if (local_stridey == local_sizex &&
                    local_stridez == local_sizex*local_sizey) {
                    std::size_t i = 0;
                    for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                        for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                            y2 = y2_cnt % local_sizey;
                            for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, i++, x2_cnt--) {
                                x2 = x2_cnt % local_sizex;
                                local_A[i] = local_fcn_conjugate_complex<local_TYPE>(v_data[z2*local_stridez + y2*local_stridey + x2]);
                            }
                        }
                    }
                } else {
                    local_stridez -= local_sizey*local_stridey;
                    for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                        for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                            y2 = y2_cnt % local_sizey;
                            for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, x2_cnt--) {
                                x2 = x2_cnt % local_sizex;
                                local_A[x] = local_fcn_conjugate_complex<local_TYPE>(v_data[z2*local_stridez + y2*local_stridey + x2]);
                            }
                            local_A += local_stridey;
                        }
                        local_A += local_stridez;
                    }
                }
            } else {
                local_stridez -= local_sizey*local_stridey;
                local_stridey -= local_sizex*local_stridex;
                for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                    for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                        y2 = y2_cnt % local_sizey;
                        for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, x2_cnt--) {
                            x2 = x2_cnt % local_sizex;
                            local_A[0] = local_fcn_conjugate_complex<local_TYPE>(v_data[z2*local_stridez + y2*local_stridey + x2*local_stridex]);
                            local_A += local_stridex;
                        }
                        local_A += local_stridey;
                    }
                    local_A += local_stridez;
                }
            }
        } else {
            local_stridez -= local_sizey*local_stridey;
            local_stridey -= local_sizex*local_stridex;
            for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                    y2 = y2_cnt % local_sizey;
                    for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, x2_cnt--) {
                        x2 = x2_cnt % local_sizex;
                        *local_A = local_fcn_conjugate_complex<local_TYPE>(* reinterpret_cast<local_TYPE *>((reinterpret_cast<char *>(v_data)) + z2*local_stridez + y2*local_stridey + x2*local_stridex));
                        tom::loop::ptr_add_byte_offset(local_A, local_stridex);
                    }
                    tom::loop::ptr_add_byte_offset(local_A, local_stridey);
                }
                tom::loop::ptr_add_byte_offset(local_A, local_stridez);
            }
        }
    }
}








namespace {
template<typename T, typename TIDX>
inline void for_each__tom__update_correlation_volume_(const T &v, T &vmax, TIDX &vindex, TIDX index) {
    if (vmax < v) {
        vmax = v;
        vindex = index;
    }
}
template<> inline void for_each__tom__update_correlation_volume_<std::complex<float >, int>(const std::complex<float > &v, std::complex<float > &vmax, int &vindex, int index) { throw std::runtime_error("Elementwise maximum of complex type is not defined."); }
template<> inline void for_each__tom__update_correlation_volume_<std::complex<double>, int>(const std::complex<double> &v, std::complex<double> &vmax, int &vindex, int index) { throw std::runtime_error("Elementwise maximum of complex type is not defined."); }
template<typename T, typename TIDX>
struct for_each__tom__update_correlation_volume {
    for_each__tom__update_correlation_volume(TIDX index): index(index) { }
    const TIDX index;
    inline void operator()(const T &v, T &vmax, TIDX & vindex, std::size_t, std::size_t, std::size_t) {
        ::for_each__tom__update_correlation_volume_<T, TIDX>(v, vmax, vindex, this->index);
    }
};
}
/****************************************************************************//**
 * \brief Take the elementwise maximum of two volumes.
 *******************************************************************************/
template<typename T, typename TIDX>
void tom::update_correlation_volume(const tom::Volume<T> &v, tom::Volume<T> &cc, tom::Volume<TIDX> &vindex, TIDX index) {
    tom::loop::for_each<    const T, const tom::Volume<  T> &,
                                  T,       tom::Volume<  T> &,
                               TIDX,       tom::Volume<TIDX> &,
                            ::for_each__tom__update_correlation_volume<T, TIDX> >(
                                v, cc, vindex, ::for_each__tom__update_correlation_volume<T, TIDX>(index));
}















/****************************************************************************//**
 * \brief First rotate the volume around a point and then shifts it.
 *
 * Uses the tom-convection for rotation. It does not rotate the object,
 * but its coordinate system with angles phi, theta, psi (in that order).
 *******************************************************************************/
template<typename T>
void tom::rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post) {

    if (&v_rot==&src) {
        throw std::invalid_argument("Self assignment is not allowed.");
    }
    if (v_rot.getSizeX()!=src.getSizeX() || v_rot.getSizeY()!=src.getSizeY() || v_rot.getSizeZ()!=src.getSizeZ()) {
        throw std::invalid_argument("The volumes must have the same size.");
    }

    if (src.getStrideX()!=sizeof(T) || src.getStrideY()%sizeof(T) || src.getStrideZ()%sizeof(T) ||
        v_rot.getStrideX()!=sizeof(T) || v_rot.getStrideY()%sizeof(T) || v_rot.getStrideZ()%sizeof(T)) {
        throw std::invalid_argument("The current implementation allows only strides of whole multiples of the element size.");
    }
    if (!tom::is_double<T>() && !tom::is_float<T>()) {
        throw std::invalid_argument("Rotation currently is only implemented for float and double type.");
    }

    void (*intfcn)() = NULL;
    double optparam_double = 0.;
    double optparam_float = 0.;
    void *optparam;
    if (!intfcn) {
        /* Set default function. */
        if (tom::is_double<T>()) {
            intfcn = (void (*)())&tom_transf_interpol3d_trilinear_double;
            optparam = &optparam_double;
        } else if (tom::is_float<T>()){
            intfcn = (void (*)())&tom_transf_interpol3d_trilinear_float;
            optparam = &optparam_float;
        } else {
            assert(0);
        }
    }


    double P[16] = { 0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 1. };
    const int axes[3]      = {   2,     0,   2 };
    const double angles[3] = { phi, theta, psi };
    tom_transf_init_rotmat_3d(P, 4, 3, axes, angles);
    tom_transf_invert(P, P, 4);

    P[ 3] = shiftx_post + centerx + P[ 0]*(-centerx+shiftx_pre) + P[ 1]*(-centery+shifty_pre) + P[ 2]*(-centerz+shiftz_pre);
    P[ 7] = shifty_post + centery + P[ 4]*(-centerx+shiftx_pre) + P[ 5]*(-centery+shifty_pre) + P[ 6]*(-centerz+shiftz_pre);
    P[11] = shiftz_post + centerz + P[ 8]*(-centerx+shiftx_pre) + P[ 9]*(-centery+shifty_pre) + P[10]*(-centerz+shiftz_pre);



    #if 0
    int res;
    if (tom::is_double<T>()) {
        res = tom_transf_transf3d_double((const double *)&src.get(), src.getSizeX(), src.getSizeY(), src.getSizeZ(), src.getStrideY()/sizeof(T), src.getStrideZ()/sizeof(T),
                                               (double *)&v_rot.get(), v_rot.getSizeX(), v_rot.getSizeY(), v_rot.getSizeZ(), v_rot.getStrideY()/sizeof(T), v_rot.getStrideZ()/sizeof(T),
                                               P, true, 0., (double (*)(const tom_transf_st_interpol_param*))intfcn, optparam);
        assert(!res);
    } else if (tom::is_float<T>()) {
        res = tom_transf_transf3d_float((const float *)&src.get(), src.getSizeX(), src.getSizeY(), src.getSizeZ(), src.getStrideY()/sizeof(T), src.getStrideZ()/sizeof(T),
                                              (float *)&v_rot.get(), v_rot.getSizeX(), v_rot.getSizeY(), v_rot.getSizeZ(), v_rot.getStrideY()/sizeof(T), v_rot.getStrideZ()/sizeof(T),
                                               P, true, 0., (float (*)(const tom_transf_st_interpol_param*))intfcn, optparam);
        assert(!res);
    } else {
        assert(0);
    }
    #else
    tom::transform<T, tom::InterpolTriLinear<T> >(src, v_rot, P, true, 0, tom::InterpolTriLinear<T>(0));
    #endif

}







namespace {
template<typename T>
struct for_each__tom__make_binary {
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) {
        a = a ? 1 : 0;
    }
};
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::make_binary(tom::Volume<T> &v) {
    tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__make_binary<T> >(v, ::for_each__tom__make_binary<T>());
}









namespace {
template<typename T>
struct for_each__tom__element_wise_power_sqrt {
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { a = tom::math::sqrt<T>(a); }
};
template<typename T>
struct for_each__tom__element_wise_power {
	for_each__tom__element_wise_power(T exponent): exponent(exponent) { assert(exponent>0); }
	T exponent;
	inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { a = tom::math::exp<T>(this->exponent*tom::math::log<T>(a)); }
};
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::element_wise_power(tom::Volume<T> &v, T exponent) {
	if (exponent == 0.5) {
		tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__element_wise_power_sqrt<T> >(v, ::for_each__tom__element_wise_power_sqrt<T>());
	} else {
		if (exponent <= 0) {
			throw std::invalid_argument("The exponent must be greater 0");
		}
		tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__element_wise_power<T> >(v, ::for_each__tom__element_wise_power<T>(exponent));
	}
}

namespace {
template<typename T>
struct for_each__tom__element_wise_max {
	for_each__tom__element_wise_max(T val): val(val) { }
	const T val;
	inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { if (a < this->val) { a = this->val; } }
};
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::element_wise_max(tom::Volume<T> &v, T val) {
	tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__element_wise_max<T> >(v, ::for_each__tom__element_wise_max<T>(val));
}








namespace {
template<typename T>
struct for_each__tom__peak {
    for_each__tom__peak(T maxinit): max(maxinit) { }
    std::vector<tom::st_idx> res;
    T max;
    inline void operator()(const T &a, std::size_t x, std::size_t y, std::size_t z) {
        if (!(a < this->max)) {
            tom::st_idx idx;
            idx.x = x;
            idx.y = y;
            idx.z = z;
            if (a > this->max) {
                this->max = a;
                res.clear();
            }
            this->res.push_back(idx);
        }
    }
};
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::vector<tom::st_idx> tom::peak(const tom::Volume<T> &v) {

    ::for_each__tom__peak<T> p(v.get());

    tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__peak<T> &>(v, p);
    return p.res;
}







namespace {
template<typename TPRECISION, typename T>
struct for_each_step__tom__init_spheremask_sigma {
    for_each_step__tom__init_spheremask_sigma(TPRECISION centerx, TPRECISION centery, TPRECISION centerz,
                                         TPRECISION radius, TPRECISION sigma, TPRECISION max_radius)
                                         : centerx(centerx), centery(centery), centerz(centerz),
                                           radius(radius), sigma(sigma), max_radius(max_radius),
                                           current_z(0), current_y(0) { }
    TPRECISION centerx, centery, centerz, radius, sigma, max_radius;
    TPRECISION current_z, current_y;
    inline void stepz(TPRECISION z) { current_z = z - centerz; current_z *= current_z; }
    inline void stepy(TPRECISION y) { current_y = y - centery; current_y  = current_y*current_y + current_z; }
    inline void operator()(T &a, TPRECISION x) {
        x -= this->centerx;
        x = tom::math::sqrt<TPRECISION>(x*x + this->current_y);
        if (x <= this->radius) {
            a = 1;
        } else if (x > this->max_radius) {
            a = 0;
        } else {
            x = ((x - this->radius)/this->sigma);
            a = tom::math::exp<TPRECISION>(- x*x);
        }
    }
};
template<typename TPRECISION>
struct for_each_step__tom__init_spheremask_sigma<TPRECISION, char> {
    for_each_step__tom__init_spheremask_sigma(TPRECISION centerx, TPRECISION centery, TPRECISION centerz,
                                         TPRECISION radius, TPRECISION sigma, TPRECISION max_radius_input)
                                         : centerx(centerx), centery(centery), centerz(centerz),
                                           radius(radius), sigma(sigma), max_radius(max_radius_input*max_radius_input),
                                           current_z(0), current_y(0) { }
    TPRECISION centerx, centery, centerz, radius, sigma, max_radius;
    TPRECISION current_z, current_y;
    inline void stepz(TPRECISION z) {
        this->current_z = z - this->centerz;
        this->current_z *= this->current_z;
    }
    inline void stepy(TPRECISION y) {
        this->current_y = y - this->centery;
        this->current_y = this->current_y*this->current_y + this->current_z;
    }
    inline void operator()(char &a, TPRECISION x) {
        x -= this->centerx;
        // Sigma with an integer type is only the same as using max_radius.
        a = x*x + this->current_y <= this->max_radius;
    }
};
template<typename TPRECISION, typename T>
struct for_each_step__tom__init_spheremask {
    for_each_step__tom__init_spheremask(TPRECISION centerx, TPRECISION centery, TPRECISION centerz, TPRECISION radius)
                                         : centerx_(centerx), centery_(centery), centerz_(centerz), radius_(radius*radius), z_(0), y_(0) { }
    TPRECISION centerx_, centery_, centerz_, radius_;
    TPRECISION z_, y_;
    inline void stepz(TPRECISION z) {
        z_ = z - centerz_;
        z_ = z_*z_; }
    inline void stepy(TPRECISION y) {
        y_ = y - centery_;
        y_ = y_*y_ + z_; }
    inline void operator()(T &a, TPRECISION x) {
        x -= centerx_;
        TPRECISION dist = x*x + y_;
        a = (dist <= radius_);
        //a = (tom::math::sqrt<TPRECISION>(x*x + this->current_y) <= this->radius);
    }
};
}
/****************************************************************************//**
 * \brief initialises volume with sphere mask.
 * \param[out] mask The mask which will be initialized.
 * \param[in] radius The radius of the sphere. Inside this radius the mask has
 *   value 0.
 * \param[in] sigma If > 0 there is a gaussian descent of the values with deviation
 *   sigma.
 * \param[in] max_radius Ignored if sigma <= 0. Otherwise it specifies, that outside
 *   of max_radius all elements are set to 0. Setting max_radius smaller than
 *   radius defaults to max_radius = radius+sqrt(2*sigma).
 * \param[in] centerx
 * \param[in] centery
 * \param[in] centerz
 *******************************************************************************/
template<typename T>
void tom::init_spheremask(tom::Volume<T> &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz) {
    typedef float TPRECISION;
    if (radius < 0.) {
        radius = 0.;
    }
    if (sigma > 0) {
        if (max_radius < radius) {
            max_radius = radius + sqrt(2.*sigma);
        }
        tom::loop::for_each_step<T, tom::Volume<T> &, ::for_each_step__tom__init_spheremask_sigma<TPRECISION, T> >(mask, ::for_each_step__tom__init_spheremask_sigma<TPRECISION, T>(centerx, centery, centerz, radius, sigma, max_radius));
    } else if (radius > 0.) {
        tom::loop::for_each_step<T, tom::Volume<T> &, ::for_each_step__tom__init_spheremask<TPRECISION, T> >(mask, ::for_each_step__tom__init_spheremask<TPRECISION, T>(centerx, centery, centerz, radius));
    } else {
        mask.setValues(0);
    }
}




namespace {
template<typename T>
struct for_each_step__tom__fourier_shell_correlation_ring {
    typedef double TPRECISION;
    typedef float  TPRECISION_RING;

    for_each_step__tom__fourier_shell_correlation_ring(std::size_t n_shells, std::size_t size):
        n_shells(n_shells),
        ampl1(NULL),
        ampl2(NULL),
        ampl_diff(NULL),
        phares1(NULL),
        phares2(NULL),
        n_shell(NULL),
        center((size+1)/2),
        sy(0),
        sz(0) {
        assert(size>0 && n_shells>0);
        try {
            ampl1       = new double[n_shells];
            ampl2       = new double[n_shells];
            ampl_diff   = new double[n_shells];
            phares1     = new double[n_shells];
            phares2     = new double[n_shells];
            n_shell     = new std::size_t[n_shells];
        } catch (...) {
            delete[] ampl1;
            delete[] ampl2;
            delete[] ampl_diff;
            delete[] phares1;
            delete[] phares2;
            delete[] n_shell;
            throw;
        }
        for (std::size_t i=0; i<n_shells; i++) {
            ampl1[i]        = 0;
            ampl2[i]        = 0;
            ampl_diff[i]    = 0;
            phares1[i]      = 0;
            phares2[i]      = 0;
            n_shell[i]      = 0;
        }
        shell_thickness = static_cast<TPRECISION_RING>(size-center) / static_cast<TPRECISION_RING>(n_shells);
        size3 = size*size*size;
    }
    ~for_each_step__tom__fourier_shell_correlation_ring() {
        delete[] ampl1;
        delete[] ampl2;
        delete[] ampl_diff;
        delete[] phares1;
        delete[] phares2;
        delete[] n_shell;
    }
    inline void stepz(signed long z) {
        sz = z - center;
        sz = sz * sz;
    }
    inline void stepy(signed long y) {
        sy = y - center;
        sy = sz + sy*sy;
    }
    inline void operator()(const std::complex<T> &f1, const std::complex<T> &f2, signed long x) {
        x -= center;
        if ((x = sy + x*x)) {
            const TPRECISION_RING shelld = tom::math::sqrt<TPRECISION_RING>(x) / shell_thickness;
            if (shelld < n_shells) {
                const std::size_t shell = static_cast<std::size_t>(shelld);
                assert(shelld>=0 && shell<n_shells);

                n_shell[shell]++;

                const TPRECISION f1_real = f1.real();
                const TPRECISION f1_imag = f1.imag();
                const TPRECISION f2_real = f2.real();
                const TPRECISION f2_imag = f2.imag();
                const TPRECISION amplitude1 = (f1_real*f1_real + f1_imag*f1_imag) / size3;
                const TPRECISION amplitude2 = (f2_real*f2_real + f2_imag*f2_imag) / size3;
                const TPRECISION f12_real = f1_real - f2_real;
                const TPRECISION f12_imag = f1_imag - f2_imag;
                const TPRECISION amplitude_diff = (f12_real*f12_real + f12_imag*f12_imag) / size3;
                const TPRECISION phares_v = tom::math::sqrt<TPRECISION>(amplitude1) + tom::math::sqrt<TPRECISION>(amplitude2); // Here there is a difference from tom_comparec.c

                ampl1[shell] += amplitude1;
                ampl2[shell] += amplitude2;
                ampl_diff[shell] += amplitude_diff;
                phares1[shell] += phares_v;

                TPRECISION arg = 2*tom::math::sqrt<TPRECISION>(amplitude1*amplitude2); /* calculation of cos(theta) for Phase Residuum */
                if (arg > 0) {
                    arg = (amplitude1 + amplitude2 - amplitude_diff)/arg;
                    if (arg>1) {
                        arg = 1;
                    } else if (arg<-1) {
                        arg = -1;
                    }
                }
                const TPRECISION delta = tom::math::acos<TPRECISION>(arg);    /* Phaseshift */
                phares2[shell] += phares_v*delta*delta;      /* Phase Residuum */
            }
        }
    }


    std::size_t n_shells;
    double *ampl1, *ampl2, *ampl_diff, *phares1, *phares2;
    std::size_t *n_shell;
    TPRECISION_RING shell_thickness;
    signed long center;
    signed long sy, sz;
    TPRECISION size3;
};
}
/****************************************************************************//**
 * The function work fastest, it the volumes are already fft-shifted and
 * reduced complex. Otherwise temporary copies have to be made.
 *******************************************************************************/
template<typename T>
void tom::fourier_shell_correlation(const tom::Volume<std::complex<T> > &F1_, bool is_fft_shifted1, const tom::Volume<std::complex<T> > &F2_, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &out) {

    if (F1_.getSizeX()!=size || F2_.getSizeX()!=size || F1_.getSizeY()!=size || F2_.getSizeY()!=size || (F1_.getSizeZ()!=size && F1_.getSizeZ()!=size/2+1) || (F2_.getSizeZ()!=size && F2_.getSizeZ()!=size/2+1)) {
        std::invalid_argument("The input volumes must have the same size (the z-dimension can opionally be reduced).");
    }
    const tom::Volume<std::complex<T> > *F1 = &F1_;
    const tom::Volume<std::complex<T> > *F2 = &F2_;

    // Prepare the input volumes so that they are reduced complex and fftshifted
    std::auto_ptr<tom::Volume<std::complex<T> > > F1p, F1p_half, F2p, F2p_half, vtmp;
    if (F1_.getSizeZ() == size/2+1) {
        if (!is_fft_shifted1) {
            vtmp.reset(    new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F1p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F1p_half.reset(new tom::Volume<std::complex<T> >(*F1p, NULL, size, size, size/2+1, F1p->getStrideX(), F1p->getStrideY(), F1p->getStrideZ()));
            tom::hermitian_symmetry_to_full(F1_, *vtmp);
            tom::fftshift(*vtmp, *F1p, false);
            F1 = F1p_half.get();
        }
    } else {
        if (is_fft_shifted1) {
            F1p_half.reset(new tom::Volume<std::complex<T> >(const_cast<std::complex<T> *>(&F1_.get()), size, size, size/2+1, F1_.getStrideX(), F1_.getStrideY(), F1_.getStrideZ(), false, NULL));
        } else {
            F1p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F1p_half.reset(new tom::Volume<std::complex<T> >(*F1p, NULL, size, size, size/2+1, F1p->getStrideX(), F1p->getStrideY(), F1p->getStrideZ()));
            tom::fftshift(F1_, *F1p, false);
        }
        F1 = F1p_half.get();
    }
    if (F2_.getSizeZ() == size/2+1) {
        if (!is_fft_shifted2) {
            if (!vtmp.get()) { vtmp.reset(    new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL)); }
            F2p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F2p_half.reset(new tom::Volume<std::complex<T> >(*F2p, NULL, size, size, size/2+1, F2p->getStrideX(), F2p->getStrideY(), F2p->getStrideZ()));
            tom::hermitian_symmetry_to_full(F2_, *vtmp);
            tom::fftshift(*vtmp, *F2p, false);
            F2 = F2p_half.get();
        }
    } else {
        if (is_fft_shifted2) {
            F2p_half.reset(new tom::Volume<std::complex<T> >(const_cast<std::complex<T> *>(&F2_.get()), size, size, size/2+1, F2_.getStrideX(), F2_.getStrideY(), F2_.getStrideZ(), false, NULL));
        } else {
            F2p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F2p_half.reset(new tom::Volume<std::complex<T> >(*F2p, NULL, size, size, size/2+1, F2p->getStrideX(), F2p->getStrideY(), F2p->getStrideZ()));
            tom::fftshift(F2_, *F2p, false);
        }
        F2 = F2p_half.get();
    }
    vtmp.reset(NULL);

    ::for_each_step__tom__fourier_shell_correlation_ring<T> s(n_shells, size);
    tom::loop::for_each_step<const std::complex<T>, const tom::Volume<std::complex<T> > &, const std::complex<T>, const tom::Volume<std::complex<T> > &, ::for_each_step__tom__fourier_shell_correlation_ring<T> &>(*F1, *F2, s);

    std::size_t ps = 0;
    out.resize(10 * n_shells);
    for (std::size_t i=0; i<n_shells; i++) {
        double ccc, mean;

        ccc = s.ampl1[i] * s.ampl2[i];
        assert(ccc >= 0);
        ccc  = (ccc>0)          ? (s.ampl1[i] + s.ampl2[i] - s.ampl_diff[i]) / sqrt(ccc) / 2.
                                : 1.;
        mean = (s.n_shell[i])   ? (s.ampl1[i] + s.ampl2[i] - s.ampl_diff[i]) / (s.n_shell[i]*2)
                                : 0.;
        if(mean < 0) {
            mean = -mean;
        }
        mean = sqrt(mean);
        const double rmsd = tom::math::sqrt<double>(s.ampl_diff[i] / (2*(s.ampl1[i] + s.ampl2[i]) - s.ampl_diff[i]));
        s.ampl1[i] = tom::math::sqrt<double>(s.ampl1[i] / s.n_shell[i]);
        s.ampl2[i] = tom::math::sqrt<double>(s.ampl2[i] / s.n_shell[i]);
        s.ampl_diff[i] = tom::math::sqrt<double>(s.ampl_diff[i]/s.n_shell[i]);
        const double sqrn = 2/tom::math::sqrt<double>(s.n_shell[i]);                 /* 2 / sqrt(N) */
        const double phares = tom::math::sqrt<double>(s.phares2[i]/s.phares1[i])*57.29577951308232087665461840231273527024; /* 180/pi */      /* Phase Residuum in degrees */

        out[i+n_shells*0] = s.n_shell[i];
        out[i+n_shells*1] = 2.*static_cast<double>(n_shells)/static_cast<double>(i+1);               /* pixel resolution */
        out[i+n_shells*2] = s.ampl_diff[i];
        out[i+n_shells*3] = rmsd;
        out[i+n_shells*4] = s.ampl1[i];
        out[i+n_shells*5] = s.ampl2[i];
        out[i+n_shells*6] = mean;
        out[i+n_shells*7] = ccc;
        out[i+n_shells*8] = sqrn;
        out[i+n_shells*9] = phares;

        ps += s.n_shell[i];
    }
}








// template instantiation
template int tom::get_tom_io_type<int                  >();
template int tom::get_tom_io_type<double               >();
template int tom::get_tom_io_type<fftw_complex         >();
template int tom::get_tom_io_type<std::complex<double> >();


template void tom::norm_mask<float , float , double>(tom::Volume<float > &v, const tom::Volume<float > &mask, tom::norm::ntype stddev_type, double *variance, bool is_boolean_mask);
template void tom::norm_mask<double, double, double>(tom::Volume<double> &v, const tom::Volume<double> &mask, tom::norm::ntype stddev_type, double *variance, bool is_boolean_mask);

template void tom::fftshift<char                 >(const Volume<char                 > &vsrc, Volume<char                 > &v, bool is_ifftshift);
template void tom::fftshift<float                >(const Volume<float                > &vsrc, Volume<float                > &v, bool is_ifftshift);
template void tom::fftshift<double               >(const Volume<double               > &vsrc, Volume<double               > &v, bool is_ifftshift);
template void tom::fftshift<std::complex<float > >(const Volume<std::complex<float > > &vsrc, Volume<std::complex<float > > &v, bool is_ifftshift);
template void tom::fftshift<std::complex<double> >(const Volume<std::complex<double> > &vsrc, Volume<std::complex<double> > &v, bool is_ifftshift);

template std::vector<tom::st_idx> tom::peak<float >(const tom::Volume<float > &v);
template std::vector<tom::st_idx> tom::peak<double>(const tom::Volume<double> &v);
template std::vector<tom::st_idx> tom::peak<float , char>(const tom::Volume<float > &v, const tom::Volume<char> &mask);
template std::vector<tom::st_idx> tom::peak<double, char>(const tom::Volume<double> &v, const tom::Volume<char> &mask);
template std::vector<tom::st_idx> tom::peak<float , float >(const tom::Volume<float > &v, const tom::Volume<float > &mask);
template std::vector<tom::st_idx> tom::peak<double, double>(const tom::Volume<double> &v, const tom::Volume<double> &mask);


template void tom::hermitian_symmetry_to_full<float                >(const tom::Volume<float                > &vsrc, tom::Volume<float                > &v);
template void tom::hermitian_symmetry_to_full<double               >(const tom::Volume<double                > &vsrc, tom::Volume<double              > &v);
template void tom::hermitian_symmetry_to_full<std::complex<float > >(const tom::Volume<std::complex<float > > &vsrc, tom::Volume<std::complex<float > > &v);
template void tom::hermitian_symmetry_to_full<std::complex<double> >(const tom::Volume<std::complex<double> > &vsrc, tom::Volume<std::complex<double> > &v);

template void tom::init_spheremask<char  >(tom::Volume<char  > &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);
//template void tom::init_spheremask<int   >(tom::Volume<int   > &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);
template void tom::init_spheremask<float >(tom::Volume<float > &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);
template void tom::init_spheremask<double>(tom::Volume<double> &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);

template void tom::update_correlation_volume<float , int32_t>(const Volume<float > &v, Volume<float > &cc, Volume<int32_t> &vindex, int32_t index);
template void tom::update_correlation_volume<double, int32_t>(const Volume<double> &v, Volume<double> &cc, Volume<int32_t> &vindex, int32_t index);

template void tom::rotate(const tom::Volume<float > &src, tom::Volume<float > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post);
template void tom::rotate(const tom::Volume<double> &src, tom::Volume<double> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post);


template void tom::element_wise_add<float , float >(tom::Volume<float > &b, const tom::Volume<float > &a);
template void tom::element_wise_add<double, double>(tom::Volume<double> &b, const tom::Volume<double> &a);

template void tom::element_wise_sub<float , float >(tom::Volume<float > &b, const tom::Volume<float > &a);
template void tom::element_wise_sub<double, double>(tom::Volume<double> &b, const tom::Volume<double> &a);

template void tom::element_wise_multiply<std::complex<float >, float >(tom::Volume<std::complex<float > > &b, const tom::Volume<float > &a);
template void tom::element_wise_multiply<std::complex<float >, double>(tom::Volume<std::complex<float > > &b, const tom::Volume<double> &a);
template void tom::element_wise_multiply<std::complex<double>, float >(tom::Volume<std::complex<double> > &b, const tom::Volume<float > &a);
template void tom::element_wise_multiply<std::complex<double>, double>(tom::Volume<std::complex<double> > &b, const tom::Volume<double> &a);

template void tom::element_wise_multiply<float,float>(tom::Volume<float > &b, const tom::Volume<float > &a);
template void tom::element_wise_multiply<double,double>(tom::Volume<double > &b, const tom::Volume<double > &a);

template void tom::element_wise_multiply<std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float > > &a);
template void tom::element_wise_multiply<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a);


template void tom::element_wise_div<std::complex<float >, float >(tom::Volume<std::complex<float > > &b, const tom::Volume<float > &a, std::complex<float > inf_value);
template void tom::element_wise_div<std::complex<double>, double>(tom::Volume<std::complex<double> > &b, const tom::Volume<double> &a, std::complex<double> inf_value);
template void tom::element_wise_div<float , float>(tom::Volume<float>  &b, const tom::Volume<float> &a, float inf_value);
template void tom::element_wise_div<double , double>(tom::Volume<double> &b, const tom::Volume<double> &a, double inf_value);

template void tom::element_wise_conj_multiply<std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float > > &a);
template void tom::element_wise_conj_multiply<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a);

template void tom::element_wise_conj_multiply<std::complex<float >, std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &v1, const tom::Volume<std::complex<float > > &v2, const tom::Volume<std::complex<float > > &v3);
template void tom::element_wise_conj_multiply<std::complex<double>, std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &v1, const tom::Volume<std::complex<double> > &v2, const tom::Volume<std::complex<double> > &v3);


template void tom::element_wise_power<float >(tom::Volume<float > &v, float  exponent);
template void tom::element_wise_power<double>(tom::Volume<double> &v, double exponent);

template void tom::element_wise_max<float >(tom::Volume<float > &v, float  val);
template void tom::element_wise_max<double>(tom::Volume<double> &v, double val);



template void tom::make_binary<float >(tom::Volume<float > &v);
template void tom::make_binary<double>(tom::Volume<double> &v);


template void tom::fourier_shell_correlation<float >(const tom::Volume<std::complex<float > > &F1, bool is_fft_shifted1, const tom::Volume<std::complex<float > > &F2, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &r);
template void tom::fourier_shell_correlation<double>(const tom::Volume<std::complex<double> > &F1, bool is_fft_shifted1, const tom::Volume<std::complex<double> > &F2, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &r);

template void tom::element_wise_set_below_threshold<float >(tom::Volume<float > &a, float  threshold, float  value);
template void tom::element_wise_set_below_threshold<double>(tom::Volume<double> &a, double threshold, double value);



