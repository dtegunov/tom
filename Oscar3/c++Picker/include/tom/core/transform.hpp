/***********************************************************************//**
 * \file transform.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.12.2007
 **************************************************************************/
#ifndef ___INCLUDE_CORE__TRANSFORM_HPP__
#define ___INCLUDE_CORE__TRANSFORM_HPP__


#include "tom/core/volume.hpp"

#include "tom/core/interpol.h"

namespace tom {






/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
class InterpolTriLinear {

public:
    typedef double idx_type;
    InterpolTriLinear(T defaultval);
    void setvolume(const tom::Volume<T> &src);
    T interpolate(const typename tom::InterpolTriLinear<T>::idx_type &x, const typename tom::InterpolTriLinear<T>::idx_type &y, const typename tom::InterpolTriLinear<T>::idx_type &z);

private:
    const T *data;
    std::size_t sizex, sizey, sizez;
    typename tom::InterpolTriLinear<T>::idx_type sizex_m1, sizey_m1, sizez_m1;
    std::size_t stridey, stridez;

    T defaultval;


};




/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
class InterpolNearestNeighbour {

public:
    typedef double idx_type;
    InterpolNearestNeighbour(T defaultval);
    void setvolume(const tom::Volume<T> &src);
    T interpolate(typename tom::InterpolNearestNeighbour<T>::idx_type x, typename tom::InterpolNearestNeighbour<T>::idx_type y, typename tom::InterpolNearestNeighbour<T>::idx_type z);

private:
    const T *data;
    std::size_t sizex, sizey, sizez;
    std::size_t stridey, stridez;

    T defaultval;


};


/***********************************************************************//**
 *
 **************************************************************************/
class InterpolInterpolate3D {

public:
    typedef float val_type;
    typedef float idx_type;
    InterpolInterpolate3D(tom::InterpolInterpolate3D::val_type defaultval, short method);
    void setvolume(const tom::Volume<tom::InterpolInterpolate3D::val_type> &src);
    tom::InterpolInterpolate3D::val_type interpolate(tom::InterpolInterpolate3D::idx_type x, tom::InterpolInterpolate3D::idx_type y, tom::InterpolInterpolate3D::idx_type z);

private:
    const float *data;
    short method;
    std::size_t sizex, sizey, sizez;
    tom::InterpolInterpolate3D::idx_type sizex_m1, sizey_m1, sizez_m1;
    tom::InterpolInterpolate3D::idx_type borderx1, bordery1, borderz1, borderx2, bordery2, borderz2;


    tom::InterpolInterpolate3D::val_type defaultval;
};




template<typename T, typename TINTERP>
void transform(const tom::Volume<T> &src, tom::Volume<T> &dst, const double *P, int is_affinity, T valinf, TINTERP interp);





}







// Inline functions






/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
inline tom::InterpolTriLinear<T>::InterpolTriLinear(T defaultval):
    data(NULL),
    sizex(0), sizey(0), sizez(0), sizex_m1(0), sizey_m1(0), sizez_m1(0), stridey(0), stridez(0),
    defaultval(defaultval) {
}


/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
inline tom::InterpolNearestNeighbour<T>::InterpolNearestNeighbour(T defaultval):
    data(NULL),
    sizex(0), sizey(0), sizez(0), stridey(0), stridez(0),
    defaultval(defaultval) {
}

/***********************************************************************//**
 *
 **************************************************************************/
inline tom::InterpolInterpolate3D::InterpolInterpolate3D(tom::InterpolInterpolate3D::val_type defaultval, short method):
    data(NULL),
    method(method),
    sizex(0), sizey(0), sizez(0),
    sizex_m1(0), sizey_m1(0), sizez_m1(0),
    borderx1(0), bordery1(0), borderz1(0),
    borderx2(0), bordery2(0), borderz2(0),
    defaultval(defaultval) {

    switch (method) {
        case __TOM_INTERPOL_NEAREST:
        case __TOM_INTERPOL_LINEAR:
        case __TOM_INTERPOL_CUBIC:
        case __TOM_INTERPOL_CSPLINE:
            break;
        default:
            throw std::invalid_argument("Non implemented method for Interpolate3D");
    }
}



/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
inline T tom::InterpolTriLinear<T>::interpolate(const typename tom::InterpolTriLinear<T>::idx_type &x, const typename tom::InterpolTriLinear<T>::idx_type &y, const typename tom::InterpolTriLinear<T>::idx_type &z) {

    typedef typename tom::InterpolTriLinear<T>::idx_type FLOATTYPE;

    if (x<0. || y<0. || z<0. ||
        x > (this->sizex_m1) ||
        y > (this->sizey_m1) ||
        z > (this->sizez_m1)) {
        return this->defaultval;
    }
    const size_t floorx = static_cast<std::size_t>(x);
    const size_t floory = static_cast<std::size_t>(y);
    const size_t floorz = static_cast<std::size_t>(z);
    const FLOATTYPE xoffseth = x - floorx;
    const FLOATTYPE yoffseth = y - floory;
    const FLOATTYPE zoffseth = z - floorz;

    const T *src = &this->data[floorz * this->stridez + floory*this->stridey + floorx];

    switch (static_cast<int>(xoffseth > 0.) | (static_cast<int>(yoffseth > 0.)<<1) | (static_cast<int>(zoffseth > 0.)<<2)) {
        case 0x07: {
                FLOATTYPE x00x, x01x, x10x, x11x, x0yx;
                x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridey;
                x01x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridez;
                x11x = src[0] + (src[1]-src[0])*xoffseth;
                src -= this->stridey;
                x10x = src[0] + (src[1]-src[0])*xoffseth;
                x0yx = x00x + (x01x - x00x) * yoffseth;
                return x0yx + ((x10x + (x11x - x10x) * yoffseth) - x0yx) * zoffseth;
            }
        case 0x06: {
                const FLOATTYPE x0y0 = src[0] + (src[this->stridey]-src[0])*yoffseth;
                src += this->stridez;
                return x0y0 + (src[0] + (src[this->stridey]-src[0])*yoffseth - x0y0)*zoffseth;
            }
        case 0x05: {
                const FLOATTYPE x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridez;
                return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)*zoffseth;
            }
        case 0x03: {
                const FLOATTYPE x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridey;
                return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)* yoffseth;
            }
        case 0x04:
            return src[0] + (src[this->stridez]-src[0])*zoffseth;
        case 0x02:
            return src[0] + (src[this->stridey]-src[0])*yoffseth;
        case 0x01:
            return src[0] + (src[1]           -src[0])*xoffseth;
    }
    return src[0];
}




/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
inline T tom::InterpolNearestNeighbour<T>::interpolate(typename tom::InterpolNearestNeighbour<T>::idx_type x, typename tom::InterpolNearestNeighbour<T>::idx_type y, typename tom::InterpolNearestNeighbour<T>::idx_type z) {

    x += 0.5;
    y += 0.5;
    z += 0.5;
    if (x<0 || y<0 || z<0 || x>=this->sizex || y>=this->sizey || z>=this->sizez) {
        return this->defaultval;
    }
    return this->data[static_cast<std::size_t>(z)*this->stridez +
                       static_cast<std::size_t>(y)*this->stridey +
                       static_cast<std::size_t>(x)];

}



/***********************************************************************//**
 *
 **************************************************************************/
inline tom::InterpolInterpolate3D::val_type tom::InterpolInterpolate3D::interpolate(
                            tom::InterpolInterpolate3D::idx_type x,
                            tom::InterpolInterpolate3D::idx_type y,
                            tom::InterpolInterpolate3D::idx_type z) {

    if (x<0 || x>this->sizex_m1 || y<0 || y>this->sizey_m1 || z<0 || z>this->sizez_m1) {
        return this->defaultval;
    }

    switch (this->method) {
        case __TOM_INTERPOL_NEAREST:
        case __TOM_INTERPOL_LINEAR:
            return ::interpolate3D(this->data, this->sizex, this->sizey, this->sizez, x, y, z, this->method);
        //case __TOM_INTERPOL_CUBIC:
        //case __TOM_INTERPOL_CSPLINE:
        default:
            if (x<this->borderx1 || x>this->borderx2 || y<this->bordery1 || y>this->bordery2 || z<this->borderz1 || z>this->borderz2) {
                return interpolate3D(this->data, this->sizex, this->sizey, this->sizez, x, y, z, __TOM_INTERPOL_LINEAR);
            }
            return ::interpolate3D(this->data, this->sizex, this->sizey, this->sizez, x, y, z, this->method);
    }
}


#endif



