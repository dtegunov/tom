/***********************************************************************//**
 * \file transform.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    18.12.2007
 **************************************************************************/
#include "tom/core/transform.hpp"




#include "tom/core/volume.hpp"
#include "tom/core/transform.h"





namespace {
template<typename T> inline bool is_simple_interpolationfcn                                 () { return false; }
template<>           inline bool is_simple_interpolationfcn<tom::InterpolTriLinear<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::InterpolTriLinear<double> >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::InterpolNearestNeighbour<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::InterpolNearestNeighbour<double> >() { return true;  }
}







/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
void tom::InterpolTriLinear<T>::setvolume(const tom::Volume<T> &src) {

    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();

    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw std::invalid_argument("The volume of InterpolTriLinear must be aligned with sizeof(T).");
    }

    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->sizex_m1 = this->sizex - 1;
    this->sizey_m1 = this->sizey - 1;
    this->sizez_m1 = this->sizez - 1;
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();

}





/***********************************************************************//**
 *
 **************************************************************************/
template<typename T>
void tom::InterpolNearestNeighbour<T>::setvolume(const tom::Volume<T> &src) {

    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();

    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw std::invalid_argument("The volume of InterpolNearestNeighbour must be aligned with sizeof(T).");
    }

    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();

}






/***********************************************************************//**
 *
 **************************************************************************/
void tom::InterpolInterpolate3D::setvolume(const tom::Volume<tom::InterpolInterpolate3D::val_type> &src) {

    if (!src.isContiguous()) {
        throw std::invalid_argument("Interpolate3D is only implemented for contiguous data.");
    }
    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->sizex_m1 = this->sizex - 1;
    this->sizey_m1 = this->sizey - 1;
    this->sizez_m1 = this->sizez - 1;

    switch (this->method) {
        case __TOM_INTERPOL_CUBIC:
            this->borderx1 = 1;
            this->bordery1 = 1;
            this->borderz1 = 2;
            this->borderx2 = this->sizex - 1;
            this->bordery2 = this->sizey - 1;
            this->borderz2 = this->sizez - 2;
            break;
        case __TOM_INTERPOL_CSPLINE:
            this->borderx1 = 2;
            this->bordery1 = 2;
            this->borderz1 = 2;
            this->borderx2 = this->sizex - 2;
            this->bordery2 = this->sizey - 2;
            this->borderz2 = this->sizez - 2;
            break;
        default:
            this->borderx1 = 0;
            this->bordery1 = 0;
            this->borderz1 = 0;
            this->borderx2 = 0;
            this->bordery2 = 0;
            this->borderz2 = 0;
            break;
    }


    this->data = &src.get();
}


/***********************************************************************//**
 *
 **************************************************************************/
template<typename T, typename TINTERP>
void tom::transform(const tom::Volume<T> &src, tom::Volume<T> &dst, const double *P, int is_affinity, T valinf, TINTERP interp) {


    std::size_t src_stridex = src.getStrideX();
    std::size_t src_stridey = src.getStrideY();
    std::size_t src_stridez = src.getStrideZ();
    std::size_t dst_stridex = dst.getStrideX();
    std::size_t dst_stridey = dst.getStrideY();
    std::size_t dst_stridez = dst.getStrideZ();

    if (&src.get() == &dst.get()) {
        throw std::invalid_argument("Self assignment not allowed.");
    }

    {
        // If the volumes are not alligned to the multiples of sizeof(T) create a local copy and work with that.
        const bool src_aligned = !(src_stridex!=sizeof(T) || src_stridey%sizeof(T) || src_stridez%sizeof(T));
        const bool dst_aligned = !(dst_stridex!=sizeof(T) || dst_stridey%sizeof(T) || dst_stridez%sizeof(T));
        if (!src_aligned || !dst_aligned) {
            std::auto_ptr<tom::Volume<T> > vsrc_aligned_, vdst_aligned_;
            const tom::Volume<T> *psrc = &src;
            tom::Volume<T> *pdst = &dst;
            if (!src_aligned) {
                vsrc_aligned_.reset(new tom::Volume<T>(src));
                psrc = vsrc_aligned_.get();
            }
            if (!dst_aligned) {
                vdst_aligned_.reset(new tom::Volume<T>(dst));
                pdst = vdst_aligned_.get();
            }
            tom::transform<T, TINTERP>(*psrc, *pdst, P, is_affinity, valinf, interp);
        }
    }

    src_stridex /= sizeof(T); src_stridey /= sizeof(T); src_stridez /= sizeof(T);
    dst_stridex /= sizeof(T); dst_stridey /= sizeof(T); dst_stridez /= sizeof(T);

    double Pinv[16];

    {
        // Compute the inverse of the transformation matrix.
        {
            Pinv[ 0] = P[ 0];   Pinv[ 1] = P[ 1];   Pinv[ 2] = P[ 2];   Pinv[ 3] = P[ 3];
            Pinv[ 4] = P[ 4];   Pinv[ 5] = P[ 5];   Pinv[ 6] = P[ 6];   Pinv[ 7] = P[ 7];
            Pinv[ 8] = P[ 8];   Pinv[ 9] = P[ 9];   Pinv[10] = P[10];   Pinv[11] = P[11];
        }
        if (is_affinity) {
            Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
        } else {
            Pinv[12] = P[12];   Pinv[13] = P[13];   Pinv[14] = P[14];   Pinv[15] = P[15];
            is_affinity = Pinv[12]==0. && Pinv[13]==0. && Pinv[14]==0. && Pinv[15]!=0.;
            if (is_affinity && Pinv[15]!=1.) {
                for (int i=0; i<12; i++) {
                    Pinv[i] /= Pinv[15];
                }
            }
        }
        if (!tom_transf_invert(Pinv, Pinv, 4)) {
            throw std::invalid_argument("The transformation matrix is non invertible");
        }
        if (is_affinity) {
            Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
        }
    }

    const std::size_t dst_sizex = dst.getSizeX();
    const std::size_t dst_sizey = dst.getSizeY();
    const std::size_t dst_sizez = dst.getSizeZ();

    if (is_affinity  &&
        Pinv[ 0]==1. && Pinv[ 1]==0. && Pinv[ 2]==0. && Pinv[ 3]==0. &&
        Pinv[ 4]==0. && Pinv[ 5]==1. && Pinv[ 6]==0. && Pinv[ 7]==0. &&
        Pinv[ 8]==0. && Pinv[ 9]==0. && Pinv[10]==1. && Pinv[11]==0. &&
        is_simple_interpolationfcn<TINTERP>()) {
        /* Projection is identity. Do only copy
           Only in case of certain interpolation classes (is_simple_interpolationfcn).
           Otherwise it may be desired to do the "interpolation" nonetheless,
           for example, if the interpolation-function is a filter in time domain over
           several neighbours in the volume space. */
        dst.setValues(tom::Volume<T>(const_cast<T *>(&src.get()), dst_sizex,dst_sizey,dst_sizez, src.getStrideX(), src.getStrideY(), src.getStrideZ(), false, NULL));
    } else {

        typedef typename TINTERP::idx_type FLOATTYPE;


        T *pdst = &dst.get();

        /* Copy the projective transformation into single variables. */
        const FLOATTYPE p00 = Pinv[ 0];
        const FLOATTYPE p01 = Pinv[ 1];
        const FLOATTYPE p02 = Pinv[ 2];
        const FLOATTYPE p03 = Pinv[ 3];
        const FLOATTYPE p10 = Pinv[ 4];
        const FLOATTYPE p11 = Pinv[ 5];
        const FLOATTYPE p12 = Pinv[ 6];
        const FLOATTYPE p13 = Pinv[ 7];
        const FLOATTYPE p20 = Pinv[ 8];
        const FLOATTYPE p21 = Pinv[ 9];
        const FLOATTYPE p22 = Pinv[10];
        const FLOATTYPE p23 = Pinv[11];
        const FLOATTYPE p30 = Pinv[12];
        const FLOATTYPE p31 = Pinv[13];
        const FLOATTYPE p32 = Pinv[14];
        const FLOATTYPE p33 = Pinv[15];
        FLOATTYPE tmpz_00, tmpz_01, tmpz_02;
        FLOATTYPE tmpy_00, tmpy_01, tmpy_02;
        std::size_t dstx, dsty, dstz;

        interp.setvolume(src);

        dst_stridez -= dst_sizey*dst_stridey;
        dst_stridey -= dst_sizex;
        if (is_affinity) {
            for (dstz=0; dstz<dst_sizez; dstz++) {
                tmpz_00 = p02*dstz + p03;
                tmpz_01 = p12*dstz + p13;
                tmpz_02 = p22*dstz + p23;
                for (dsty=0; dsty<dst_sizey; dsty++) {
                    tmpy_00 = p01*dsty + tmpz_00;
                    tmpy_01 = p11*dsty + tmpz_01;
                    tmpy_02 = p21*dsty + tmpz_02;
                    for (dstx=0; dstx<dst_sizex; dstx++) {
                        *pdst++ = interp.interpolate(p00*dstx + tmpy_00, p10*dstx + tmpy_01, p20*dstx + tmpy_02);
                    }
                    pdst += dst_stridey;
                }
                pdst += dst_stridez;
            }
        } else {
            FLOATTYPE pw, tmpz_03, tmpy_03;
            for (dstz=0; dstz<dst_sizez; dstz++) {
                tmpz_00 = p02*dstz + p03;
                tmpz_01 = p12*dstz + p13;
                tmpz_02 = p22*dstz + p23;
                tmpz_03 = p32*dstz + p33;
                for (dsty=0; dsty<dst_sizey; dsty++) {
                    tmpy_00 = p01*dsty + tmpz_00;
                    tmpy_01 = p11*dsty + tmpz_01;
                    tmpy_02 = p21*dsty + tmpz_02;
                    tmpy_03 = p31*dsty + tmpz_03;
                    for (dstx=0; dstx<dst_sizex; dstx++) {
                        if((pw = p30*dstx + tmpy_03) == 0.) {
                            *pdst++ = valinf;
                        } else {
                            *pdst++ = interp.interpolate((p00*dstx + tmpy_00) / pw, (p10*dstx + tmpy_01) / pw, (p20*dstx + tmpy_02) / pw);
                        }
                    }
                    pdst += dst_stridey;
                }
                pdst += dst_stridez;
            }
        }
    }

}




template void tom::transform<float , tom::InterpolNearestNeighbour<float > >(const tom::Volume<float > &src, tom::Volume<float > &dst, const double *P, int is_affinity, float  valinf, tom::InterpolNearestNeighbour<float > interp);
template void tom::transform<double, tom::InterpolNearestNeighbour<double> >(const tom::Volume<double> &src, tom::Volume<double> &dst, const double *P, int is_affinity, double valinf, tom::InterpolNearestNeighbour<double> interp);
template void tom::transform<float , tom::InterpolTriLinear<float > >(const tom::Volume<float > &src, tom::Volume<float > &dst, const double *P, int is_affinity, float  valinf, tom::InterpolTriLinear<float > interp);
template void tom::transform<double, tom::InterpolTriLinear<double> >(const tom::Volume<double> &src, tom::Volume<double> &dst, const double *P, int is_affinity, double valinf, tom::InterpolTriLinear<double> interp);
template void tom::transform<float , tom::InterpolInterpolate3D>(const tom::Volume<float > &src, tom::Volume<float > &dst, const double *P, int is_affinity, float  valinf, tom::InterpolInterpolate3D interp);







