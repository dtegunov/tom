/***********************************************************************//**
 * \file wedge.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    02.12.2007
 **************************************************************************/
#include "tom/core/wedge.hpp"


#include <iostream>
#include <assert.h>
#include <typeinfo>

#include <boost/lexical_cast.hpp>


#include "tom/core/transform.h"
#include "tom/core/volume_fcn.hpp"



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::SimpleWedge<T>::SimpleWedge(const tom::SimpleWedge<T> &v) {
    assert(0);
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::SimpleWedge<T>::init_wedge_volume(std::size_t sizex, std::size_t sizey, std::size_t sizez, bool reduced) const {

    std::auto_ptr<tom::Volume<T> > old_wedge_reduced, old_wedge_full;

    if (this->phi!=this->current_phi || this->psi!=this->current_psi || this->theta!=this->current_theta) {
        old_wedge_reduced.reset(this->wedge_reduced);
        this->wedge_reduced = NULL;
        old_wedge_full.reset(this->wedge_full);
        this->wedge_full = NULL;
    }
    this->current_phi = this->phi;
    this->current_psi = this->psi;
    this->current_theta = this->theta;

    // Check if the value is already initialized. if so, return.
    if (reduced) {
        if (this->wedge_reduced) {
            if (this->wedge_reduced->getSizeX()==sizex && this->wedge_reduced->getSizeY()==sizey && this->wedge_reduced->getSizeZ()==sizez/2+1) {
                return;
            }
            old_wedge_reduced.reset(this->wedge_reduced);
            this->wedge_reduced = NULL;
        }
    } else {
        if (this->wedge_full) {
            if (this->wedge_full->getSizeX()==sizex && this->wedge_full->getSizeY()==sizey && this->wedge_full->getSizeZ()==sizez) {
                return;
            }
            old_wedge_full.reset(this->wedge_full);
            this->wedge_full = NULL;
        }
    }

    if (reduced) {
        if (this->wedge_full) {
            if (this->wedge_full->getSizeX()==sizex && this->wedge_full->getSizeY()==sizey && this->wedge_full->getSizeZ()==sizez) {
                this->wedge_reduced = new tom::Volume<T>(*this->wedge_full, NULL, sizex, sizey, sizez/2+1, this->wedge_full->getStrideX(), this->wedge_full->getStrideY(), this->wedge_full->getStrideZ());
                return;
            }
            delete this->wedge_full;
            this->wedge_full = NULL;
        }
        if (old_wedge_reduced.get() &&  old_wedge_reduced.get()->getSizeX()==sizex && old_wedge_reduced.get()->getSizeY()==sizey && old_wedge_reduced.get()->getSizeZ()==sizez/2+1) {
            this->wedge_reduced = old_wedge_reduced.release();
        } else {
            this->wedge_reduced = new tom::Volume<T>(sizex, sizey, sizez/2+1, NULL, NULL);
        }
        this->init_wedge_volume(*this->wedge_reduced, sizez);
    } else {
        if (old_wedge_full.get() &&  old_wedge_full.get()->getSizeX()==sizex && old_wedge_full.get()->getSizeY()==sizey && old_wedge_full.get()->getSizeZ()==sizez) {
            this->wedge_full = old_wedge_full.release();
        } else {
            this->wedge_full = new tom::Volume<T>(sizex, sizey, sizez, NULL, NULL);
        }
        if (this->wedge_reduced) {
            if (this->wedge_reduced->getSizeX()==sizex && this->wedge_reduced->getSizeY()==sizey && this->wedge_reduced->getSizeZ()==sizez/2+1) {
                tom::Volume<T>(*this->wedge_full, NULL, sizex, sizey, sizez/2+1, this->wedge_full->getStrideX(), this->wedge_full->getStrideY(), this->wedge_full->getStrideZ()
                    ).setValues(*this->wedge_reduced);
                ::tom::hermitian_symmetry_to_full<T>(*this->wedge_full);
                delete this->wedge_reduced;
                this->wedge_reduced = NULL;
                return;
            }
            delete this->wedge_reduced;
            this->wedge_reduced = NULL;
        }
        this->init_wedge_volume(*this->wedge_full, sizez);
    }
}







/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::SimpleWedge<T>::rotate(double phi, double psi, double theta) {
    this->phi = phi;
    this->psi = psi;
    this->theta = theta;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline tom::SimpleWedge<T>::SimpleWedge(double angle, double cutoff_radius)
    :   wedge_reduced(NULL),
        wedge_full(NULL),
        current_phi(0.),
        current_psi(0.),
        current_theta(0.),
        angle(angle),
        phi(0.),
        psi(0.),
        theta(0.),
        cutoff_radius(cutoff_radius) {
}





/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::SimpleWedge<T>::init_wedge_volume(tom::Volume<T> &v, std::size_t sizez) const {


    assert(sizez==v.getSizeZ() || sizez/2+1==v.getSizeZ());

    T *vdata = &v.get();
    const std::size_t sizex = v.getSizeX();
    const std::size_t sizey = v.getSizeY();
    const std::size_t v_sizez = v.getSizeZ();
    const std::size_t offsetx = sizex/2 + sizex%2;
    const std::size_t offsety = sizey/2 + sizey%2;
    const std::size_t offsetz = sizez/2 + sizez%2; // For the fftshift.

    assert(v.getStrideX()==sizeof(T) && v.getStrideY()==sizex*v.getStrideX() && v.getStrideZ()==sizey*v.getStrideY());

    typedef float TCOMP;

    TCOMP centerx = ceil((sizex-1) / 2.);
    TCOMP centery = ceil((sizey-1) / 2.);
    TCOMP centerz = ceil((sizez-1) / 2.);

    TCOMP P[16] = {   0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 1. };
    {
        double Plocal[16] = { 0., 0., 0., 0.,
                              0., 0., 0., 0.,
                              0., 0., 0., 0.,
                              0., 0., 0., 1. };
        const int axes[3]      = {   2,     0,   2 };
        const double angles[3] = { this->phi, this->theta, this->psi };
        tom_transf_init_rotmat_3d(Plocal, 4, 3, axes, angles);
        Plocal[ 3] = - (Plocal[ 0]*centerx + Plocal[ 1]*centery + Plocal[ 2]*centerz);
        Plocal[ 7] = - (Plocal[ 4]*centerx + Plocal[ 5]*centery + Plocal[ 6]*centerz);
        Plocal[11] = - (Plocal[ 8]*centerx + Plocal[ 9]*centery + Plocal[10]*centerz);
        for (int i=0; i<16; i++) {
            P[i] = Plocal[i];
        }
    }
    const TCOMP maxdist = sqrt(centerx*centerx + centery*centery + centerz*centerz);

    const TCOMP tan_angle = tan(this->angle);
    const TCOMP cutoff_radius = this->cutoff_radius>0&&this->cutoff_radius<maxdist ? this->cutoff_radius : 0.;

    const TCOMP z0_threshold = 1e-4;

    std::size_t x, y, z;
    std::size_t xs, ys, zs;
    TCOMP x_, y_, z_;
    TCOMP x_tmpz, y_tmpz, z_tmpz;
    TCOMP x_tmpy, y_tmpy, z_tmpy;
    if (cutoff_radius) {
        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%sizez;
            x_tmpz = P[ 2]*zs + P[ 3];
            y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<sizey; y++) {
                ys = (y + offsety)%sizey;
                x_tmpy = P[ 1]*ys + x_tmpz;
                y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<sizex; x++) {
                    xs = (x + offsetx)%sizex;
                    x_ =                P[ 0]*xs + x_tmpy;
                    y_ =                P[ 4]*xs + y_tmpy;
                    z_ = tom::math::abs(P[ 8]*xs + z_tmpy);
                    *vdata++ = (sqrt(x_*x_ + y_*y_ + z_*z_) <= cutoff_radius) &&
                               ((z_<z0_threshold) || (tan_angle <= (tom::math::abs<TCOMP>(x_) / z_)));
                }
            }
        }
    } else {
        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%sizez;
            x_tmpz = P[ 2]*zs + P[ 3];
            //y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<sizey; y++) {
                ys = (y + offsety)%sizey;
                x_tmpy = P[ 1]*ys + x_tmpz;
                //y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<sizex; x++) {
                    xs = (x + offsetx)%sizex;
                    x_ = P[ 0]*xs + x_tmpy;
                    //y_ = P[ 4]*xs + y_tmpy;
                    z_ = P[ 8]*xs + z_tmpy;
                    *vdata++ = (z_<z0_threshold) || (tan_angle <= (tom::math::abs<TCOMP>(x_) / tom::math::abs<TCOMP>(z_)));
                }
            }
        }
    }
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
bool tom::SimpleWedge<T>::apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez) const {

    bool reduced = false;
    if (sizez != vsrc.getSizeZ()) {
        if (sizez/2+1 != vsrc.getSizeZ()) {
            throw std::invalid_argument("either sizeZ or sizeZ/2+1 must be equal to the size of the volume");
        }
        reduced = true;
    }

    if (this->angle>0 || this->cutoff_radius<=tom::math::sqrt<double>((vsrc.getSizeX()/2)*(vsrc.getSizeX()/2) + (vsrc.getSizeY()/2)*(vsrc.getSizeY()/2) + (vsrc.getSizeZ()/2)*(vsrc.getSizeZ()/2)) + 1) {

        this->init_wedge_volume(vsrc.getSizeX(), vsrc.getSizeY(), sizez, reduced);

        if (reduced) {
            tom::element_wise_multiply<std::complex<T>, T>(vsrc, *this->wedge_reduced);
            //std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //this->wedge_reduced->write_to_em("/fs/home/haller/haller/DA/data/outputdir/wedge.em", NULL);
        } else {
            tom::element_wise_multiply<std::complex<T>, T>(vsrc, *this->wedge_full);
        }
        return true;

    }
    return false;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::SimpleWedge<T>::toString() const {
    std::stringstream ss;
    ss << "SimpleWedge(" << (this->getAngle()*(180./3.141592653589793238512808959406186204433)) << "°, " << this->getCutoffRadius() << ")";
    return ss.str();
}

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::VolumeWedge<T>::toString() const {
    std::stringstream ss;
    ss << "VolumeWedge(" << this->dims[0] << "x" << this->dims[1] << "x" << this->dims[2] << (this->active ? "" : "not active") << ")";
    return ss.str();
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::VolumeEMWedge<T>::toString() const {
    std::stringstream ss;
    ss << "VolumeEMWedge(" << this->filename << ":" << this->binning << ", " << this->dims[0] << "x" << this->dims[1] << "x" << this->dims[2] << (this->active ? "" : "not active") << ")";
    return ss.str();
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::SimpleWedge<T>::get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const {
    bool reduced = false;
    if (sizez != size_z_vol) {
        if (sizez/2+1 != size_z_vol) {
            throw std::invalid_argument("either sizez or sizez/2+1 must be equal to the size of the wedge-volume");
        }
        reduced = true;
    }

    this->init_wedge_volume(sizex, sizey, sizez, reduced);

    return reduced ? this->wedge_reduced : this->wedge_full;
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::VolumeWedge<T>::VolumeWedge():
    w(NULL),
    wedge_full(NULL),
    wedge_reduced(NULL),
    current_phi(0),
    current_psi(0),
    current_theta(0),
    phi(0),
    psi(0),
    theta(0) {

    this->dims[0] = 0;
    this->dims[1] = 0;
    this->dims[2] = 0;
}

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T> template<typename T2>
tom::VolumeWedge<T>::VolumeWedge(const tom::Volume<T2> &w, bool copy):
    w(NULL),
    wedge_full(NULL),
    wedge_reduced(NULL),
    current_phi(0),
    current_psi(0),
    current_theta(0),
    phi(0),
    psi(0),
    theta(0) {

    this->dims[0] = w.getSizeX();
    this->dims[1] = w.getSizeY();
    this->dims[2] = w.getSizeZ();

    T2 w_min, w_max;
    w.minmax(w_min, w_max);

    // The wedge must contain only positive values.
    if (w_min < 0) {
        throw std::invalid_argument("Wedge can not have negative values.");
    }
    if (w_max == 0) {
        throw std::invalid_argument("The entire wedge is zero.");
    }
    this->active = !(w == w_max);
    if (this->active) {
        std::auto_ptr<tom::Volume<T> > autov(new tom::Volume<T>(this->dims[0], this->dims[1], this->dims[2], NULL, NULL));
        autov->setValues(w);
        // Scale so that the used wedge has always values between 0..1
        autov->template shift_scale<double>(0, 1./static_cast<double>(w_max));
        this->w = autov.release();
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T> template<typename T2>
tom::VolumeWedge<T>::VolumeWedge(tom::Volume<T2> &w, bool copy):
    w(NULL),
    wedge_full(NULL),
    wedge_reduced(NULL),
    current_phi(0),
    current_psi(0),
    current_theta(0),
    phi(0),
    psi(0),
    theta(0) {

    this->dims[0] = w.getSizeX();
    this->dims[1] = w.getSizeY();
    this->dims[2] = w.getSizeZ();

    this->active = !(w == 1);

    if (this->active) {
        std::auto_ptr<tom::Volume<T> > autov;
        if (copy || typeid(T) != typeid(T2)) {
            autov.reset(new tom::Volume<T>(this->dims[0], this->dims[1], this->dims[2], NULL, NULL));
            autov->setValues(w);
        } else {
            autov.reset(new tom::Volume<T>(w, NULL, this->dims[0], this->dims[1], this->dims[2], w.getStrideX(), w.getStrideY(), w.getStrideZ()));
        }
        this->w = autov.release();
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
bool tom::VolumeWedge<T>::is_active() const {
    return this->active;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::VolumeWedge<T>::rotate(double phi, double psi, double theta) {
    this->phi = phi;
    this->psi = psi;
    this->theta = theta;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::VolumeWedge<T>::get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const {

    if (!this->active) {
        return NULL;
    }

    if (sizex!=this->dims[0] || sizey!=this->dims[1] || sizez!=this->dims[2] || (size_z_vol!=sizez && size_z_vol!=sizez/2+1)) {
        throw std::invalid_argument("The size of the wedge does not match the saved volume");
    }

    bool reduced = sizez!=size_z_vol;

    this->init_wedge_volume();

    if (reduced) {
        return this->wedge_reduced;
    } else {
        return this->wedge_full;
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
bool tom::VolumeWedge<T>::apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez) const {

    if (vsrc.getSizeX()!=this->dims[0] || vsrc.getSizeY()!=this->dims[1] || sizez!=this->dims[2]) {
        throw std::invalid_argument("The size of the wedge does not match the saved volume");
    }

    bool reduced = false;
    if (sizez != vsrc.getSizeZ()) {
        if (sizez/2+1 != vsrc.getSizeZ()) {
            throw std::invalid_argument("either sizeZ or sizeZ/2+1 must be equal to the size of the volume");
        }
        reduced = true;
    }

    if (!this->active) {
        return false;
    }

    this->init_wedge_volume();

    if (reduced) {
        tom::element_wise_multiply<std::complex<T>, T>(vsrc, *this->wedge_reduced);
    } else {
        tom::element_wise_multiply<std::complex<T>, T>(vsrc, *this->wedge_full);
    }

    return true;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::VolumeWedge<T>::VolumeWedge(const tom::VolumeWedge<T> &v) {
    assert(0);
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::VolumeEMWedge<T>::VolumeEMWedge(const tom::VolumeEMWedge<T> &v) {
    assert(0);
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::VolumeWedge<T>::init_wedge_volume() const {


    assert(this->active && this->w);

    if (this->wedge_full && this->current_phi==this->phi && this->current_psi==this->psi && this->current_theta==this->theta) {

    } else {
        if (!this->wedge_full) {
            assert(!this->wedge_reduced);
            this->wedge_full = new tom::Volume<T>(this->dims[0], this->dims[1], this->dims[2], NULL, NULL);
            this->wedge_reduced = new tom::Volume<T>(*this->wedge_full, NULL, this->dims[0], this->dims[1], this->dims[2]/2+1, this->wedge_full->getStrideX(), this->wedge_full->getStrideY(), this->wedge_full->getStrideZ());
        }

        tom::Volume<T> v(this->dims[0], this->dims[1], this->dims[2], NULL, NULL);
        tom::rotate(*this->w, v, this->phi, this->psi, this->theta);
        tom::fftshift(v, *this->wedge_full, true);

        this->current_phi = phi;
        this->current_psi = psi;
        this->current_theta = theta;

    }
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::VolumeEMWedge<T>::VolumeEMWedge(const std::string &filename, std::size_t binning_)
    : VolumeWedge<T>(), filename(filename), binning(binning_) {

    if (this->binning < 1) { this->binning = 1; }

    std::auto_ptr<tom::Volume<T> > v;
    {
        const uint32_t vbinning[3] = { this->binning, this->binning, this->binning };
        tom::Volume<T> *pvol;
        tom::read_from_em<T>(pvol, filename, NULL,NULL,vbinning, NULL, NULL,NULL);
        v.reset(pvol);
    }

    this->dims[0] = v->getSizeX();
    this->dims[1] = v->getSizeY();
    this->dims[2] = v->getSizeZ();

    this->active = !(*v == 1);

    if (this->active) {
        this->w = v.release();
    }

}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::auto_ptr<tom::Wedge<T> > tom::WedgeDescriptor_NoWedge<T>::createWedge() const {
    std::auto_ptr<tom::Wedge<T> > res(new tom::NoWedge<T>());
    return res;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::auto_ptr<tom::Wedge<T> > tom::WedgeDescriptor_SimpleWedge<T>::createWedge() const {
    std::auto_ptr<tom::Wedge<T> > res(new tom::SimpleWedge<T>(this->angle, this->cutoff_radius));
    return res;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::auto_ptr<tom::Wedge<T> > tom::WedgeDescriptor_EMWedge<T>::createWedge() const {
    std::auto_ptr<tom::Wedge<T> > res(new tom::VolumeEMWedge<T>(this->filename, this->binning));
    return res;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::WedgeDescriptor<T> *tom::WedgeDescriptor<T>::fromConfigEntry(const std::string s, std::size_t binning) {

    std::istringstream iss(s);
    std::string s1, s2, s3;
    iss >> std::ws >> s1 >> std::ws >> s2;

    std::auto_ptr<tom::WedgeDescriptor<T> > res;

    if (!binning) { binning = 1; }


    //std::cout << std::setw(10) << linecnt << ": \"" << s1 << "\" \"" << s2 << "\"";

    if (s1 == "" || s1 == "nowedge") {
        if (!s2.empty()) {
            throw std::runtime_error("nowedge can not have other parameters. (WedgeDescriptor::fromConfigEntry)");
        }
        res.reset(new tom::WedgeDescriptor_NoWedge<T>());
    } else if (s1 == "emfile") {
        if (s2.size() > 2 && s2[0]=='"' && s2[s2.size()-1]=='"') {
            s2 = s2.substr(1, s2.size()-2);
        }
        if (s2.empty()) {
            throw std::runtime_error("emfile wedge expects the filename as parameter. (WedgeDescriptor::fromConfigEntry)");
        }
        res.reset(new tom::WedgeDescriptor_EMWedge<T>(s2, binning));
    } else if (s1 == "simple") {
        double angle, cutoff_radius;
        iss >> std::ws >> s3;
        try {
            angle = boost::lexical_cast<double>(s2);
            cutoff_radius = boost::lexical_cast<double>(s3);
        } catch (boost::bad_lexical_cast &e) {
            throw std::runtime_error("simple wedge expects two floating point arguments (i.e. the angle in DEG and the cutoff-radius). (WedgeDescriptor::fromConfigEntry)");
        }
        res.reset(new tom::WedgeDescriptor_SimpleWedge<T>(angle*0.01745329251994329576913914624236578987393, cutoff_radius, binning));
    } else {
        throw std::runtime_error("Unknown desciption for wedge (\"" + s + "\"). (WedgeDescriptor::fromConfigEntry)");
    }
    return res.release();
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::WedgeDescriptor_NoWedge<T>::getConfigEntry() const {
    return "nowedge";
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::WedgeDescriptor_SimpleWedge<T>::getConfigEntry() const {
    std::ostringstream ss;
    ss << "simple " << (this->angle*57.29577951308232087665461840231273527024) << " " << this->cutoff_radius_nobinning;
    return ss.str();
}

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::WedgeDescriptor_EMWedge<T>::getConfigEntry() const {
    std::ostringstream ss;
    ss << "emfile  \"" << this->filename << "\"";
    return ss.str();
}



// template instantiations.
template class tom::Wedge<float >;
template class tom::Wedge<double>;
template class tom::NoWedge<float >;
template class tom::NoWedge<double>;
template class tom::SimpleWedge<float >;
template class tom::SimpleWedge<double>;
template class tom::VolumeWedge<float >;
template class tom::VolumeWedge<double>;
template class tom::VolumeEMWedge<float >;
template class tom::VolumeEMWedge<double>;

template tom::VolumeWedge<float >::VolumeWedge(tom::Volume<float > &w, bool copy);
template tom::VolumeWedge<float >::VolumeWedge(tom::Volume<double> &w, bool copy);
template tom::VolumeWedge<double>::VolumeWedge(tom::Volume<float > &w, bool copy);
template tom::VolumeWedge<double>::VolumeWedge(tom::Volume<double> &w, bool copy);

template tom::VolumeWedge<float >::VolumeWedge(const tom::Volume<float > &w, bool copy);
template tom::VolumeWedge<float >::VolumeWedge(const tom::Volume<double> &w, bool copy);
template tom::VolumeWedge<double>::VolumeWedge(const tom::Volume<float > &w, bool copy);
template tom::VolumeWedge<double>::VolumeWedge(const tom::Volume<double> &w, bool copy);

template class tom::WedgeDescriptor_NoWedge<float >;
template class tom::WedgeDescriptor_NoWedge<double>;
template class tom::WedgeDescriptor_SimpleWedge<float >;
template class tom::WedgeDescriptor_SimpleWedge<double>;
template class tom::WedgeDescriptor_EMWedge<float >;
template class tom::WedgeDescriptor_EMWedge<double>;

template tom::WedgeDescriptor<float > *tom::WedgeDescriptor<float >::fromConfigEntry(const std::string s, std::size_t binning);
template tom::WedgeDescriptor<double> *tom::WedgeDescriptor<double>::fromConfigEntry(const std::string s, std::size_t binning);






