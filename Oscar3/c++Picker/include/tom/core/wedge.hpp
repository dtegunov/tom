/****************************************************************************//**
 * \file wedge.hpp
 * \brief The header file for wedge classes.
 * \author  Thomas Haller
 * \version 0.1
 * \date    02.12.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__WEDGE_HPP__
#define ___INCLUDE_CORE__WEDGE_HPP__




#include <complex>
#include <iosfwd>


#include "tom/core/defines.h"
#include "tom/core/volume.hpp"
#include "tom/core/fftw_plan.hpp"


namespace tom {



template<typename T> class WedgeDescriptor;
template<typename T> class WedgeDescriptor_NoWedge;
template<typename T> class WedgeDescriptor_SimpleWedge;
template<typename T> class WedgeDescriptor_EMWedge;


/****************************************************************************//**
 *
 *******************************************************************************/
template <typename T>
class Wedge {

public:

    typedef T real;
    typedef std::complex<T> complex;

    virtual ~Wedge() { };

    virtual std::auto_ptr<tom::WedgeDescriptor<T> > createWedgeDescriptor() const = 0;

    virtual void rotate(double phi, double psi, double theta) = 0;
    virtual bool is_active() const = 0;
    virtual const tom::Volume<T> *get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const = 0;

    /****************************************************************************//**
     * \brief Applay the wedge to the volume in fourier space.
     *
     * \param[in,out] vsrc The fourier transformed volume of complex elements.
     * \param[in] sizez The size of volume along Z. Since vsrc is the fourier transformed
     *   of a real volume, it has a hermitian symmetry. sizez==v.getSizeZ() means,
     *   that the volume vsrc is given full. If sizez/2+1==v.getSizeZ(), than only the
     *   upper half of the volume is given (as output by the fftw-r2c transformation).
     *   If size is not one of these two values, an exception is thrown.
     * \return A boolean volume, specifying whether vsrc was changed by the wedge.
     *   This can be usefull, to decide whether a normalisation done before has to be
     *   repeated or not.
     *******************************************************************************/
    virtual bool apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez ) const = 0;


    virtual std::string toString() const  = 0;

};





/****************************************************************************//**
 *
 *******************************************************************************/
template <typename T>
class NoWedge: public Wedge<T> {

public:

    NoWedge();
    virtual ~NoWedge();

    virtual std::auto_ptr<tom::WedgeDescriptor<T> > createWedgeDescriptor() const { return std::auto_ptr<tom::WedgeDescriptor<T> >(new tom::WedgeDescriptor_NoWedge<T>()); }

    virtual void rotate(double phi, double psi, double theta);
    virtual bool apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez) const;
    virtual bool is_active() const;
    virtual const tom::Volume<T> *get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const;


    virtual std::string toString() const { return "NoWedge"; };

};



/****************************************************************************//**
 *
 *******************************************************************************/
template <typename T>
class SimpleWedge: public Wedge<T> {

public:

    SimpleWedge(double angle, double cutoff_radius);
    virtual ~SimpleWedge();

    virtual std::auto_ptr<tom::WedgeDescriptor<T> > createWedgeDescriptor() const { return std::auto_ptr<tom::WedgeDescriptor<T> >(new tom::WedgeDescriptor_SimpleWedge<T>(this->angle, this->cutoff_radius, 0)); }

    virtual void rotate(double phi, double psi, double theta);
    virtual bool apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez) const;
    virtual bool is_active() const;
    virtual const tom::Volume<T> *get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const ;

    double getAngle() const         { return this->angle;           }
    double getCutoffRadius() const { return this->cutoff_radius;   }

    virtual std::string toString() const;


private:


    SimpleWedge(const SimpleWedge<T> &c);


    mutable tom::Volume<T> *wedge_reduced;
    mutable tom::Volume<T> *wedge_full;
    mutable double current_phi;
    mutable double current_psi;
    mutable double current_theta;

    double angle;
    double phi;
    double psi;
    double theta;
    double cutoff_radius;

    void init_wedge_volume(std::size_t sizex, std::size_t sizey, std::size_t sizez, bool reduced) const;
    void init_wedge_volume(tom::Volume<T> &v, std::size_t sizez) const;

};






/****************************************************************************//**
 *
 *******************************************************************************/
template <typename T>
class VolumeWedge: public Wedge<T> {

public:

    template<typename T2> VolumeWedge(tom::Volume<T2> &w, bool copy);
    template<typename T2> VolumeWedge(const tom::Volume<T2> &w, bool copy);
    virtual ~VolumeWedge();

    virtual std::auto_ptr<tom::WedgeDescriptor<T> > createWedgeDescriptor() const { throw std::runtime_error("The WedgeDescriptor for VolumeWedge is not yet implemented (as not needed)."); }

    virtual void rotate(double phi, double psi, double theta);
    virtual bool apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez) const;
    virtual bool is_active() const;
    virtual const tom::Volume<T> *get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const ;

    virtual std::string toString() const;

protected:
    tom::Volume<T> *w;
    std::size_t dims[3];
    VolumeWedge();
    bool active;


private:
    VolumeWedge(const VolumeWedge<T> &c);
    void init_wedge_volume() const;

    mutable tom::Volume<T> *wedge_full;
    mutable tom::Volume<T> *wedge_reduced;
    mutable double current_phi;
    mutable double current_psi;
    mutable double current_theta;

    double phi;
    double psi;
    double theta;


};

/****************************************************************************//**
 *
 *******************************************************************************/
template <typename T>
class VolumeEMWedge: public VolumeWedge<T> {

public:
    VolumeEMWedge(const std::string &filename, std::size_t binning);
    virtual ~VolumeEMWedge();
    virtual std::auto_ptr<tom::WedgeDescriptor<T> > createWedgeDescriptor() const { return std::auto_ptr<tom::WedgeDescriptor<T> >(new tom::WedgeDescriptor_EMWedge<T>(this->filename, this->binning)); }

    virtual std::string toString() const;
private:
    VolumeEMWedge(const VolumeEMWedge<T> &c);

    std::string filename;
    std::size_t binning;
};










/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
class WedgeDescriptor {
public:
    virtual std::auto_ptr<tom::Wedge<T> > createWedge() const = 0;
    virtual ~WedgeDescriptor() { };
    virtual int getTypeID() const = 0;

    virtual std::string getConfigEntry() const = 0;
    static WedgeDescriptor<T> *fromConfigEntry(const std::string s, std::size_t binning);
};




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
class WedgeDescriptor_NoWedge: public WedgeDescriptor<T> {
public:
    virtual std::auto_ptr<tom::Wedge<T> > createWedge() const;
    virtual ~WedgeDescriptor_NoWedge() { };
    static int getClassTypeID() { return __LINE__; }
    virtual int getTypeID() const { return this->getClassTypeID(); };
    virtual std::string getConfigEntry() const;
};




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
class WedgeDescriptor_SimpleWedge: public WedgeDescriptor<T> {
public:
    WedgeDescriptor_SimpleWedge(double angle, double cutoff_radius, std::size_t binning)
        : binning(binning?binning:1),
          cutoff_radius_nobinning(cutoff_radius),
          angle(angle),
          cutoff_radius(cutoff_radius/(binning?binning:1)) { }

    virtual std::auto_ptr<tom::Wedge<T> > createWedge() const;
    virtual ~WedgeDescriptor_SimpleWedge() { };
    double getAngle() const { return this->angle; }
    double getCutoff_radius() const { return this->cutoff_radius; }
    static int getClassTypeID() { return __LINE__; }
    virtual int getTypeID() const { return this->getClassTypeID(); };
    virtual std::string getConfigEntry() const;
private:
    std::size_t binning;
    double cutoff_radius_nobinning;
    double angle;
    double cutoff_radius;
};




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
class WedgeDescriptor_EMWedge: public WedgeDescriptor<T> {
public:
    WedgeDescriptor_EMWedge(const std::string &filename, std::size_t binning): filename(filename), binning(binning) { };
    virtual std::auto_ptr<tom::Wedge<T> > createWedge() const;
    virtual ~WedgeDescriptor_EMWedge() { };
    const std::string &getFilename() const { return this->filename; }
    std::size_t getBinning() const { return this->binning; }
    static int getClassTypeID() { return __LINE__; }
    virtual int getTypeID() const { return this->getClassTypeID(); };
    virtual std::string getConfigEntry() const;
private:
    const std::string filename;
    const std::size_t binning;
};


}


template<typename T>
std::istream &operator>>(std::istream &is, const tom::Wedge<T> &w) { is << w.toString(); return is; }






// Inline functions.

/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline tom::NoWedge<T>::NoWedge() {
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline tom::NoWedge<T>::~NoWedge() {
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline void tom::NoWedge<T>::rotate(double phi, double psi, double theta) {
    // Do nothing :)
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline bool tom::NoWedge<T>::apply(tom::Volume<typename Wedge<T>::complex> &vsrc, std::size_t sizez) const {
    return false;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline bool tom::NoWedge<T>::is_active() const {
    return false;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline const tom::Volume<T> *tom::NoWedge<T>::get_wedge(std::size_t sizex, std::size_t sizey, std::size_t size_z_vol, std::size_t sizez) const {
    if (size_z_vol != sizez && sizez/2+1 != size_z_vol) {
        throw std::invalid_argument("The size of the Z dimension of the wedge must be either the reduced form (sizez/2+1) or the full.");
    }
    return NULL;
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline tom::SimpleWedge<T>::~SimpleWedge() {
    delete this->wedge_reduced;
    delete this->wedge_full;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline bool tom::SimpleWedge<T>::is_active() const {
    return  this->angle>0 || this->cutoff_radius>0;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline tom::VolumeWedge<T>::~VolumeWedge() {
    delete this->w;
    delete this->wedge_full;
    delete this->wedge_reduced;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline tom::VolumeEMWedge<T>::~VolumeEMWedge() {
}


#endif









