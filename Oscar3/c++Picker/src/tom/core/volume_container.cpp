/***********************************************************************//**
 * \file volume_container.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    02.12.2007
 **************************************************************************/
#include "tom/core/volume_container.hpp"


#include <typeinfo>


/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
void tom::VolumeContainerEM<T>::loadVolume() const {

    if (this->v) {
        return;
    }

    std::auto_ptr<tom::Volume<T> > tmp;
    tom::Volume<T> *pvol;
    std::stringstream ss;

    try {
        tom::read_from_em<T>(pvol, this->filename, NULL, NULL, NULL, NULL, NULL, NULL);
        tmp.reset(pvol);
    } catch (int &i) {
        ss << "Error reading volume from em-file " << this->filename << " (" << i << ")";
        throw std::runtime_error(ss.str());
    }

    this->clearCache();
    this->v = tmp.release();
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T> template<typename T2>
tom::VolumeContainerMem<T>::VolumeContainerMem(std::string name, const tom::Volume<T2> &v, bool copy): name(name), v(NULL) {
    std::auto_ptr<tom::Volume<T> > autov(new tom::Volume<T>(v.getSizeX(), v.getSizeY(), v.getSizeZ(), NULL, NULL));
    autov->setValues(v);
    this->v = autov.get();
}


/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T> template<typename T2>
tom::VolumeContainerMem<T>::VolumeContainerMem(std::string name, tom::Volume<T2> &v, bool copy): name(name), v(NULL) {

    std::auto_ptr<tom::Volume<T> > autov;
    if (copy || typeid(T)!=typeid(T2)) {
        autov.reset(new tom::Volume<T>(v.getSizeX(), v.getSizeY(), v.getSizeZ(), NULL, NULL));
        autov->setValues(v);
    } else {
        autov.reset(new tom::Volume<T>(const_cast<tom::Volume<T2> &>(v), NULL, v.getSizeX(), v.getSizeY(), v.getSizeZ(), v.getStrideX(), v.getStrideY(), v.getStrideZ()));
    }
    this->v = autov.release();
}


/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
tom::VolumeContainerMem<T>::~VolumeContainerMem() {
    delete this->v;
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
std::string tom::VolumeContainerMem<T>::getName() const {
    return this->name;
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
void tom::VolumeContainerMem<T>::clearCache() const {
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
const tom::Volume<T> &tom::VolumeContainerMem<T>::getVolume() const {
    return *this->v;
}



/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
tom::VolumeContainerMem<T>::VolumeContainerMem(const VolumeContainerMem<T> &c) {
    throw std::runtime_error("Don't call the cctor");
}



// Template instantiations

template class tom::VolumeContainerEM<float >;
template class tom::VolumeContainerEM<double>;
template class tom::VolumeContainerMem<float >;
template class tom::VolumeContainerMem<double>;

template tom::VolumeContainerMem<float >::VolumeContainerMem(std::string name, tom::Volume<float > &v, bool copy);
template tom::VolumeContainerMem<float >::VolumeContainerMem(std::string name, tom::Volume<double> &v, bool copy);
template tom::VolumeContainerMem<double>::VolumeContainerMem(std::string name, tom::Volume<float > &v, bool copy);
template tom::VolumeContainerMem<double>::VolumeContainerMem(std::string name, tom::Volume<double> &v, bool copy);

template tom::VolumeContainerMem<float >::VolumeContainerMem(std::string name, const tom::Volume<float > &v, bool copy);
template tom::VolumeContainerMem<float >::VolumeContainerMem(std::string name, const tom::Volume<double> &v, bool copy);
template tom::VolumeContainerMem<double>::VolumeContainerMem(std::string name, const tom::Volume<float > &v, bool copy);
template tom::VolumeContainerMem<double>::VolumeContainerMem(std::string name, const tom::Volume<double> &v, bool copy);


