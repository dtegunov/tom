/****************************************************************************//**
 * \file volume_container.hpp
 * \brief The header file for volume_container.
 * \author  Thomas Haller
 * \version 0.1
 * \date    11.12.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__VOLUME_CONTAINER_HPP__
#define ___INCLUDE_CORE__VOLUME_CONTAINER_HPP__


#include <complex>



#include "tom/core/volume.hpp"


namespace tom {


template <typename T>
class VolumeContainer {
public:

    virtual ~VolumeContainer() {};

    virtual std::string getName() const = 0;
    virtual void clearCache() const = 0;
    virtual const tom::Volume<T> &getVolume() const = 0;

};


template <typename T>
class VolumeContainerEM: public VolumeContainer<T> {

public:

    VolumeContainerEM(const std::string &filename);
    virtual ~VolumeContainerEM();

    virtual std::string getName() const;
    virtual void clearCache() const;
    virtual const tom::Volume<T> &getVolume() const;

private:
    VolumeContainerEM(const VolumeContainerEM<T> &c);

    void loadVolume() const;
    std::string filename;
    mutable tom::Volume<T> *v;
};




template <typename T>
class VolumeContainerMem: public VolumeContainer<T> {

public:

    template<typename T2> VolumeContainerMem(std::string name, const tom::Volume<T2> &v, bool copy);
    template<typename T2> VolumeContainerMem(std::string name, tom::Volume<T2> &v, bool copy);
    virtual ~VolumeContainerMem();

    virtual std::string getName() const;
    virtual void clearCache() const;
    virtual const tom::Volume<T> &getVolume() const;

private:
    VolumeContainerMem(const VolumeContainerMem<T> &c);
    std::string name;

    tom::Volume<T> *v;

};









}









// Inline functions.

/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
inline tom::VolumeContainerEM<T>::VolumeContainerEM(const std::string &filename): filename(filename), v(NULL) {
}

/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
inline tom::VolumeContainerEM<T>::~VolumeContainerEM() {
    delete this->v;
}

/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
inline void tom::VolumeContainerEM<T>::clearCache() const {
    delete this->v;
    this->v = NULL;
}

/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
inline const tom::Volume<T> &tom::VolumeContainerEM<T>::getVolume() const {
    this->loadVolume();
    return *this->v;
}


/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
inline tom::VolumeContainerEM<T>::VolumeContainerEM(const VolumeContainerEM<T> &c) {
    throw std::runtime_error("Don't call the cctor");
}

/****************************************************************************//**
 * \brief
 *******************************************************************************/
template<typename T>
inline std::string tom::VolumeContainerEM<T>::getName() const {
    return this->filename;
}

#endif
