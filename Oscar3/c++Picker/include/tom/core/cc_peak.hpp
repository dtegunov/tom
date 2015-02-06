/****************************************************************************//**
 * \file cc_peak.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    05.02.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORR__CC_PEAK_HPP__
#define ___INCLUDE_CORR__CC_PEAK_HPP__







namespace tom {

/****************************************************************************//**
 * Structure to save the data of a 3D-correlation peak. Namely its shift,
 * the angle and the correlation value.
 *******************************************************************************/
template<typename T>
class cc_peak {
public:
    typedef double idx_type; // make it double here, because in case of binned volumes, the peak can be a non integer number.
    typedef int32_t angle_idx_type;
    typedef T val_type;

    typename cc_peak<T>::idx_type x, y, z;
    typename cc_peak<T>::angle_idx_type angle_idx;
    typename cc_peak<T>::val_type val;
    cc_peak(): x(0), y(0), z(0), angle_idx(-1), val(0) {};
};


}





#endif


