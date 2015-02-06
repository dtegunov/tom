/****************************************************************************//**
 * \file volume_loop.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    17.01.2008
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__VOLUME_LOOP_HPP__
#define ___INCLUDE_CORE__VOLUME_LOOP_HPP__





namespace tom {
namespace loop {

template<typename TPTR, typename TVOL, typename TOP> void for_each      (TVOL v, TOP local_OPERATOR);
template<typename TPTR, typename TVOL, typename TOP> void for_each_while(TVOL v, TOP local_OPERATOR);
template<typename TPTR, typename TVOL, typename TOP> void for_each_step (TVOL v, TOP local_OPERATOR);

template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TOP> void for_each(TVOL1 v1, TVOL2 v2, TOP local_OPERATOR);
template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TOP> void for_each_step(TVOL1 v1, TVOL2 v2, TOP local_OPERATOR);
template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TOP> void for_each_while(TVOL1 v1, TVOL2 v2, TOP local_OPERATOR);

template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TPTR3, typename TVOL3, typename TOP> void for_each(TVOL1 v1, TVOL2 v2, TVOL3 v3, TOP local_OPERATOR);


template<typename T> void ptr_add_byte_offset(T *&ptr, std::ptrdiff_t offset);




}
}






// INLINE FUNCTIONS




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline void tom::loop::ptr_add_byte_offset(T *&ptr, std::ptrdiff_t offset) {
    ptr = (T *)((char *)ptr + offset);
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR, typename TVOL, typename TOP>
inline void tom::loop::for_each(TVOL v, TOP local_OPERATOR) {
    typedef TPTR local_TYPE;
    local_TYPE *local_A = &v.get();
    std::size_t local_sizex = v.getSizeX();
    std::size_t local_sizey = v.getSizeY();
    std::size_t local_sizez = v.getSizeZ();
    std::size_t local_stridex = v.getStrideX();
    std::size_t local_stridey = v.getStrideY();
    std::size_t local_stridez = v.getStrideZ();

    {
        std::size_t x, y, z;

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
                    for (z=0; z<local_sizez; z++) {
                        for (y=0; y<local_sizey; y++) {
                            for (x=0; x<local_sizex; x++, i++) {
                                local_OPERATOR(local_A[i], x, y, z);
                            }
                        }
                    }
                } else {
                    local_stridez -= local_sizey*local_stridey;
                    for (z=0; z<local_sizez; z++) {
                        for (y=0; y<local_sizey; y++) {
                            for (x=0; x<local_sizex; x++) {
                                local_OPERATOR(local_A[x], x, y, z);
                            }
                            local_A += local_stridey;
                        }
                        local_A += local_stridez;
                    }
                }
            } else {
                local_stridez -= local_sizey*local_stridey;
                local_stridey -= local_sizex*local_stridex;
                for (z=0; z<local_sizez; z++) {
                    for (y=0; y<local_sizey; y++) {
                        for (x=0; x<local_sizex; x++) {
                            local_OPERATOR(local_A[x], x, y, z);
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

            for (z=0; z<local_sizez; z++) {
                for (y=0; y<local_sizey; y++) {
                    for (x=0; x<local_sizex; x++) {
                        local_OPERATOR(local_A[x], x, y, z);
                        tom::loop::ptr_add_byte_offset(local_A, local_stridex);
                    }
                    tom::loop::ptr_add_byte_offset(local_A, local_stridey);
                }
                tom::loop::ptr_add_byte_offset(local_A, local_stridez);
            }
        }
    }
}





/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR, typename TVOL, typename TOP>
inline void tom::loop::for_each_step(TVOL v, TOP local_OPERATOR) {
    typedef TPTR local_TYPE;
    local_TYPE *local_A = &v.get();
    std::size_t local_sizex = v.getSizeX();
    std::size_t local_sizey = v.getSizeY();
    std::size_t local_sizez = v.getSizeZ();
    std::size_t local_stridex = v.getStrideX();
    std::size_t local_stridey = v.getStrideY();
    std::size_t local_stridez = v.getStrideZ();

    {
        std::size_t x, y, z;

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
                    for (z=0; z<local_sizez; z++) {
                        local_OPERATOR.stepz(z);
                        for (y=0; y<local_sizey; y++) {
                            local_OPERATOR.stepy(y);
                            for (x=0; x<local_sizex; x++, i++) {
                                local_OPERATOR(local_A[i], x);
                            }
                        }
                    }
                } else {
                    local_stridez -= local_sizey*local_stridey;
                    for (z=0; z<local_sizez; z++) {
                        local_OPERATOR.stepz(z);
                        for (y=0; y<local_sizey; y++) {
                            local_OPERATOR.stepy(y);
                            for (x=0; x<local_sizex; x++) {
                                local_OPERATOR(local_A[x], x);
                            }
                            local_A += local_stridey;
                        }
                        local_A += local_stridez;
                    }
                }
            } else {
                local_stridez -= local_sizey*local_stridey;
                local_stridey -= local_sizex*local_stridex;
                for (z=0; z<local_sizez; z++) {
                    local_OPERATOR.stepz(z);
                    for (y=0; y<local_sizey; y++) {
                        local_OPERATOR.stepy(y);
                        for (x=0; x<local_sizex; x++) {
                            local_OPERATOR(local_A[x], x);
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

            for (z=0; z<local_sizez; z++) {
                local_OPERATOR.stepz(z);
                for (y=0; y<local_sizey; y++) {
                    local_OPERATOR.stepy(y);
                    for (x=0; x<local_sizex; x++) {
                        local_OPERATOR(local_A[x], x);
                        tom::loop::ptr_add_byte_offset(local_A, local_stridex);
                    }
                    tom::loop::ptr_add_byte_offset(local_A, local_stridey);
                }
                tom::loop::ptr_add_byte_offset(local_A, local_stridez);
            }
        }
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR, typename TVOL, typename TOP>
inline void tom::loop::for_each_while(TVOL v, TOP local_OPERATOR) {
    typedef TPTR local_TYPE;
    local_TYPE *local_A = &v.get();
    std::size_t local_sizex = v.getSizeX();
    std::size_t local_sizey = v.getSizeY();
    std::size_t local_sizez = v.getSizeZ();
    std::size_t local_stridex = v.getStrideX();
    std::size_t local_stridey = v.getStrideY();
    std::size_t local_stridez = v.getStrideZ();

    {
        std::size_t x, y, z;

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
                    for (z=0; z<local_sizez; z++) {
                        for (y=0; y<local_sizey; y++) {
                            for (x=0; x<local_sizex; x++, i++) {
                                if (!local_OPERATOR(local_A[i], x, y, z)) {
                                    return;
                                }
                            }
                        }
                    }
                } else {
                    local_stridez -= local_sizey*local_stridey;
                    for (z=0; z<local_sizez; z++) {
                        for (y=0; y<local_sizey; y++) {
                            for (x=0; x<local_sizex; x++) {
                                if (!local_OPERATOR(local_A[x], x, y, z)) {
                                    return;
                                }
                            }
                            local_A += local_stridey;
                        }
                        local_A += local_stridez;
                    }
                }
            } else {
                local_stridez -= local_sizey*local_stridey;
                local_stridey -= local_sizex*local_stridex;
                for (z=0; z<local_sizez; z++) {
                    for (y=0; y<local_sizey; y++) {
                        for (x=0; x<local_sizex; x++) {
                            if (!local_OPERATOR(local_A[x], x, y, z)) {
                                return;
                            }
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

            for (z=0; z<local_sizez; z++) {
                for (y=0; y<local_sizey; y++) {
                    for (x=0; x<local_sizex; x++) {
                        if (!local_OPERATOR(local_A[x], x, y, z)) {
                            return;
                        }
                        tom::loop::ptr_add_byte_offset(local_A, local_stridex);
                    }
                    tom::loop::ptr_add_byte_offset(local_A, local_stridey);
                }
                tom::loop::ptr_add_byte_offset(local_A, local_stridez);
            }
        }
    }
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TOP>
inline void tom::loop::for_each(TVOL1 v1, TVOL2 v2, TOP local_OPERATOR) {

    if (!v1.is_equal_size(v2)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }

    TPTR1 *ptr1 = &v1.get();
    TPTR2 *ptr2 = &v2.get();
    std::size_t sizex = v1.getSizeX();
    std::size_t sizey = v1.getSizeY();
    std::size_t sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    std::size_t stridex2 = v2.getStrideX();
    std::size_t stridey2 = v2.getStrideY();
    std::size_t stridez2 = v2.getStrideZ();

    std::size_t x, y, z;


    if (!(stridex1%sizeof(TPTR1)) &&
        !(stridey1%sizeof(TPTR1)) &&
        !(stridez1%sizeof(TPTR1)) &&
        !(stridex2%sizeof(TPTR2)) &&
        !(stridey2%sizeof(TPTR2)) &&
        !(stridez2%sizeof(TPTR2))) {

        stridex1 /= sizeof(TPTR1);
        stridey1 /= sizeof(TPTR1);
        stridez1 /= sizeof(TPTR1);
        stridex2 /= sizeof(TPTR2);
        stridey2 /= sizeof(TPTR2);
        stridez2 /= sizeof(TPTR2);

        if (stridex1 == 1 && stridex2 == 1) {
            if (stridey1 == sizex &&
                stridez1 == sizex*sizey &&
                stridey1 == stridey2 &&
                stridez1 == stridez2) {
                std::size_t i = 0;
                for (z=0; z<sizez; z++) {
                    for (y=0; y<sizey; y++) {
                        for (x=0; x<sizex; x++, i++) {
                            local_OPERATOR(ptr1[i], ptr2[i], x, y, z);
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                stridez2 -= sizey*stridey2;

                for (z=0; z<sizez; z++) {
                    for (y=0; y<sizey; y++) {
                        for (x=0; x<sizex; x++) {
                            local_OPERATOR(ptr1[x], ptr2[x], x, y, z);
                        }
                        ptr1 += stridey1;
                        ptr2 += stridey2;
                    }
                    ptr1 += stridez1;
                    ptr2 += stridez2;
                }
            }
        } else {
            stridez1 -= sizey*stridey1;
            stridez2 -= sizey*stridey2;
            stridey1 -= sizex*stridex1;
            stridey2 -= sizex*stridex2;
            for (z=0; z<sizez; z++) {
                for (y=0; y<sizey; y++) {
                    for (x=0; x<sizex; x++) {
                        local_OPERATOR(ptr1[x], ptr2[x], x, y, z);
                        ptr1 += stridex1;
                        ptr2 += stridex2;
                    }
                    ptr1 += stridey1;
                    ptr2 += stridey2;
                }
                ptr1 += stridez1;
                ptr2 += stridez2;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridez2 -= sizey*stridey2;
        stridey1 -= sizex*stridex1;
        stridey2 -= sizex*stridex2;
        for (z=0; z<sizez; z++) {
            for (y=0; y<sizey; y++) {
                for (x=0; x<sizex; x++) {
                    local_OPERATOR(ptr1[x], ptr2[x], x, y, z);
                    tom::loop::ptr_add_byte_offset(ptr1, stridex1);
                    tom::loop::ptr_add_byte_offset(ptr2, stridex2);
                }
                tom::loop::ptr_add_byte_offset(ptr1, stridey1);
                tom::loop::ptr_add_byte_offset(ptr2, stridey2);
            }
            tom::loop::ptr_add_byte_offset(ptr1, stridez1);
            tom::loop::ptr_add_byte_offset(ptr2, stridez2);
        }
    }
}


#include <iostream>


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TOP>
inline void tom::loop::for_each_step(TVOL1 v1, TVOL2 v2, TOP local_OPERATOR) {

    if (!v1.is_equal_size(v2)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }

    TPTR1 *ptr1 = &v1.get();
    TPTR2 *ptr2 = &v2.get();
    std::size_t sizex = v1.getSizeX();
    std::size_t sizey = v1.getSizeY();
    std::size_t sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    std::size_t stridex2 = v2.getStrideX();
    std::size_t stridey2 = v2.getStrideY();
    std::size_t stridez2 = v2.getStrideZ();

    std::size_t x, y, z;


    if (!(stridex1%sizeof(TPTR1)) &&
        !(stridey1%sizeof(TPTR1)) &&
        !(stridez1%sizeof(TPTR1)) &&
        !(stridex2%sizeof(TPTR2)) &&
        !(stridey2%sizeof(TPTR2)) &&
        !(stridez2%sizeof(TPTR2))) {

        stridex1 /= sizeof(TPTR1);
        stridey1 /= sizeof(TPTR1);
        stridez1 /= sizeof(TPTR1);
        stridex2 /= sizeof(TPTR2);
        stridey2 /= sizeof(TPTR2);
        stridez2 /= sizeof(TPTR2);

        if (stridex1 == 1 && stridex2 == 1) {
            if (stridey1 == sizex &&
                stridez1 == sizex*sizey &&
                stridey1 == stridey2 &&
                stridez1 == stridez2) {
                std::size_t i = 0;
                for (z=0; z<sizez; z++) {
                    local_OPERATOR.stepz(z);
                    for (y=0; y<sizey; y++) {
                        local_OPERATOR.stepy(y);
                        for (x=0; x<sizex; x++, i++) {
                            local_OPERATOR(ptr1[i], ptr2[i], x);
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                stridez2 -= sizey*stridey2;

                for (z=0; z<sizez; z++) {
                    local_OPERATOR.stepz(z);
                    for (y=0; y<sizey; y++) {
                        local_OPERATOR.stepy(y);
                        for (x=0; x<sizex; x++) {
                            local_OPERATOR(ptr1[x], ptr2[x], x);
                        }
                        ptr1 += stridey1;
                        ptr2 += stridey2;
                    }
                    ptr1 += stridez1;
                    ptr2 += stridez2;
                }
            }
        } else {
            stridez1 -= sizey*stridey1;
            stridez2 -= sizey*stridey2;
            stridey1 -= sizex*stridex1;
            stridey2 -= sizex*stridex2;
            for (z=0; z<sizez; z++) {
                local_OPERATOR.stepy(z);
                for (y=0; y<sizey; y++) {
                    local_OPERATOR.stepy(y);
                    for (x=0; x<sizex; x++) {
                        local_OPERATOR(ptr1[x], ptr2[x], x);
                        ptr1 += stridex1;
                        ptr2 += stridex2;
                    }
                    ptr1 += stridey1;
                    ptr2 += stridey2;
                }
                ptr1 += stridez1;
                ptr2 += stridez2;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridez2 -= sizey*stridey2;
        stridey1 -= sizex*stridex1;
        stridey2 -= sizex*stridex2;
        for (z=0; z<sizez; z++) {
            local_OPERATOR.stepz(z);
            for (y=0; y<sizey; y++) {
                local_OPERATOR.stepy(y);
                for (x=0; x<sizex; x++) {
                    local_OPERATOR(ptr1[x], ptr2[x], x);
                    tom::loop::ptr_add_byte_offset(ptr1, stridex1);
                    tom::loop::ptr_add_byte_offset(ptr2, stridex2);
                }
                tom::loop::ptr_add_byte_offset(ptr1, stridey1);
                tom::loop::ptr_add_byte_offset(ptr2, stridey2);
            }
            tom::loop::ptr_add_byte_offset(ptr1, stridez1);
            tom::loop::ptr_add_byte_offset(ptr2, stridez2);
        }
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TOP>
inline void tom::loop::for_each_while(TVOL1 v1, TVOL2 v2, TOP local_OPERATOR) {

    if (!v1.is_equal_size(v2)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }

    TPTR1 *ptr1 = &v1.get();
    TPTR2 *ptr2 = &v2.get();
    std::size_t sizex = v1.getSizeX();
    std::size_t sizey = v1.getSizeY();
    std::size_t sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    std::size_t stridex2 = v2.getStrideX();
    std::size_t stridey2 = v2.getStrideY();
    std::size_t stridez2 = v2.getStrideZ();

    std::size_t x, y, z;

    if (!(stridex1%sizeof(TPTR1)) &&
        !(stridey1%sizeof(TPTR1)) &&
        !(stridez1%sizeof(TPTR1)) &&
        !(stridex2%sizeof(TPTR2)) &&
        !(stridey2%sizeof(TPTR2)) &&
        !(stridez2%sizeof(TPTR2))) {

        stridex1 /= sizeof(TPTR1);
        stridey1 /= sizeof(TPTR1);
        stridez1 /= sizeof(TPTR1);
        stridex2 /= sizeof(TPTR2);
        stridey2 /= sizeof(TPTR2);
        stridez2 /= sizeof(TPTR2);

        if (stridex1 == 1 && stridex2 == 1) {
            if (stridey1 == sizex &&
                stridez1 == sizex*sizey &&
                stridey1 == stridey2 &&
                stridez1 == stridez2) {
                std::size_t i = 0;
                for (z=0; z<sizez; z++) {
                    for (y=0; y<sizey; y++) {
                        for (x=0; x<sizex; x++, i++) {
                            if (!local_OPERATOR(ptr1[i], ptr2[i], x, y, z)) {
                                return;
                            }
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                stridez2 -= sizey*stridey2;
                for (z=0; z<sizez; z++) {
                    for (y=0; y<sizey; y++) {
                        for (x=0; x<sizex; x++) {
                            if (!local_OPERATOR(ptr1[x], ptr2[x], x, y, z)) {
                                return;
                            }
                        }
                        ptr1 += stridey2;
                        ptr2 += stridey2;
                    }
                    ptr1 += stridez1;
                    ptr2 += stridez2;
                }
            }
        } else {
            stridez1 -= sizey*stridey1;
            stridez2 -= sizey*stridey2;
            stridey1 -= sizex*stridex1;
            stridey2 -= sizex*stridex2;
            for (z=0; z<sizez; z++) {
                for (y=0; y<sizey; y++) {
                    for (x=0; x<sizex; x++) {
                        if (!local_OPERATOR(ptr1[x], ptr2[x], x, y, z)) {
                            return;
                        }
                        ptr1 += stridex1;
                        ptr2 += stridex2;
                    }
                    ptr1 += stridey1;
                    ptr2 += stridey2;
                }
                ptr1 += stridez1;
                ptr2 += stridez2;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridez2 -= sizey*stridey2;
        stridey1 -= sizex*stridex1;
        stridey2 -= sizex*stridex2;
        for (z=0; z<sizez; z++) {
            for (y=0; y<sizey; y++) {
                for (x=0; x<sizex; x++) {
                    if (!local_OPERATOR(ptr1[x], ptr2[x], x, y, z)) {
                        return;
                    }
                    tom::loop::ptr_add_byte_offset(ptr1, stridex1);
                    tom::loop::ptr_add_byte_offset(ptr2, stridex2);
                }
                tom::loop::ptr_add_byte_offset(ptr1, stridey1);
                tom::loop::ptr_add_byte_offset(ptr2, stridey2);
            }
            tom::loop::ptr_add_byte_offset(ptr1, stridez1);
            tom::loop::ptr_add_byte_offset(ptr2, stridez2);
        }
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename TPTR1, typename TVOL1, typename TPTR2, typename TVOL2, typename TPTR3, typename TVOL3, typename TOP>
void tom::loop::for_each(TVOL1 v1, TVOL2 v2, TVOL3 v3, TOP local_OPERATOR) {

    if (!v1.is_equal_size(v2) || !v1.is_equal_size(v3)) {
        throw std::invalid_argument("All volumes must have the same size.");
    }

    TPTR1 *ptr1 = &v1.get();
    TPTR2 *ptr2 = &v2.get();
    TPTR3 *ptr3 = &v3.get();
    std::size_t sizex = v1.getSizeX();
    std::size_t sizey = v1.getSizeY();
    std::size_t sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    std::size_t stridex2 = v2.getStrideX();
    std::size_t stridey2 = v2.getStrideY();
    std::size_t stridez2 = v2.getStrideZ();
    std::size_t stridex3 = v3.getStrideX();
    std::size_t stridey3 = v3.getStrideY();
    std::size_t stridez3 = v3.getStrideZ();

    std::size_t x, y, z;

    if (!(stridex1%sizeof(TPTR1)) &&
        !(stridey1%sizeof(TPTR1)) &&
        !(stridez1%sizeof(TPTR1)) &&
        !(stridex2%sizeof(TPTR2)) &&
        !(stridey2%sizeof(TPTR2)) &&
        !(stridez2%sizeof(TPTR2)) &&
        !(stridex3%sizeof(TPTR3)) &&
        !(stridey3%sizeof(TPTR3)) &&
        !(stridez3%sizeof(TPTR3))) {

        stridex1 /= sizeof(TPTR1);
        stridey1 /= sizeof(TPTR1);
        stridez1 /= sizeof(TPTR1);
        stridex2 /= sizeof(TPTR2);
        stridey2 /= sizeof(TPTR2);
        stridez2 /= sizeof(TPTR2);
        stridex3 /= sizeof(TPTR3);
        stridey3 /= sizeof(TPTR3);
        stridez3 /= sizeof(TPTR3);

        if (stridex1 == 1 && stridex2 == 1 && stridex3 == 1) {
            if (stridey1 == sizex &&
                stridez1 == sizex*sizey &&
                stridey1 == stridey2 &&
                stridez1 == stridez2 &&
                stridey1 == stridey3 &&
                stridez1 == stridez3) {
                std::size_t i = 0;
                for (z=0; z<sizez; z++) {
                    for (y=0; y<sizey; y++) {
                        for (x=0; x<sizex; x++, i++) {
                            local_OPERATOR(ptr1[i], ptr2[i], ptr3[i], x, y, z);
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                stridez2 -= sizey*stridey2;
                stridez3 -= sizey*stridey3;
                for (z=0; z<sizez; z++) {
                    for (y=0; y<sizey; y++) {
                        for (x=0; x<sizex; x++) {
                            local_OPERATOR(ptr1[x], ptr2[x], ptr3[x], x, y, z);
                        }
                        ptr1 += stridey1;
                        ptr2 += stridey2;
                        ptr3 += stridey3;
                    }
                    ptr1 += stridez1;
                    ptr2 += stridez2;
                    ptr3 += stridez3;
                }
            }
        } else {
            stridez1 -= sizey*stridey1;
            stridez2 -= sizey*stridey2;
            stridez3 -= sizey*stridey3;
            stridey1 -= sizex*stridex1;
            stridey2 -= sizex*stridex2;
            stridey3 -= sizex*stridex3;
            for (z=0; z<sizez; z++) {
                for (y=0; y<sizey; y++) {
                    for (x=0; x<sizex; x++) {
                        local_OPERATOR(ptr1[x], ptr2[x], ptr3[x], x, y, z);
                        ptr1 += stridex1;
                        ptr2 += stridex2;
                        ptr3 += stridex3;
                    }
                    ptr1 += stridey1;
                    ptr2 += stridey2;
                    ptr3 += stridey3;
                }
                ptr1 += stridez1;
                ptr2 += stridez2;
                ptr3 += stridey3;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridez2 -= sizey*stridey2;
        stridez3 -= sizey*stridey3;
        stridey1 -= sizex*stridex1;
        stridey2 -= sizex*stridex2;
        stridey3 -= sizex*stridex3;
        for (z=0; z<sizez; z++) {
            for (y=0; y<sizey; y++) {
                for (x=0; x<sizex; x++) {
                    local_OPERATOR(ptr1[x], ptr2[x], ptr3[x], x, y, z);
                    tom::loop::ptr_add_byte_offset(ptr1, stridex1);
                    tom::loop::ptr_add_byte_offset(ptr2, stridex2);
                    tom::loop::ptr_add_byte_offset(ptr3, stridex3);
                }
                tom::loop::ptr_add_byte_offset(ptr1, stridey1);
                tom::loop::ptr_add_byte_offset(ptr2, stridey2);
                tom::loop::ptr_add_byte_offset(ptr3, stridey3);
            }
            tom::loop::ptr_add_byte_offset(ptr1, stridez1);
            tom::loop::ptr_add_byte_offset(ptr2, stridez2);
            tom::loop::ptr_add_byte_offset(ptr3, stridez3);
        }
    }
}



#endif


