/****************************************************************************//**
 * \file _transform.c
 * \author Thomas Haller
 * \date Oct. 26, 2007
 * \brief Contains implementations of tom_transform.c
 *
 * Include this C-file(!) with defining the names of the functions and the
 * datatype
 *******************************************************************************/




#undef FLOATTYPE
#undef tom_transf_interpol3d_nearest_double
#undef tom_transf_interpol3d_trilinear_double
#undef tom_transf_proj3d_double
#undef tom_transf_proj3d_solid_zxz_double
#undef tom_transf_rotate3d_double
#undef tom_transf_transf3d_double

#if FLOATTYPE_ISDOUBLE
    #define FLOATTYPE                                       double
#elif FLOATTYPE_ISFLOAT
    #define FLOATTYPE                                       float
    #define tom_transf_interpol3d_nearest_double            tom_transf_interpol3d_nearest_float
    #define tom_transf_interpol3d_trilinear_double          tom_transf_interpol3d_trilinear_float
    #define tom_transf_proj3d_double                        tom_transf_proj3d_float
    #define tom_transf_proj3d_solid_zxz_double              tom_transf_proj3d_solid_zxz_float
    #define tom_transf_rotate3d_double                      tom_transf_rotate3d_float
    #define tom_transf_transf3d_double                      tom_transf_transf3d_float
#else
    #error FLOATTYPE not appropriate defined.
#endif



/****************************************************************************//**
 * \brief Rotate a 3D volume.
 *
 * Rotates a Volume around a point in space with angles psi,phi,theta.
 * The form of the rotation is the same as in the TOM-toolbox for the
 * function tom_rotate. That means that it performs a rotation of the
 * coordinate system with the angles phi, psi and theta. The order of
 * the rotation is the following: first it rotates around the Z-axis with
 * angle phi, then around the X-Axis with angle theta and finally again
 * around Z with angle psi.
 *
 * \param[in]  src The input volume.
 * \param[out] dst A pointer to the output volume. It must be the same
 *   size and the same datatype as src (non overlapping).
 * \param[in] sizex The Number of elements along the X-direction.
 * \param[in] sizey The Number of elements along the Y-direction.
 * \param[in] sizez The Number of elements along the Z-direction.
 * \param[in] phi Angle to rotate around the Z-axis, in radians.
 * \param[in] psi Angle to rotate around the Z-axis, in radians.
 * \param[in] theta Angle to rotate around the X-axis, in radians.
 * \param[in] centerx The x-coordinate of the center of the rotation.
 * \param[in] centery The y-coordinate of the center of the rotation.
 * \param[in] centerz The z-coordinate of the center of the rotation.
 * \param[in] intfcn The interpolation function. If set to null, the
 *   standard interpolation function interpol3d_trilinear_GENERIC is used.
 * \param[in] optparam Parameters passed to the interpolation function.
 *   This parameter is ignored if intfcn==NULL. It depends entirely on
 *   intfcn, what is expected here.
 * \return 0 in case of success. This value is the same as returned
 *   from transf3d_GENERIC.
 *
 * This function is only a wrapper for tom_transf_transf3d_double where the in- and
 * output volume have the same dimension, are continous in memory and the
 * transformation is a rotation of the coordinate system around a point.
 *******************************************************************************/
int tom_transf_rotate3d_double(const FLOATTYPE *src, FLOATTYPE *dst, size_t sizex, size_t sizey, size_t sizez,
              double phi, double psi, double theta, double centerx, double centery, double centerz,
              FLOATTYPE (*intfcn)(const tom_transf_st_interpol_param *param), void *optparam) {

    double P[16] = { 0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 1. };
    int axes[3] = { 2, 0, 2 };
    double angles[3];
    FLOATTYPE defaultval;

    angles[0] = phi;
    angles[1] = theta;
    angles[2] = psi;
    tom_transf_init_rotmat_3d(P, 4, 3, axes, angles);

    tom_transf_invert(P, P, 4);

    P[ 3] = centerx - (P[ 0]*centerx + P[ 1]*centery + P[ 2]*centerz);
    P[ 7] = centery - (P[ 4]*centerx + P[ 5]*centery + P[ 6]*centerz);
    P[11] = centerz - (P[ 8]*centerx + P[ 9]*centery + P[10]*centerz);

    if (!intfcn) {
        /* Set default function. */
        intfcn = tom_transf_interpol3d_trilinear_double;
        defaultval = 0.;
        optparam = &defaultval;
    }

    return tom_transf_transf3d_double(  src, sizex, sizey, sizez, sizex, sizex*sizey,
                                        dst, sizex, sizey, sizez, sizex, sizex*sizey,
                                        P, 1, 0., intfcn, optparam);
}





/****************************************************************************//**
 * \brief One of the 3d-interpolation functions to do trilinear interpolation.
 * \param param A pointer to a structure tom_transf_st_interpol_param which
 *   is initialized by the calling function with the parameters of the current
 *   interpolation.
 * \returns The interpolated value.
 *******************************************************************************/
FLOATTYPE tom_transf_interpol3d_trilinear_double(const tom_transf_st_interpol_param *param) {

#if 0
    const size_t sizex_1 = param->sizex - 1;
    const size_t sizey_1 = param->sizey - 1;
    const size_t sizez_1 = param->sizez - 1;
    if (param->x < 0 || param->x > sizex_1 ||
        param->y < 0 || param->y > sizey_1 ||
        param->z < 0 || param->z > sizez_1) {
        return ((const FLOATTYPE *)param->optparam)[0];
    }

    {

        FLOATTYPE img0 = 1.;
        FLOATTYPE img1 = 1.;
        FLOATTYPE img2 = 1.;
        FLOATTYPE img3 = 1.;
        FLOATTYPE img4 = 1.;
        FLOATTYPE img5 = 1.;
        FLOATTYPE img6 = 1.;
        FLOATTYPE AA,BB,CC,DD;

        const size_t index1 = param->stepy;
        const size_t index2 = index1 + 1;
        const size_t index3 = param->stepz;
        const size_t index4 = index3 + 1;
        const size_t index5 = index1 + index3;
        const size_t index6 = index5 + 1;

        const FLOATTYPE x = param->x;
        const FLOATTYPE y = param->y;
        const FLOATTYPE z = param->z;

        const size_t floorx = (size_t)x;
        const size_t floory = (size_t)y;
        const size_t floorz = (size_t)z;

        /* Interpolation */
        const FLOATTYPE vx2 = x - floorx;
        const FLOATTYPE vy2 = y - floory;
        const FLOATTYPE vy1 = 1 - vy2;
        const FLOATTYPE vz2 = z - floorz;
        const FLOATTYPE vz1 = 1 - vz2;

        const FLOATTYPE *src = &((const FLOATTYPE *)param->src)[floorz * param->stepz + floory * param->stepy + floorx];

        /* the following section detects border pixels to avoid exceeding dimensions */
        if (!vx2) {
            img0 = img2 = img4 = img6 = 0.;
        }
        if (!vy2) {
            img1 = img2 = img5 = img6 = 0.;
        }
        if (!vz2) {
            img3 = img4 = img5 = img6 = 0.;
        }

        if (img0) {
            img0 = src[1];
        }
        if (img1) {
            img1 = src[index1];
        }
        if (img2) {
            img2 = src[index2];
        }
        if (img3) {
            img3 = src[index3];
        }
        if (img4) {
            img4 = src[index4];
        }
        if (img5) {
            img5 = src[index5];
        }
        if (img6) {
            img6 = src[index6];
        }

        /* interpolation */
        AA = src[0] + (img0 - src[0]) * vx2;
        BB = img1   + (img2 - img1)   * vx2;
        CC = img3   + (img4 - img3)   * vx2;
        DD = img5   + (img6 - img5)   * vx2;
        return  (AA * vy1 + BB * vy2) * vz1 + (CC * vy1 + DD * vy2) * vz2;
    }
#else
    if (param->x<0. || param->y<0. || param->z<0. ||
        param->x > (param->sizex-1) ||
        param->y > (param->sizey-1) ||
        param->z > (param->sizez-1)) {
        return ((const FLOATTYPE *)param->optparam)[0];
    }
    {
        const size_t floorx = (size_t)param->x;
        const size_t floory = (size_t)param->y;
        const size_t floorz = (size_t)param->z;
        const FLOATTYPE xoffseth = param->x - floorx;
        const FLOATTYPE yoffseth = param->y - floory;
        const FLOATTYPE zoffseth = param->z - floorz;

        const FLOATTYPE *src = &((const FLOATTYPE *)param->src)[floorz * param->stepz + floory * param->stepy + floorx];

        switch ((xoffseth > 0.) | ((yoffseth > 0.)<<1) | ((zoffseth > 0.)<<2)) {
            case 0x07: {
                    FLOATTYPE x00x, x01x, x10x, x11x, x0yx;
                    x00x = src[0] + (src[1]-src[0])*xoffseth;
                    src += param->stepy;
                    x01x = src[0] + (src[1]-src[0])*xoffseth;
                    src += param->stepz;
                    x11x = src[0] + (src[1]-src[0])*xoffseth;
                    src -= param->stepy;
                    x10x = src[0] + (src[1]-src[0])*xoffseth;
                    x0yx = x00x + (x01x - x00x) * yoffseth;
                    return x0yx + ((x10x + (x11x - x10x) * yoffseth) - x0yx) * zoffseth;
                }
            case 0x06: {
                    const FLOATTYPE x0y0 = src[0] + (src[param->stepy]-src[0])*yoffseth;
                    src += param->stepz;
                    return x0y0 + (src[0] + (src[param->stepy]-src[0])*yoffseth - x0y0)*zoffseth;
                }
            case 0x05: {
                    const FLOATTYPE x00x = src[0] + (src[1]-src[0])*xoffseth;
                    src += param->stepz;
                    return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)*zoffseth;
                }
            case 0x03: {
                    const FLOATTYPE x00x = src[0] + (src[1]-src[0])*xoffseth;
                    src += param->stepy;
                    return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)* yoffseth;
                }
            case 0x04:
                return src[0] + (src[param->stepz]-src[0])*zoffseth;
            case 0x02:
                return src[0] + (src[param->stepy]-src[0])*yoffseth;
            case 0x01:
                return src[0] + (src[1]           -src[0])*xoffseth;
            default:
                return src[0];
        }
    }
#endif

}






/****************************************************************************//**
 * \brief Nearest neighbour interpolation.
 *
 * \param[in] Pointer to the current interpolation position.
 * \return The interpolated value. If the coordinate is outside the volume,
 *   *((const FLOATTYPE *)param->optparam) is returned (i.e.
 *   *((const FLOATTYPE *)param->optparam) is the default value and must
 *   as a valid pointer to FLOATTYPE).
 *******************************************************************************/
FLOATTYPE tom_transf_interpol3d_nearest_double(const tom_transf_st_interpol_param *param) {
    const double x = param->x + 0.5;
    const double y = param->y + 0.5;
    const double z = param->z + 0.5;
    if (x<0 || y<0 || z<0 || x>=param->sizex || y>=param->sizey || z>=param->sizez) {
        return *((const FLOATTYPE *)param->optparam);
    }
    return ((const FLOATTYPE *)param->src)[((size_t)z)*param->stepz +
                                           ((size_t)y)*param->stepy +
                                           ((size_t)x)*1];
}








/****************************************************************************//**
 * \brief Applay a projective transformation to an volume.
 *
 * \param[in] src The source array
 * \param[in] srcsizex The size of the volume in x-direction.
 * \param[in] srcsizey The size of the volume in y-direction.
 * \param[in] srcsizez The size of the volume in z-direction.
 * \param[in] srcstepx How distant two neighbouring values of the source volume
 *   are separated in memory in x-direction. For example suppose the voxel
 *   src[z,y,x] is located at position &src[z,y,x] in memory. src[z,y,x+1]
 *   lies at position (&src[z,y,x] * srcstepx*sizeof(FLOATTYPE)). In most cases
 *   the memory will be continous. In that case srcstepx must be 1.
 * \param[in] srcstepy analog to srcstepx this is the distance along y-direction.
 *   The adresse of src[z,y+1,x] will be (&src[z,y,x] * srcstepy*sizeof(FLOATTYPE)).
 *   In case of continous memory, this value will be set to srcsizex.
 * \param[in] srcstepz analog to srcstepx this is the distance along z-direction.
 *   The adresse of src[z+1,y,x] will be (&src[z,y,x] * srcstepz*sizeof(FLOATTYPE)).
 *   In case of continous memory, this value will be set to srcsizex*srcsizey.
 * \param[out] dst The destination array. It must already be allocated according
 *   to the parameters dstsizex,dstsizey,dstsizez,dststepx,dststepy,dststepz.
 * \param[in] dstsizex Analog to srcsizex it is the resulting x-size of the output volume.
 * \param[in] dstsizey Analog to srcsizey it is the resulting y-size of the output volume.
 * \param[in] dstsizez Analog to srcsizez it is the resulting z-size of the output volume.
 * \param[in] dststepx Analog to srcstepx the step in memory from one x-value to the next.
 * \param[in] dststepy Analog to srcstepy the step in memory from one y-value to the next.
 * \param[in] dststepz Analog to srcstepz the step in memory from one z-value to the next.
 * \param[in] P The projection matrix. It must be a 4x4 matrix (non singlar).
 * \param[in] valinf If the P is a true projectivity (not an affine transformation, i.e.
 *   the last row does not equal [0.,0.,0.,x]) then some voxels in the destination volume
 *   can derive from an point at infinity. In that case the interpolation function is not
 *   called, but this value is assigned. In most cases this value is irrelevant.
 * \param[in] intfcn A pointer to the interpolation function. It returns the interplated value.
 * \param[in] optparam The optioal parameters of the interpolation function. This value is
 *   passed to the interpolation function through the global variable interpolparam_optparam.
 *   See the implementation of the used interplation function intfcn, what is expected here.
 *   For example it is usefull for the interpolation function to use this param, to get an
 *   default value if the value is outside the src-volume.
 * \return 0 if no error occured. Otherwise an error status is returned.
 *
 * This function fills the voxels of dst with interpolatates valued from src.
 * The output volume dst has (virtual) voxels at coordinates [xd,yd,zd] where
 * xd=0:dstsizex, yd=0:dstsizey, zd=0:dstsizez. These voxels are filled with
 * interpolated values of the source volume src at position
 * [w*xs, w*ys, w*zs, w]^T = P^(-1) * [xd,yd,zd,1]^T. The coordinates of the
 * given voxels of the source volume range from xs=0:(srcsizex-1), ys=0:(srcsizey-1),
 * zs=0:(srcsizez-1). Every element of the destination volume is set to valinf if
 * w==0 or to intfcn() otherwise. It depends entirely on the interpolation function
 * to calculate the destination voxel. Before calling intfcn the global variables
 * interpolparam_src, interpolparam_srcsizex, interpolparam_srcsizey,
 * interpolparam_srcsizez, 1, interpolparam_srcstepy,
 * interpolparam_srcstepz, interpolparam_px, interpolparam_py, interpolparam_pz,
 * interpolparam_optparam are set proper and can be used by the interpolation function
 * to determine the voxel-value. In this case interpolparam_px, interpolparam_py and
 * interpolparam_pz is set to xs, ys and zs respectively. The interpolation function
 * should not(!) modifiey these global variables :).
 * Note that by setting the projection matrix P to a specific function you can do
 * EVERY linear transformation (i.e. scaling, rotation, warping...). with
 * the step parameters, you can cut out a subvolume of a bigger source volume,
 * transform it, and return the result in an other subvolume. Note further
 * that the values outside the source volume are undefined and it depends
 * entirely on the interpolation function to handle that.
 *******************************************************************************/
int tom_transf_transf3d_double( const FLOATTYPE *src,
                                size_t srcsizex, size_t srcsizey, size_t srcsizez,
                                size_t srcstepy, size_t srcstepz,
                                FLOATTYPE *dst,
                                size_t dstsizex, size_t dstsizey, size_t dstsizez,
                                size_t dststepy, size_t dststepz,
                                const double *P,
                                int is_affinity,
                                FLOATTYPE valinf,
                                FLOATTYPE (*intfcn)(const tom_transf_st_interpol_param *param),
                                void *optparam) {

    double Pinv[16];
    size_t dstx, dsty, dstz;





    /* Check input parameters. */
    if (!src || !dst || !P || !intfcn ||
        srcsizex<1 || srcsizey<1 || srcsizez<1 ||
        srcstepy<srcsizex ||
        srcstepz<srcsizey*srcstepy ||
        dstsizex<1 || dstsizey<1 || dstsizez<1 ||
        dststepy<dstsizex ||
        dststepz<dstsizey*dststepy) {
        return -1;
    }

    /* Compute the inverse of the transformation matrix. */
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
            for (dstx=0; dstx<12; dstx++) {
                Pinv[dstx] /= Pinv[15];
            }
        }
    }
    if (!tom_transf_invert(Pinv, Pinv, 4)) { return -2; }
    if (is_affinity) {
        Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
    }

    dststepz -= dstsizey*dststepy;

    if (is_affinity &&
        Pinv[ 0]==1. && Pinv[ 1]==0. && Pinv[ 2]==0. && Pinv[ 3]==0. &&
        Pinv[ 4]==0. && Pinv[ 5]==1. && Pinv[ 6]==0. && Pinv[ 7]==0. &&
        Pinv[ 8]==0. && Pinv[ 9]==0. && Pinv[10]==1. && Pinv[11]==0. &&
        (intfcn==&tom_transf_interpol3d_trilinear_double ||
         intfcn==&tom_transf_interpol3d_nearest_double)) {
        /* Projection is identity. Do only copy
           Only in case of trilinear or nearest neighbour interpolation.
           Otherwise it may be desired to do the "interpolation" nonetheless,
           for example, if the interpolation-function is a low pass filter over
           several neighbours in the volume space. */
        size_t row_length = sizeof(FLOATTYPE)*dstsizex;

        srcstepz -= srcsizey*srcstepy;
        for (dstz=0; dstz<dstsizez; dstz++) {
            for (dsty=0; dsty<dstsizey; dsty++) {
                memcpy(dst, src, row_length);
                dst += dststepy;
                src += srcstepy;
            }
            dst += dststepz;
            src += srcstepz;
        }
    } else {

        /* Copy the projective transformation into single variables. */
        const double p00 = Pinv[ 0];
        const double p01 = Pinv[ 1];
        const double p02 = Pinv[ 2];
        const double p03 = Pinv[ 3];
        const double p10 = Pinv[ 4];
        const double p11 = Pinv[ 5];
        const double p12 = Pinv[ 6];
        const double p13 = Pinv[ 7];
        const double p20 = Pinv[ 8];
        const double p21 = Pinv[ 9];
        const double p22 = Pinv[10];
        const double p23 = Pinv[11];
        const double p30 = Pinv[12];
        const double p31 = Pinv[13];
        const double p32 = Pinv[14];
        const double p33 = Pinv[15];
        double tmpz_00, tmpz_01, tmpz_02;
        double tmpy_00, tmpy_01, tmpy_02;

        tom_transf_st_interpol_param intpar;


        /* Set parameters of the interpolation function. */
        intpar.src = src;
        intpar.sizex = srcsizex;
        intpar.sizey = srcsizey;
        intpar.sizez = srcsizez;
        intpar.stepy = srcstepy;
        intpar.stepz = srcstepz;
        intpar.x = 0.;
        intpar.y = 0.;
        intpar.z = 0.;
        intpar.optparam = optparam;


        dststepy -= dstsizex;
        if (is_affinity) {
            for (dstz=0; dstz<dstsizez; dstz++) {
                tmpz_00 = p02*dstz + p03;
                tmpz_01 = p12*dstz + p13;
                tmpz_02 = p22*dstz + p23;
                for (dsty=0; dsty<dstsizey; dsty++) {
                    tmpy_00 = p01*dsty + tmpz_00;
                    tmpy_01 = p11*dsty + tmpz_01;
                    tmpy_02 = p21*dsty + tmpz_02;
                    for (dstx=0; dstx<dstsizex; dstx++) {
                        intpar.x = p00*dstx + tmpy_00;
                        intpar.y = p10*dstx + tmpy_01;
                        intpar.z = p20*dstx + tmpy_02;

                        *dst++ = intfcn(&intpar);
                    }
                    dst += dststepy;
                }
                dst += dststepz;
            }
        } else {
            double pw, tmpz_03, tmpy_03;
            for (dstz=0; dstz<dstsizez; dstz++) {
                tmpz_00 = p02*dstz + p03;
                tmpz_01 = p12*dstz + p13;
                tmpz_02 = p22*dstz + p23;
                tmpz_03 = p32*dstz + p33;
                for (dsty=0; dsty<dstsizey; dsty++) {
                    tmpy_00 = p01*dsty + tmpz_00;
                    tmpy_01 = p11*dsty + tmpz_01;
                    tmpy_02 = p21*dsty + tmpz_02;
                    tmpy_03 = p31*dsty + tmpz_03;
                    for (dstx=0; dstx<dstsizex; dstx++) {
                        if((pw = p30*dstx + tmpy_03) == 0.) {
                            *dst++ = valinf;
                        } else {
                            intpar.x = (p00*dstx + tmpy_00) / pw;
                            intpar.y = (p10*dstx + tmpy_01) / pw;
                            intpar.z = (p20*dstx + tmpy_02) / pw;
                            *dst++ = intfcn(&intpar);
                        }
                    }
                    dst += dststepy;
                }
                dst += dststepz;
            }
        }
    }
    return 0;
}







/****************************************************************************//**
 * \brief Projects a 3D volume to a plane where the volume is rotated.
 *
 * Projects a 3D volume onto a plane. The volume lies continous in memory and
 * is rotated with angles phi and theta around its center. Then it is projected
 * onto the image which lies in the XY-plane. This is the same as tom_proj3d and
 * tom_proj3d2c from the matlab toolbox TOM. There is also an mex-interface for
 * it in tom_mex_proj3d.c.
 * \param[in] src The source volume.
 * \param[in] srcsizex The size of the volume.
 * \param[in] srcsizey The size of the volume.
 * \param[in] srcsizez The size of the volume.
 * \param[in] phi Rotate with angle phi. The angle is given in degrees.
 * \param[in] theta Rotate with angle theta. The angle is given in degrees.
 * \param[in] dst The destination image. It has (implicitly) the size
 *   srcsizex*srcsizey.
 * \param[in] normalize If true, each pixel of the resulting image is diveded by
 *   its weight. Thus it contains the mean, rather then the sum of all voxels which
 *   are projected on it.
 * \param[in] defval This parameter is the default value to set on a pixel
 *   if no voxel projected on it. It is ignored if !normalize.
 * \return Passes the result from tom_transf_proj3d_double. That means 0 in case of
 *   invalid arguments and 1 otherwise.
 *
 * The rotation (around the center of the volume) is done in the following
 * order: first with (90-phi) degrees around the Z-axis, then with -theta
 * degrees around the X-axis and finally again with -(90-phi) degrees around
 * Z. This function is only a wrapper for \link tom_transf_proj3d_double \endlink.
 *******************************************************************************/
int tom_transf_proj3d_solid_zxz_double(const FLOATTYPE *src, size_t srcsizex, size_t srcsizey, size_t srcsizez, double phi, double theta, FLOATTYPE *dst, int normalize, FLOATTYPE defval) {

    const double center_x = ((double)srcsizex)/2.;
    const double center_y = ((double)srcsizey)/2.;
    const double center_z = ((double)srcsizez)/2.;
    double P[12];
    const double deg2rad = __PI__ / 180.;
    const int axes[3] = { 2, 0, 2 };
    double angles[3];
    double *weight = NULL;
    int res;

    /* Initialize the projection matrix. */
    angles[0] =  (90. - phi) * deg2rad;
    angles[1] = -    theta  * deg2rad;
    angles[2] = -angles[0];
    tom_transf_init_rotmat_3d(P, 4, 3, axes, angles);

    P[ 3] = center_x - (P[ 0]*center_x + P[ 1]*center_y + P[ 2]*center_z);
    P[ 7] = center_y - (P[ 4]*center_x + P[ 5]*center_y + P[ 6]*center_z);
    P[11] = center_z - (P[ 8]*center_x + P[ 9]*center_y + P[10]*center_z);

    if (normalize) {
        weight = (double *)malloc(srcsizey*srcsizex*sizeof(double));
        if (!weight) {
            return 0;
        }
    }

    res = tom_transf_proj3d_double(src, srcsizex, srcsizey, srcsizez, srcsizex, srcsizex*srcsizey, P, dst, srcsizex, srcsizey, srcsizex, weight);

    if (normalize) {
        size_t ix, iy;
        FLOATTYPE *pdst = dst;
        double *pweight = weight;
        double w;
        for (iy=0; iy<srcsizey; iy++) {
            for (ix=0; ix<srcsizex; ix++) {
                w = *pweight++;
                if (w) {
                    *pdst++ /= w;
                } else {
                    *pdst++ = defval;
                }
            }
        }
        free(weight);
    }

    return res;
}



/****************************************************************************//**
 * \brief projects a 3d-volume onto a plane.
 *
 * Takes an input volume and projects it onto a 2D-plane. This is done by
 * walking throug every voxel of the source volume, calculate its projection
 * and distribute it to the 4 closest pixels in the destination image
 * (somewhat of the inverse operation of bilinear interpolation).
 *
 * \param[in] src A pointer to the source volume.
 * \param[in] srcsizex The size of the input volume along X.
 * \param[in] srcsizey The size of the input volume along Y.
 * \param[in] srcsizez The size of the input volume along Z.
 * \param[in] srcstepy The number of elements in one row.
 *   Usually, in case of continous memory, this equals the size
 *   of the volume along X (i.e. srcsizex).
 * \param[in] srcstepz The number of elements in one "level" of the
 *   volume. In case of continous memory it is srcsizey*srcsizex.
 * \param[in] P A pointer to the projection matrix. It must contain
 *   8 elements and is the upper 2x4 part of an 3x4 affine projection
 *   matrix (with last row equal [0,0,0,1]).
 * \param[out] dst A pointer to the image. It contains the sum of all
 *   projected voxels.
 * \param[in] dstsizex The size of the output image.
 * \param[in] dstsizey The size of the output image.
 * \param[in] dststepy The number of elements in one row of the image.
 *   The memory pointed by dst must be at least dststepy*dstsizey*sizeof(*dst)
 *   bytes large.
 * \param[out] weight An optional output parameter with the sum of how
 *   many voxels contributed to each pixel. If not NULL, it must be
 *   dstsizex*dstsizey*sizeof(*weight) bytes large with no gap between
 *   the rows (stepy=dstsizex). This parameter may be usefull to
 *   to devide each resulting pixel at the end to get its mean.
 * \return In case of invalid input parameters it returns 0, otherwise 1.
 *
 * Setting the projection matrix P to the right values, every affine
 * projection (i.e. parallel projection) can be done. A voxel has
 * 3D coordinates [xv,yv,zv] where xv ranges from 0 to srcsizex-1 (yv,zv analog).
 * A 2D destination pixel has image coodinates [xi,yi] where xi ranges from
 * 0 to dstsizex-1 (yi analog). The 3D point is projected to the image point
 * P * [xv,yv,zv,1]^T, where P is a 2x4 matrix. See proj3d_solid_zxz_GENERIC
 * how the matrix is specialized to a rotation around the center and projection
 * along Z.
 * \todo Is there no weighting necessary? Rethink about it :)
 *
 * \par Credits
 * This is an reimplementation of an old fortran source:
\verbatim
SUBROUTINE PROJ (A,JX,NA,EM,FR)
C-----------------------------------------------------------------------
C     BERECHNET DIE PROJEKTIONEN EINES ZWEI- ODER DREIDIMENSIONALEN
C     ARRAYS IN BELIEBIGER RICHTUNG.
C     APRIL 16, 1980, HEGERL
C     Last update:  16-Apr-1991,  R. Hegerl
C     29-JUL-1992: Option "U" permits the projection of stacks with
C                  arbitrary height. In this case
C                  FR(3) = total height of the input array (=NZ)
C                  FR(4) = running number of slice to be projected
C                             (R. Hegerl)
C     23-SEP-1992: The output array may have arbitrary dimensions,
C                  defined in FR(3),FR(4)  (R. Hegerl)
C-----------------------------------------------------------------------
      INCLUDE '../common/emsys.f'
      INCLUDE '../common/emfile.f'
      INCLUDE '../common/emary.f'
C
      REAL      A(1), EM(40,1), FR(1)
      INTEGER*4 JX(1), NA(16,1)
C
C      REAL      SR             ! REAL FUNCTION
      INTEGER*4 JCHECK, MA     ! INTEGER FUNCTION
C
C  Local variables
      REAL      CPHI, CTHE, DIPX, DIPY, DP(3), DPNEW(2), DM(3,3),
     $          DX, DXP, DY, DYP, DZ, H, PHIR, RLX, RLXNEW,
     $          RLY, RLYNEW, SPHI, STHE, THETAR, W,
     $          X(3), XP(2), XFACT, YFACT, ZWX, ZWY
      INTEGER*4 I, I1, I2, I3, I4, IC(3), INDX, INDP, IPX, IPY,
     $          JX1, JX2, KA, NI, NX, NXP, NXYZ, NY, NYP, NZ
      LOGICAL   PRU, LIST
C
C CHECK THE INPUT
      PRU = JCHECK('U').GT.0
      LIST = JCHECK('!').GT.0
      JX1 = JX(1)
      JX2 = JX(2)
      IF (NA(1,JX1).NE.4) CALL ERREND(63,*800,'PROJ * ARRAY NOT REAL')
      NX = NA(2,JX1)
      NY = NA(3,JX1)
      NZ = NA(4,JX1)
      NXYZ = NX*NY*NZ
      IF (PRU) THEN
        NZ = FR(3)
        NXYZ = NX*NY
      END IF
      IF (.NOT.PRU .AND. JCHECK('V').EQ.0) THEN
        CALL PURGE(JX2)
        CALL TRANS(JX1,JX2)
        IF (FR(3).NE.0) NA(2,JX2) = FR(3)
        IF (FR(4).NE.0) NA(3,JX2) = FR(4)
      END IF
      NA(4,JX2) = 1
      IF (NZ.EQ.1) NA(3,JX2) = 1
      KA = MA(NA(1,JX2))
      NA(6,JX2) = KA
      NXP = NA(2,JX2)
      NYP = NA(3,JX2)
C
C  Check memory setup
      IC(1) = NXYZ
      IF (NA(9,JX1).NE.0) IC(1) = IBLK
      IC(2) = KA
      CALL CHADR (JX,2,NA,IC,'PROJ',*800)
C
      IF (EM(1,JX1).EQ.0) THEN
C       WRITE (NTOUT,907)
C907    FORMAT (/' WARNING: PROJ * EM-DATA MISSING')
        RLX = NX
        H = 0.
        IF (NZ.GT.1) H = NZ
      ELSE
        RLX = EM(10,JX1)*1.E4/(EM(4,JX1)*EM(5,JX1))
        H = EM(33,JX1)*1.E4/(EM(4,JX1)*EM(5,JX1))
      END IF
      RLY = RLX*NY/NX
      DX = RLX/NX
      DY = RLY/NY
      DZ = H/NZ
      NI = 2
      IF (DZ.GT.0.) NI = 3
      IF (NZ.GT.1.AND.NI.EQ.2.OR.NZ.EQ.1.AND.NI.EQ.3) CALL
     $   ERREND(66,*800,'PROJ * WRONG EM-DATA: H MISSING OR WRONG')
C
      DP(1) = (NX/2)*DX
      DP(2) = (NY/2)*DY
      DP(3) = (NZ/2)*DZ
      XFACT = 1.0
      YFACT = 1.0
      IF (FR(5).GT.0.0) XFACT = FR(5)
      IF (FR(6).GT.0.0) YFACT = FR(6)
      DXP = DX*XFACT
      DYP = DY*YFACT
      DPNEW(1) = (NXP/2)*DXP
      DPNEW(2) = (NYP/2)*DYP
      RLXNEW = RLX*FLOAT(NXP)/FLOAT(NX)
      RLYNEW = RLXNEW*NYP/NXP
C
C BERECHNUNG DER PROJEKTIONSRICHTUNG UND -EBENE
      PHIR = FR(1)*PI/180
      THETAR = FR(2)*PI/180
      CPHI=COS(PHIR)
      SPHI=SIN(PHIR)
      CTHE=COS(THETAR)
      STHE=SIN(THETAR)
      IF (NI.EQ.2) THEN
        DM(1,1)=CPHI
        DM(2,1)=-SPHI
        DM(1,2)=SPHI
        DM(2,2)=CPHI
      ELSE
        DM(1,1) = CTHE*CPHI*CPHI+SPHI*SPHI
        DM(2,1) = CTHE*CPHI*SPHI-CPHI*SPHI
        DM(3,1) = -STHE*CPHI
        DM(1,2) = DM(2,1)
        DM(2,2) = CTHE*SPHI*SPHI+CPHI*CPHI
        DM(3,2) = -STHE*SPHI
        DM(1,3) = -DM(3,1)
        DM(2,3) = -DM(3,2)
        DM(3,3) = CTHE
      END IF
      IF (LIST .AND. (.NOT.PRU .OR. (PRU.AND.FR(4).EQ.1))) THEN
        CALL LINE_COUNT(1)
        WRITE (NTOUT,901) (DM(I,NI), I=1,NI)
 901    FORMAT ('   PROJECTION DIRECTION:',3F8.3)
      END IF
C
C INITIALISIERUNG
      IF (.NOT.PRU) THEN
        DO 40   I1 = 1,KA
  40      A(I1) = 0.
      END IF
C
C BERECHNUNG DES PROJEKTIONSPUNKTES
      X(3) = -DZ
C     INDX = 1
      INDX = NA(15,JX1) + 1
      IF (PRU) THEN
        NZ = 1
        X(3) = -DZ + (FR(4)-1)*DZ
      END IF
      DO 115   I1=1,NZ
      X(3)=X(3)+DZ
      X(2)=-DY
      DO 110   I2=1,NY
        X(2) = X(2) + DY
        X(1) = -DX
        DO 105   I3=1,NX
          X(1) = X(1) + DX
          XP(1) = 0.0
          XP(2) = 0.0
          DO 100   I4 = 1,NI
            XP(1)=XP(1)+(X(I4)-DP(I4))*DM(I4,1)
 100        XP(2)=XP(2)+(X(I4)-DP(I4))*DM(I4,2)
          XP(1) = XP(1) + DPNEW(1)
          XP(2) = XP(2) + DPNEW(2)
          IF (NI.EQ.2) XP(2) = 0.
          IF  (XP(1).LT.0.0.OR.XP(1).GT.RLXNEW-DXP.OR.XP(2).LT.0.0
     $         .OR.XP(2).GT.RLYNEW-DYP)  GO TO 105
C
C INTERPOLATION
          ZWX=XP(1)/DXP
          ZWY=XP(2)/DYP
          IPX=IFIX(ZWX)
          IPY=IFIX(ZWY)
          DIPX=ZWX-FLOAT(IPX)
          DIPY=ZWY-FLOAT(IPY)
C         W = SR(INDX,JX1)
          W = ARRAY(INDX)
          INDP = IPX+1+IPY*NXP
          A(INDP) = A(INDP) + (1.0-DIPX)*(1.0-DIPY)*W
          A(INDP+1) = A(INDP+1) + DIPX*(1.0-DIPY)*W
          IF (NI.EQ.3) THEN
            INDP = INDP + NXP
            A(INDP) = A(INDP) + (1.0-DIPX)*DIPY*W
            A(INDP+1) = A(INDP+1) + DIPX*DIPY*W
          END IF
 105      INDX = INDX + 1
 110    CONTINUE
 115  CONTINUE
      IF (NI.EQ.3) EM(33,JX2) = 0.
 800  RETURN
      END
\endverbatim
 *******************************************************************************/
int tom_transf_proj3d_double(const FLOATTYPE *src, size_t srcsizex, size_t srcsizey, size_t srcsizez, size_t srcstepy, size_t srcstepz,
      const double *P, FLOATTYPE *dst, size_t dstsizex, size_t dstsizey, size_t dststepy, double *weight) {


    size_t ix, iy , iz;

    if (!src || !dst ||
        srcsizex<1 || srcsizey<1 || srcsizez<1 ||
        srcstepy<srcsizex || srcstepz<srcstepy*srcsizey ||
        dstsizex<1 || dstsizey<1 ||
        dststepy<dstsizex) {
        return 0;
    }

    /*  initialization output projection */
    {
        FLOATTYPE *pdst = dst;
        for (iy=0; iy<dstsizey; iy++) {
            for (ix=0; ix<dstsizex; ix++) {
                pdst[ix] = 0.;
            }
            pdst += dststepy;
        }
    }


    {
        const double p00 = P[ 0];
        const double p01 = P[ 1];
        const double p02 = P[ 2];
        const double p03 = P[ 3];
        const double p10 = P[ 4];
        const double p11 = P[ 5];
        const double p12 = P[ 6];
        const double p13 = P[ 7];
        const size_t dstsizex_m_1 = dstsizex-1;
        const size_t dstsizey_m_1 = dstsizey-1;
        double p0_zw, p1_zw, p0_yzw, p1_yzw;
        FLOATTYPE value;
        struct {
            double x;
            double y;
        } p, dist_p, dist_m_1_p;
        struct {
            size_t x;
            size_t y;
        } pi;
        size_t idx_dst;
        srcstepz -= srcstepy*srcsizey;

        if (weight) {
            size_t idx_weight;
            {
                double *pweight = weight;
                for (iy=0; iy<dstsizey; iy++) {
                    for (ix=0; ix<dstsizex; ix++) {
                        *pweight++ = 0.;
                    }
                }
            }
            for (iz=0; iz<srcsizez; iz++) {
                p0_zw = p02*iz + p03;
                p1_zw = p12*iz + p13;
                for (iy=0; iy<srcsizey; iy++) {
                    p0_yzw = p01*iy + p0_zw;
                    p1_yzw = p11*iy + p1_zw;
                    for (ix=0; ix<srcsizex; ix++) {
                        p.x = p00*ix + p0_yzw;
                        if (p.x<0 || p.x>dstsizex_m_1) {
                            continue;
                        }
                        p.y = p10*ix + p1_yzw;
                        if (p.y<0 || p.y>dstsizey_m_1) {
                            continue;
                        }
                        pi.x = (size_t)floor(p.x);
                        pi.y = (size_t)floor(p.y);
                        dist_p.x = p.x - (double)pi.x;
                        dist_p.y = p.y - (double)pi.y;
                        dist_m_1_p.x = 1.0 - dist_p.x;
                        dist_m_1_p.y = 1.0 - dist_p.y;
                        idx_dst = pi.y*dststepy + pi.x;
                        idx_weight = pi.y * dstsizex + pi.x;
                        value = src[ix];
                        dst   [idx_dst   ] += (FLOATTYPE)(dist_m_1_p.x * dist_m_1_p.y * value);
                        weight[idx_weight] +=             dist_m_1_p.x * dist_m_1_p.y;
                        if (dist_p.x) {
                            dst[idx_dst      +1] += (FLOATTYPE)(dist_p.x * dist_m_1_p.y * value);
                            weight[idx_weight+1] +=             dist_p.x * dist_m_1_p.y;
                            if (dist_p.y) {
                                dst   [idx_dst   +=dststepy] += (FLOATTYPE)(dist_m_1_p.x * dist_p.y * value);
                                weight[idx_weight+=dstsizex] +=             dist_m_1_p.x * dist_p.y;
                                dst   [idx_dst   + 1       ] += (FLOATTYPE)(dist_p.x     * dist_p.y * value);
                                weight[idx_weight+ 1       ] +=             dist_p.x     * dist_p.y;
                            }
                        } else if (dist_p.y) {
                            dst   [idx_dst   + dststepy] += (FLOATTYPE)(dist_m_1_p.x * dist_p.y * value);
                            weight[idx_weight+ dstsizex] +=             dist_m_1_p.x * dist_p.y;
                        }
                    }
                    src += srcstepy;
                }
                src += srcstepz;
            }
        } else {
            for (iz=0; iz<srcsizez; iz++) {
                p0_zw = p02*iz + p03;
                p1_zw = p12*iz + p13;
                for (iy=0; iy<srcsizey; iy++) {
                    p0_yzw = p01*iy + p0_zw;
                    p1_yzw = p11*iy + p1_zw;
                    for (ix=0; ix<srcsizex; ix++) {
                        p.x = p00*ix + p0_yzw;
                        if (p.x<0 || p.x>dstsizex_m_1) {
                            continue;
                        }
                        p.y = p10*ix + p1_yzw;
                        if (p.y<0 || p.y>dstsizey_m_1) {
                            continue;
                        }
                        pi.x = (size_t)floor(p.x);
                        pi.y = (size_t)floor(p.y);
                        dist_p.x = p.x - (double)pi.x;
                        dist_p.y = p.y - (double)pi.y;
                        dist_m_1_p.x = 1.0 - dist_p.x;
                        dist_m_1_p.y = 1.0 - dist_p.y;
                        idx_dst = pi.y*dststepy + pi.x;
                        value = src[ix];
                        dst[idx_dst] += (FLOATTYPE)(dist_m_1_p.x * dist_m_1_p.y * value);
                        if (dist_p.x) {
                            dst[idx_dst+1] += (FLOATTYPE)(dist_p.x * dist_m_1_p.y * value);
                            if (dist_p.y) {
                                dst[idx_dst+=dststepy] += (FLOATTYPE)(dist_m_1_p.x * dist_p.y * value);
                                dst[idx_dst+ 1       ] += (FLOATTYPE)(dist_p.x     * dist_p.y * value);
                            }
                        } else if (dist_p.y) {
                            dst[idx_dst+dststepy] += (FLOATTYPE)(dist_m_1_p.x * dist_p.y * value);
                        }
                    }
                    src += srcstepy;
                }
                src += srcstepz;
            }
        }
    }

    return 1;

}









