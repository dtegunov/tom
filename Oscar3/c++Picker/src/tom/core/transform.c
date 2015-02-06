/****************************************************************************//**
 * \file transform.c
 * \author Thomas Haller
 * \date Oct. 26, 2007
 * \brief Implementation of the methods in transform.h
 *******************************************************************************/
#include "tom/core/transform.h"


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>



#define __PI__              3.14159265358979323846264338327950288419716939937510


#if defined _WIN32
    #define dgetri_ dgetri
    #define dgetrf_ dgetrf
#endif


#ifdef __cplusplus
    #define __EXTERN_LAPACK extern "C"
#else
    #define __EXTERN_LAPACK extern
#endif
/****************************************************************************//**
 * \brief Prototype for LAPACK functions needed for matrix inversion.
 *
 * See http://www.netlib.org/lapack/double/dgetri.f
 *******************************************************************************/
__EXTERN_LAPACK int dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO );
/** \copydoc dgetri_ */
__EXTERN_LAPACK int dgetrf_(int *N, int *M, double *A, int *LDA, int *IPIV, int *INFO);




/****************************************************************************//**
 * \brief Constructs a 3D rotation matrix of euler angles.
 *
 * \param[out] P Memory to the rotation matrix. The necessary size depends on
 *   the parameter homogen.
 * \param[in] stepy The number of matrix elements in P, how distante two
 *   rows of P are in memory. Obviously it must be greater or equal 3.
 *   For example the first element in the second row is at position &P[stepy].
 *   Usually the elements lie continous in memory, one element next the other,
 *   then stepy==3.
 * \param[in] nrots The number of rotations to do.
 * \param[in] axes An array of length nrots which specifies around
 *   which axis the i-th rotation has to be done. It means,
 *   first it is rotated around axes[0], then axes[1] and so on.
 *   The x axes corresponds to 0, y to 1 and z to 2.
 *   The orientation of the rotation is the following:
 *     \li Around X-axes (axes[i]==0):
 *         \verbatim
    [[              1,              0,              0];
     [              0, cos(angles[i]),-sin(angles[i])];
     [              0, sin(angles[i]), cos(angles[i])]]
\endverbatim
 *     \li Around Y-axes (axes[i]==1):
 *         \verbatim
    [[ cos(angles[i]),              0, sin(angles[i])];
     [              0,              1,              0];
     [-sin(angles[i]),              0, cos(angles[i])]]
\endverbatim
 *     \li Around Z-axes (axes[i]==2):
 *         \verbatim
    [[ cos(angles[i]),-sin(angles[i]),              0];
     [ sin(angles[i]), cos(angles[i]),              0];
     [              0,              0,              0]]
\endverbatim
 * \param[in] angles The rotation angles in radians!!! of each
 *   of the nrots rotations.
 * \return Returns 1 in case of success, otherwise 0.
 *   Failure can only occure, when input-parameters are wrong. In that case
 *   P stays unchanged.
 *
 * A 3D rotation has a 3x3 matrix as result. Only the 9 elements of P are touched
 * and its entrence-value is irrelevant.
 *
 * \par
 * For example to get a rotation matrix according to the used standard in TOM,
 * rotate around the axes (zxz) with angles (phi, theta, psi) (in radians).
 * The result will be the same as for example with tom_sum_rotation.m
 * \code
 * double P[9];
 * int axes[3] = { 2, 0, 2 };
 * double angles[3] = { phi, theta, psi }
 * init_rotmat_3d(P, 3, 3, axes, angles);
 * \endcode
 *******************************************************************************/
int tom_transf_init_rotmat_3d(double *P, size_t stepy, size_t nrots, const int *axes, const double *angles) {

    int i, j;
    size_t ni;
    double sina, cosa;
    double p0, p1;

    if (!P || (nrots && (!axes||!angles)) || stepy<3) {
        return 0;
    }
    for (ni=0; ni<nrots; ni++) {
        if (axes[ni] < 0 || axes[ni]>=3) {
            return 0;
        }
    }

    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            if (j==i) {
                P[i*stepy + j] = 1.;
            } else {
                P[i*stepy + j] = 0.;
            }
        }
    }
    for (ni=0; ni<nrots; ni++) {
        if (angles[ni] - floor(angles[ni]/(__PI__/2.))*(2.*__PI__) != 0.) {
            sina = sin(angles[ni]);
            cosa = cos(angles[ni]);
            switch (axes[ni]) {
                case 0:
                    for (i=0; i<3; i++) {
                        p0 = P[i*stepy+1];
                        p1 = P[i*stepy+2];
                        P[i*stepy+1] = + p0*cosa - p1*sina;
                        P[i*stepy+2] = + p0*sina + p1*cosa;
                    }
                    break;
                case 1:
                    for (i=0; i<3; i++) {
                        p0 = P[i*stepy+0];
                        p1 = P[i*stepy+2];
                        P[i*stepy+0] = + p0*cosa + p1*sina;
                        P[i*stepy+2] = - p0*sina + p1*cosa;
                    }
                    break;
                case 2:
                    for (i=0; i<3; i++) {
                        p0 = P[i*stepy+0];
                        p1 = P[i*stepy+1];
                        P[i*stepy+0] = + p0*cosa - p1*sina;
                        P[i*stepy+1] = + p0*sina + p1*cosa;
                    }
                    break;
            }
        }
    }
    return 1;
}




/****************************************************************************//**
 * \brief Compute inverse of a matrix.
 *
 * \param[in]  P the input matrix.
 * \param[out] Pinv the output matrix inverted.
 * \param[in]  N The size of the matrix.
 * \returns 0 if the input arguments are non valid or the matrix is
 *   singular.
 *
 * P and Pinv are matrices of size N and lie continous in memory
 * (for example the third element in the second row is
 * P[1,2] = P[1*N + 2]).
 * Is self assignment safe.
 * If the function returns 0 the content of Pinv is undefined.
 *******************************************************************************/
int tom_transf_invert(const double *P, double *Pinv, int N) {

#if 0
    int INFO;
    int LWORK = N*N;
    int *IPIV = NULL;
    double *WORK = NULL;

    if (!(IPIV = (int *)malloc(sizeof(int)*N)) ||
        !(WORK = (double *)malloc(sizeof(double)*LWORK))) {
        if (IPIV) { free(IPIV); }
        return 0;
    }

    if (P != Pinv) {
        size_t i;
        size_t NN = N*N;
        for (i=0; i<NN; i++) {
            Pinv[i] = P[i];
        }
    }

    dgetrf_(&N, &N, Pinv, &N, IPIV, &INFO);

    if (INFO) {
        free(IPIV);
        free(WORK);
        return 0;
    }

    dgetri_(&N, Pinv, &N, IPIV, WORK, &LWORK, &INFO);

    free(IPIV);
    free(WORK);
#else
    const double *Plocal = P;
    int i;
    double div;
    double *Pl = NULL;
    if (N==1) {
        Pinv[0] = 1./P[0];
    } else if (N==2) {
        if (P == Pinv) {
            Pl = (double *)malloc(sizeof(double) * N*N);
            if (!Pl) {
                return 0;
            }
            for (i=0; i<N*N; i++) {
                Pl[i] = P[i];
            }
            Plocal = Pl;
        }
        div = Plocal[0]*Plocal[3] - Plocal[1]*Plocal[2];
        if (!div) {
            return 0;
        }
        Pinv[0] =   Plocal[3] / div;
        Pinv[1] = - Plocal[1] / div;
        Pinv[2] = - Plocal[2] / div;
        Pinv[3] =   Plocal[0] / div;
        if (Pl) { free(Pl); }
    } else if (N==3) {
        if (P == Pinv) {
            Pl = (double *)malloc(sizeof(double) * N*N);
            if (!Pl) {
                return 0;
            }
            for (i=0; i<N*N; i++) {
                Pl[i] = P[i];
            }
            Plocal = Pl;
        }
        #if 0
        | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
        | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
        | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

        with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)
        #endif
        div =        + Plocal[0]*(Plocal[8]*Plocal[4] - Plocal[7]*Plocal[5])
                     - Plocal[3]*(Plocal[8]*Plocal[1] - Plocal[7]*Plocal[2])
                     + Plocal[6]*(Plocal[5]*Plocal[1] - Plocal[4]*Plocal[2]);
        if (!div) {
            return 0;
        }
        Pinv[0] =   (Plocal[4]*Plocal[8] - Plocal[5]*Plocal[7]) / div;
        Pinv[1] =   (Plocal[7]*Plocal[2] - Plocal[8]*Plocal[1]) / div;
        Pinv[2] =   (Plocal[1]*Plocal[5] - Plocal[2]*Plocal[4]) / div;
        Pinv[3] =   (Plocal[5]*Plocal[6] - Plocal[3]*Plocal[8]) / div;
        Pinv[4] =   (Plocal[0]*Plocal[8] - Plocal[2]*Plocal[6]) / div;
        Pinv[5] =   (Plocal[2]*Plocal[3] - Plocal[0]*Plocal[5]) / div;
        Pinv[6] =   (Plocal[3]*Plocal[7] - Plocal[4]*Plocal[6]) / div;
        Pinv[7] =   (Plocal[1]*Plocal[6] - Plocal[0]*Plocal[7]) / div;
        Pinv[8] =   (Plocal[0]*Plocal[4] - Plocal[1]*Plocal[3]) / div;
        if (Pl) { free(Pl); }
    } else if (N==4) {
        double Ainv[9];
        double D_mCAinvB__inv;
        double CAinv[3];
        double AinvB[3];

        Ainv[0] = P[ 0];
        Ainv[1] = P[ 1];
        Ainv[2] = P[ 2];
        Ainv[3] = P[ 4];
        Ainv[4] = P[ 5];
        Ainv[5] = P[ 6];
        Ainv[6] = P[ 8];
        Ainv[7] = P[ 9];
        Ainv[8] = P[10];
        if (!tom_transf_invert(Ainv, Ainv, 3)) {
            return 0;
        }
        CAinv[0] = P[12]*Ainv[0] + P[13]*Ainv[3] + P[14]*Ainv[6];
        CAinv[1] = P[12]*Ainv[1] + P[13]*Ainv[4] + P[14]*Ainv[7];
        CAinv[2] = P[12]*Ainv[2] + P[13]*Ainv[5] + P[14]*Ainv[8];
        AinvB[0] = Ainv[0]*P[ 3] + Ainv[1]*P[ 7] + Ainv[2]*P[11];
        AinvB[1] = Ainv[3]*P[ 3] + Ainv[4]*P[ 7] + Ainv[5]*P[11];
        AinvB[2] = Ainv[6]*P[ 3] + Ainv[7]*P[ 7] + Ainv[8]*P[11];
        D_mCAinvB__inv = P[15] - ( CAinv[0]*P[ 3] + CAinv[1]*P[ 7] + CAinv[2]*P[11] );
        if (!D_mCAinvB__inv) {
            return 0;
        }
        D_mCAinvB__inv = 1. / D_mCAinvB__inv;

        Pinv[ 0] = Ainv[0] + AinvB[0]*D_mCAinvB__inv*CAinv[0];
        Pinv[ 1] = Ainv[1] + AinvB[0]*D_mCAinvB__inv*CAinv[1];
        Pinv[ 2] = Ainv[2] + AinvB[0]*D_mCAinvB__inv*CAinv[2];
        Pinv[ 3] = - AinvB[0] * D_mCAinvB__inv;
        Pinv[ 4] = Ainv[3] + AinvB[1]*D_mCAinvB__inv*CAinv[0];
        Pinv[ 5] = Ainv[4] + AinvB[1]*D_mCAinvB__inv*CAinv[1];
        Pinv[ 6] = Ainv[5] + AinvB[1]*D_mCAinvB__inv*CAinv[2];
        Pinv[ 7] = - AinvB[1] * D_mCAinvB__inv;
        Pinv[ 8] = Ainv[6] + AinvB[2]*D_mCAinvB__inv*CAinv[0];
        Pinv[ 9] = Ainv[7] + AinvB[2]*D_mCAinvB__inv*CAinv[1];
        Pinv[10] = Ainv[8] + AinvB[2]*D_mCAinvB__inv*CAinv[2];
        Pinv[11] = - AinvB[2] * D_mCAinvB__inv;
        Pinv[12] = - D_mCAinvB__inv * CAinv[0];
        Pinv[13] = - D_mCAinvB__inv * CAinv[1];
        Pinv[14] = - D_mCAinvB__inv * CAinv[2];
        Pinv[15] = D_mCAinvB__inv;
    } else {
        printf("ERROR: NOT IMPLEMENTED MATRIX INVERSION OF ORDER>4 (%s %d)", __FILE__, __LINE__);
        return 0;
    }
#endif
    return 1;

}







/****************************************************************************//**
 * \brief Call the interpolation functions from tom_interpol.h
 *
 * Beware, the memory must be continous. i.e stepy==sizex and stepz==sizex*sizey
 * Set the optparam to a valid pointer of type interpol3d_continous_float_param.
 *******************************************************************************/
float tom_transf_interpol3d_continous_float(const tom_transf_st_interpol_param *param) {
    if (param->x < 0 || param->x > param->sizex-1 ||
        param->y < 0 || param->y > param->sizey-1 ||
        param->z < 0 || param->z > param->sizez-1) {
        return ((const tom_transf_interpol3d_continous_float_param *)param->optparam)->defaultval;
    }
    return interpolate3D((const float *)param->src, param->sizex, param->sizey, param->sizez, param->x, param->y, param->z, ((const tom_transf_interpol3d_continous_float_param *)param->optparam)->method);
}


#undef FLOATTYPE_ISFLOAT
#undef FLOATTYPE_ISDOUBLE
#define FLOATTYPE_ISFLOAT 1
#include "_transform.c"



#undef FLOATTYPE_ISFLOAT
#undef FLOATTYPE_ISDOUBLE
#define FLOATTYPE_ISDOUBLE 1
#include "_transform.c"










