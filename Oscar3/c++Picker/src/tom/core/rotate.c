/****************************************************************************//**
 * \file rotate.c
 * \author Thomas Hrabe
 * \date Jan. 24, 2007
 *******************************************************************************/
#include "tom/core/rotate.h"




/* rot3d.c    Performs 3D Rotation
* The syntax is:
*
*        rot3d(IN,OUT,PHI,PSI,THETA,INTERP)
*
*
*
* ---------------------------------------------
* ROTATION INFORMATION
* ---------------------------------------------
*
* This programm uses the following 3d Euler rotation
* (first index: line, second index: column):
*
* rm00=cospsi*cosphi-costheta*sinpsi*sinphi;
* rm10=sinpsi*cosphi+costheta*cospsi*sinphi;
* rm20=sintheta*sinphi;
* rm01=-cospsi*sinphi-costheta*sinpsi*cosphi;
* rm11=-sinpsi*sinphi+costheta*cospsi*cosphi;
* rm21=sintheta*cosphi;
* rm02=sintheta*sinpsi;
* rm12=-sintheta*cospsi;
* rm22=costheta;
*
* This is a 3d Euler rotation around (zxz)-axes
* with Euler angles (psi,theta,phi).
*
* For more information about 3d Euler rotation visit:
* wiki -> Definitions -> Overview 3d Euler rotation
 *
 *For detailed information about spline interpolation, see
 *http://mathworld.wolfram.com/CubicSpline.html
* ---------------------------------------------
*=================================================================*/


#include <stdio.h>
#include <time.h>
#include <stddef.h>
#include <math.h>


#include "tom/core/interpol.h"




#define  PI ((double)3.14159265358979323846264338327950288419716939937510)


/* 2D Rotation */



/*linear*/
void rot2dlin (float *image,float *rotimg,long sx,long sy,float *p_phi,float px,float py,int euler_dim)
{
float rm00, rm01, rm10, rm11;   /* rot matrix */
float pi, pj;           /* coordinates according to pivot */
float r_x, r_y;         /* rotated pixel */
long  i, j;         /* loop variables */
long  sx_1, sy_1;   /* highest pixels  | image width and height */
/*float *imgp;*/
float phi;


    phi=p_phi[0];

    phi = phi * PI / 180;


    if (euler_dim == 1)
    {
        /* rotation matrix */
        rm00 = cos(phi);
        rm01 = -sin(phi);
        rm10 = sin(phi);
        rm11 = cos(phi);
        /*printf("coeff: %f %f %f %f /n",rm00,rm01,rm10,rm11);*/
    }
    else
    {

        rm00 = p_phi[0];
        rm10 = p_phi[1];
        rm01 = p_phi[2];
        rm11 = p_phi[3];
    }

    sx_1 = sx - 1; /*width*/
    sy_1 = sy - 1; /*height*/
    for (pj=-py, j=0; j<sy; j++, pj++)
    {
        for (pi=-px, i=0; i<sx; i++, pi++)
        {


            /* transformation of coordinates (rotataion)*/
            r_x = px + rm00 * pi + rm10 * pj;
            if (r_x < 0 || r_x > sx_1 )
            {
                *rotimg++ = 0;    /* pixel not inside -> set 0*/
                continue;
            }
            r_y = py + rm01 * pi + rm11 * pj;
            if (r_y < 0 || r_y > sy_1 )
            {
                *rotimg++ = 0;
                continue;
            }

            /* Interpolation */

            *rotimg++ = interpolate2D(image,sx,sy, r_x,r_y,__TOM_INTERPOL_LINEAR);
        }
    }
}

void rot2dcubic(float *image,float *rotimg,long sx,long sy,float *p_phi,float px,float py,int euler_dim){
float rm00, rm01, rm10, rm11;   /* rot matrix */
float pi, pj;           /* coordinates according to pivot */
float r_x, r_y;         /* rotated pixel */
long  i, j;         /* loop variables */
/*long  sx_1, sy_1; *//* highest pixels  | image width and height */
/*float *imgp;*/
float phi;

long borderx1,borderx2,bordery1,bordery2;
phi=p_phi[0];

phi = phi * PI / 180;


    if (euler_dim == 1)
    {
        /* rotation matrix */
        rm00 = cos(phi);
        rm01 = -sin(phi);
        rm10 = sin(phi);
        rm11 = cos(phi);
        /*printf("coeff: %f %f %f %f /n",rm00,rm01,rm10,rm11);*/
    }
    else
    {

        rm00 = p_phi[0];
        rm10 = p_phi[1];
        rm01 = p_phi[2];
        rm11 = p_phi[3];
    }

    borderx1=1;
    bordery1=1;
    borderx2=sx-1;
    bordery2=sy-1;

    for (pj=-py, j=0; j<sy; j++, pj++)
    {

        for (pi=-px, i=0; i<sx; i++, pi++)
        {

            /*printf("x%d y%d %d %d\n",sx,sy,i,j);*/
            /* transformation of coordinates */
            /* rotate pixel of interest */
            r_x = px + rm00 * pi + rm10 * pj;
            r_y = py + rm01 * pi + rm11 * pj;

             /*if any point of interpolation square out of bounds, return linear2D interpolation or 0*/
            if (r_x < borderx1 || r_x > borderx2 || r_y < borderx1 || r_y > bordery2)
            {
                if(r_x < 1 || r_x > sx-1 || r_y < 1 || r_y > sy-1){
                    /*if not in image anymore, return 0*/
                    *rotimg++ = 0;
                    continue;
                }
                *rotimg++ = interpolate2D(image,sx,sy, r_x,r_y,__TOM_INTERPOL_LINEAR);
            }else{
                *rotimg++ = interpolate2D(image,sx,sy, r_x,r_y,__TOM_INTERPOL_CUBIC);
            }
        }
    }
}

/*cubic spline*/
void rot2dcspline (float *image,float *rotimg,long sx,long sy,float *p_phi,float px,float py,int euler_dim,float* fact)
{
float rm00, rm01, rm10, rm11;   /* rot matrix */
float pi, pj;           /* coordinates according to pivot */
float r_x, r_y;         /* rotated pixel */
long  i, j;         /* loop variables */
long  sx_1, sy_1;   /* highest pixels  | image width and height */
/*long  floorx, floory;*/       /* rotated coordinates as integer */
/*float vx1, vx2;*/     /* interp. parameter */
/*float vy1, vy2;*/     /*  ... */
/*float *imgp;*/
float phi;
float f1,f2,f3,f4;
long borderx1,borderx2,bordery1,bordery2;

phi=p_phi[0];

phi = phi * PI / 180;


if (euler_dim == 1)
{
    /* rotation matrix */
    rm00 = cos(phi);
    rm01 = -sin(phi);
    rm10 = sin(phi);
    rm11 = cos(phi);
    /*printf("coeff: %f %f %f %f /n",rm00,rm01,rm10,rm11);*/
}
else
{

    rm00 = p_phi[0];
    rm10 = p_phi[1];
    rm01 = p_phi[2];
    rm11 = p_phi[3];
}

/*cubic spline interpolation*/
    sx_1 = sx - 1;
    sy_1 = sy - 1;


    /*values for natural spline , quadratic*/

    if(fact != NULL){
            f1=fact[0];/* -0.1556; 1556 , 0940*/
            f2=fact[1];/* 0.3111; 3111 , 1880*/
            f3=fact[2];/* -0.0889; 0889 , 0342*/
            f4=fact[3];/* 0.0444; 0444 , 0171*/
    }else{
            f1= -0.1556;
            f2=0.3111;
            f3=-0.0889;
            f4=0.0444;
    }

    borderx1=2;
    bordery1=2;
    borderx2=sx-2;
    bordery2=sy-2;

    for (pj=-py, j=0; j<sy; j++, pj++)
    {

        for (pi=-px, i=0; i<sx; i++, pi++)
        {

            /*printf("x%d y%d %d %d\n",sx,sy,i,j);*/
            /* transformation of coordinates */
            /* rotate pixel of interest */
            r_x = px + rm00 * pi + rm10 * pj;
            r_y = py + rm01 * pi + rm11 * pj;

             /*if any point of interpolation square out of bounds, return linear2D interpolation or 0*/
            if (r_x < borderx1 || r_x > borderx2 || r_y < borderx1 || r_y > bordery2)
            {
                if(r_x < 1 || r_x > sx-1 || r_y < 1 || r_y > sy-1){
                    /*if not in image anymore, return 0*/
                    *rotimg++ = 0;
                    continue;
                }
                *rotimg++ = interpolate2D(image,sx,sy, r_x,r_y,__TOM_INTERPOL_LINEAR);
            }else{
                *rotimg++ = interpolate2D(image,sx,sy, r_x,r_y,__TOM_INTERPOL_CSPLINE);
            }
        }
    }
}



/* 3D Rotation */
void rot3dcspline (float *image,float *rotimg,long sx,long sy,long sz,float *p_euler_A,float px,float py,float pz,int euler_dim,float* fact)
{
long  sx_1, sy_1, sz_1; /* highest pixels */
long  sxy;
long  i, j, k;      /* loop variables */
long  pi, pj, pk;   /* loop variables */
float r_x, r_y, r_z;    /* rotated coordinates */
float rm00, rm01, rm02, rm10, rm11, rm12, rm20, rm21, rm22;     /* rotation matrix */
float sinphi=0., sinpsi=0., sintheta=0.; /* sin of rotation angles */
float cosphi=0., cospsi=0., costheta=0.; /* cos of rotation angles */
/*long  floorx, floory, floorz;*/       /* rotated coordinates as integer */
/*float vx1, vx2;*/
/*float vy1, vy2;*/         /* difference between floorx & r_x , ... */
/*float vz1, vz2;*/
/*long  iindex;*/           /* index of pixel in image vector */
float angles[] = { 0, 30, 45, 60, 90, 120, 135, 150, 180, 210, 225, 240, 270, 300, 315, 330 };
float angle_cos[16];
float angle_sin[16];
float phi,psi,theta;


float f1,f2,f3,f4;
long borderx1,borderx2,bordery1,bordery2,borderz1,borderz2;
/*long z;*/
/*float value[4];*/


phi=p_euler_A[0];
psi=p_euler_A[1];
theta=p_euler_A[2];

angle_cos[0]=1;
angle_cos[1]=sqrt(3)/2;
angle_cos[2]=sqrt(2)/2;
angle_cos[3]=0.5;
angle_cos[4]=0;
angle_cos[5]=-0.5;
angle_cos[6]=-sqrt(2)/2;
angle_cos[7]=-sqrt(3)/2;
angle_cos[8]=-1;
angle_cos[9]=-sqrt(3)/2;
angle_cos[10]=-sqrt(2)/2;
angle_cos[11]=-0.5;
angle_cos[12]=0;
angle_cos[13]=0.5;
angle_cos[14]=sqrt(2)/2;
angle_cos[15]=sqrt(3)/2;
angle_sin[0]=0;
angle_sin[1]=0.5;
angle_sin[2]=sqrt(2)/2;
angle_sin[3]=sqrt(3)/2;
angle_sin[4]=1;
angle_sin[5]=sqrt(3)/2;
angle_sin[6]=sqrt(2)/2;
angle_sin[7]=0.5;
angle_sin[8]=0;
angle_sin[9]=-0.5;
angle_sin[10]=-sqrt(2)/2;
angle_sin[11]=-sqrt(3)/2;
angle_sin[12]=-1;
angle_sin[13]=-sqrt(3)/2;
angle_sin[14]=-sqrt(2)/2;
angle_sin[15]=-0.5;

sx_1 = (long)(sx - 0.5);
sy_1 = (long)(sy - 0.5);
sz_1 = (long)(sz - 0.5);

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == phi )
    {
        cosphi = angle_cos[i];
        sinphi = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    phi = phi * PI / 180;
    cosphi=cos(phi);
    sinphi=sin(phi);
}

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == psi )
    {
        cospsi = angle_cos[i];
        sinpsi = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    psi = psi * PI / 180;
    cospsi=cos(psi);
    sinpsi=sin(psi);
}

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == theta )
    {
        costheta = angle_cos[i];
        sintheta = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    theta = theta * PI / 180;
    costheta=cos(theta);
    sintheta=sin(theta);
}

/* calculation of rotation matrix */

if (euler_dim == 1)
{
    rm00=cospsi*cosphi-costheta*sinpsi*sinphi;
    rm10=sinpsi*cosphi+costheta*cospsi*sinphi;
    rm20=sintheta*sinphi;
    rm01=-cospsi*sinphi-costheta*sinpsi*cosphi;
    rm11=-sinpsi*sinphi+costheta*cospsi*cosphi;
    rm21=sintheta*cosphi;
    rm02=sintheta*sinpsi;
    rm12=-sintheta*cospsi;
    rm22=costheta;
    /*printf("coeff: %f %f %f  %f %f %f  %f %f %f /n",rm00,rm10,rm20,rm01,rm11,rm21,rm02,rm12,rm22);*/
}
else
{
    rm00=p_euler_A[0];
    rm10=p_euler_A[1];
    rm20=p_euler_A[2];
    rm01=p_euler_A[3];
    rm11=p_euler_A[4];
    rm21=p_euler_A[5];
    rm02=p_euler_A[6];
    rm12=p_euler_A[7];
    rm22=p_euler_A[8];
}


    sx_1 = sx - 1;
    sy_1 = sy - 1;
    sz_1 = sz - 1;
    sxy = sx * sy;


    borderx1=2;
    bordery1=2;
    borderz1=2;
    borderx2=sx-2;
    bordery2=sy-2;
    borderz2=sz-2;

    /*values for natural spline , quadratic*/
    if(fact != NULL){
            f1=fact[0];/* -0.1556; 1556 , 0940*/
            f2=fact[1];/* 0.3111; 3111 , 1880*/
            f3=fact[2];/* -0.0889; 0889 , 0342*/
            f4=fact[3];/* 0.0444; 0444 , 0171*/
    }else{
            f1= -0.1556;
            f2=0.3111;
            f3=-0.0889;
            f4=0.0444;
    }
    for (k=0; k < sz; k++)
    {

        for (j=0; j < sy; j++)
        {
            for (i=0; i < sx; i++)
            {
                pi = (long)(i-px);
                pj = (long)(j-py);
                pk = (long)(k-pz);

                /* transformation of coordinates */
                r_x = px + rm00 * pi + rm10 * pj + rm20 * pk;
                if (r_x < 0 || r_x > sx_1 )
                {
                    *rotimg++ = 0;   /* this pixel was not inside the image */
                    continue;
                }
                r_y = py + rm01 * pi + rm11 * pj + rm21 * pk;
                if (r_y < 0 || r_y > sy_1 )
                {
                    *rotimg++ = 0;
                    continue;
                }
                r_z = pz + rm02 * pi + rm12 * pj + rm22 * pk;
                if (r_z < 0 || r_z > sz_1 )
                {
                    *rotimg++ = 0;
                    continue;
                }

                /*If pixel is out of spline bounds, return linear3d interpolation*/
                if (r_x < borderx1 || r_x > borderx2 || r_y < bordery1 || r_y > bordery2 || r_z < borderz1 || r_z > borderz2)
                {
                    *rotimg++ = interpolate3D(image, sx, sy, sz, r_x,r_y, r_z, __TOM_INTERPOL_LINEAR);
                    continue;
                }

                /*CSpline interpolation here*/
                *rotimg++ = interpolate3D(image, sx, sy, sz, r_x,r_y, r_z, __TOM_INTERPOL_CSPLINE);

            }
        }
    }

}

void rot3dcubic (float *image,float *rotimg,long sx,long sy,long sz,float *p_euler_A,float px,float py,float pz,int euler_dim)
{
long  sx_1, sy_1, sz_1; /* highest pixels */
long  sxy;
long  i, j, k;      /* loop variables */
long  pi, pj, pk;   /* loop variables */
float r_x, r_y, r_z;    /* rotated coordinates */
float rm00, rm01, rm02, rm10, rm11, rm12, rm20, rm21, rm22;     /* rotation matrix */
float sinphi=0., sinpsi=0., sintheta=0.; /* sin of rotation angles */
float cosphi=0., cospsi=0., costheta=0.; /* cos of rotation angles */
/*long  floorx, floory, floorz;*/       /* rotated coordinates as integer */
/*float vx1, vx2;*/
/*float vy1, vy2;*/         /* difference between floorx & r_x , ... */
/*float vz1, vz2;*/
/*long  iindex;*/           /* index of pixel in image vector */
float angles[] = { 0, 30, 45, 60, 90, 120, 135, 150, 180, 210, 225, 240, 270, 300, 315, 330 };
float angle_cos[16];
float angle_sin[16];
float phi,psi,theta;


/*float f1,f2,f3,f4;*/
long borderx1,borderx2,bordery1,bordery2,borderz1,borderz2;
/*long z;*/
/*float value4];*/


phi=p_euler_A[0];
psi=p_euler_A[1];
theta=p_euler_A[2];

angle_cos[0]=1;
angle_cos[1]=sqrt(3)/2;
angle_cos[2]=sqrt(2)/2;
angle_cos[3]=0.5;
angle_cos[4]=0;
angle_cos[5]=-0.5;
angle_cos[6]=-sqrt(2)/2;
angle_cos[7]=-sqrt(3)/2;
angle_cos[8]=-1;
angle_cos[9]=-sqrt(3)/2;
angle_cos[10]=-sqrt(2)/2;
angle_cos[11]=-0.5;
angle_cos[12]=0;
angle_cos[13]=0.5;
angle_cos[14]=sqrt(2)/2;
angle_cos[15]=sqrt(3)/2;
angle_sin[0]=0;
angle_sin[1]=0.5;
angle_sin[2]=sqrt(2)/2;
angle_sin[3]=sqrt(3)/2;
angle_sin[4]=1;
angle_sin[5]=sqrt(3)/2;
angle_sin[6]=sqrt(2)/2;
angle_sin[7]=0.5;
angle_sin[8]=0;
angle_sin[9]=-0.5;
angle_sin[10]=-sqrt(2)/2;
angle_sin[11]=-sqrt(3)/2;
angle_sin[12]=-1;
angle_sin[13]=-sqrt(3)/2;
angle_sin[14]=-sqrt(2)/2;
angle_sin[15]=-0.5;

sx_1 = (long)(sx - 0.5);
sy_1 = (long)(sy - 0.5);
sz_1 = (long)(sz - 0.5);

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == phi )
    {
        cosphi = angle_cos[i];
        sinphi = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    phi = phi * PI / 180;
    cosphi=cos(phi);
    sinphi=sin(phi);
}

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == psi )
    {
        cospsi = angle_cos[i];
        sinpsi = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    psi = psi * PI / 180;
    cospsi=cos(psi);
    sinpsi=sin(psi);
}

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == theta )
    {
        costheta = angle_cos[i];
        sintheta = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    theta = theta * PI / 180;
    costheta=cos(theta);
    sintheta=sin(theta);
}

/* calculation of rotation matrix */

if (euler_dim == 1)
{
    rm00=cospsi*cosphi-costheta*sinpsi*sinphi;
    rm10=sinpsi*cosphi+costheta*cospsi*sinphi;
    rm20=sintheta*sinphi;
    rm01=-cospsi*sinphi-costheta*sinpsi*cosphi;
    rm11=-sinpsi*sinphi+costheta*cospsi*cosphi;
    rm21=sintheta*cosphi;
    rm02=sintheta*sinpsi;
    rm12=-sintheta*cospsi;
    rm22=costheta;
    /*printf("coeff: %f %f %f  %f %f %f  %f %f %f /n",rm00,rm10,rm20,rm01,rm11,rm21,rm02,rm12,rm22);*/
}
else
{
    rm00=p_euler_A[0];
    rm10=p_euler_A[1];
    rm20=p_euler_A[2];
    rm01=p_euler_A[3];
    rm11=p_euler_A[4];
    rm21=p_euler_A[5];
    rm02=p_euler_A[6];
    rm12=p_euler_A[7];
    rm22=p_euler_A[8];
}


    sx_1 = sx - 1;
    sy_1 = sy - 1;
    sz_1 = sz - 1;
    sxy = sx * sy;


    borderx1=1;
    bordery1=1;
    borderz1=2;
    borderx2=sx-1;
    bordery2=sy-1;
    borderz2=sz-2;


    for (k=0; k < sz; k++)
    {

        for (j=0; j < sy; j++)
        {
            for (i=0; i < sx; i++)
            {
                pi = (long)(i-px);
                pj = (long)(j-py);
                pk = (long)(k-pz);

                /* transformation of coordinates */
                r_x = px + rm00 * pi + rm10 * pj + rm20 * pk;
                if (r_x < 0 || r_x > sx_1 )
                {
                    *rotimg++ = 0;   /* this pixel was not inside the image */
                    continue;
                }
                r_y = py + rm01 * pi + rm11 * pj + rm21 * pk;
                if (r_y < 0 || r_y > sy_1 )
                {
                    *rotimg++ = 0;
                    continue;
                }
                r_z = pz + rm02 * pi + rm12 * pj + rm22 * pk;
                if (r_z < 0 || r_z > sz_1 )
                {
                    *rotimg++ = 0;
                    continue;
                }

                /*If pixel is out of spline bounds, return linear3d interpolation*/
                if (r_x < borderx1 || r_x > borderx2 || r_y < bordery1 || r_y > bordery2 || r_z < borderz1 || r_z > borderz2)
                {
                    *rotimg++ = interpolate3D(image, sx, sy, sz, r_x,r_y, r_z, __TOM_INTERPOL_LINEAR);
                    continue;
                }

                /*CSpline interpolation here*/
                *rotimg++ = interpolate3D(image, sx, sy, sz, r_x,r_y, r_z, __TOM_INTERPOL_CUBIC);

            }
        }
    }
}


/* 3D Rotation */
void rot3dlin (float *image,float *rotimg,long sx,long sy,long sz,float *p_euler_A,float px,float py,float pz,int euler_dim)
{
long  sx_1, sy_1, sz_1; /* highest pixels */
long  sxy;
long  i, j, k;      /* loop variables */
long  pi, pj, pk;   /* loop variables */
float r_x, r_y, r_z;    /* rotated coordinates */
float rm00, rm01, rm02, rm10, rm11, rm12, rm20, rm21, rm22;     /* rotation matrix */
float sinphi=0., sinpsi=0., sintheta=0.; /* sin of rotation angles */
float cosphi=0., cospsi=0., costheta=0.; /* cos of rotation angles */
/*long  floorx, floory, floorz;*/       /* rotated coordinates as integer */
/*float vx1, vx2;*/
/*float vy1, vy2;*/         /* difference between floorx & r_x , ... */
/*float vz1, vz2;*/
long  index1, index2, index3, index4, index5, index6;
float angles[] = { 0, 30, 45, 60, 90, 120, 135, 150, 180, 210, 225, 240, 270, 300, 315, 330 };
float angle_cos[16];
float angle_sin[16];
float phi,psi,theta;


phi=p_euler_A[0];
psi=p_euler_A[1];
theta=p_euler_A[2];

angle_cos[0]=1;
angle_cos[1]=sqrt(3)/2;
angle_cos[2]=sqrt(2)/2;
angle_cos[3]=0.5;
angle_cos[4]=0;
angle_cos[5]=-0.5;
angle_cos[6]=-sqrt(2)/2;
angle_cos[7]=-sqrt(3)/2;
angle_cos[8]=-1;
angle_cos[9]=-sqrt(3)/2;
angle_cos[10]=-sqrt(2)/2;
angle_cos[11]=-0.5;
angle_cos[12]=0;
angle_cos[13]=0.5;
angle_cos[14]=sqrt(2)/2;
angle_cos[15]=sqrt(3)/2;
angle_sin[0]=0;
angle_sin[1]=0.5;
angle_sin[2]=sqrt(2)/2;
angle_sin[3]=sqrt(3)/2;
angle_sin[4]=1;
angle_sin[5]=sqrt(3)/2;
angle_sin[6]=sqrt(2)/2;
angle_sin[7]=0.5;
angle_sin[8]=0;
angle_sin[9]=-0.5;
angle_sin[10]=-sqrt(2)/2;
angle_sin[11]=-sqrt(3)/2;
angle_sin[12]=-1;
angle_sin[13]=-sqrt(3)/2;
angle_sin[14]=-sqrt(2)/2;
angle_sin[15]=-0.5;

sx_1 = (long)(sx - 0.5);
sy_1 = (long)(sy - 0.5);
sz_1 = (long)(sz - 0.5);

#define cpu_time_used (((double) (end - start)) / CLOCKS_PER_SEC)


for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == phi )
    {
        cosphi = angle_cos[i];
        sinphi = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    phi = phi * PI / 180;
    cosphi=cos(phi);
    sinphi=sin(phi);
}

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == psi )
    {
        cospsi = angle_cos[i];
        sinpsi = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    psi = psi * PI / 180;
    cospsi=cos(psi);
    sinpsi=sin(psi);
}

for (i=0, j=0 ; i<16; i++)
{
    if (angles[i] == theta )
    {
        costheta = angle_cos[i];
        sintheta = angle_sin[i];
        j = 1;
    }
}

if (j < 1)
{
    theta = theta * PI / 180;
    costheta=cos(theta);
    sintheta=sin(theta);
}

/* calculation of rotation matrix */

if (euler_dim == 1)
{
    rm00=cospsi*cosphi-costheta*sinpsi*sinphi;
    rm10=sinpsi*cosphi+costheta*cospsi*sinphi;
    rm20=sintheta*sinphi;
    rm01=-cospsi*sinphi-costheta*sinpsi*cosphi;
    rm11=-sinpsi*sinphi+costheta*cospsi*cosphi;
    rm21=sintheta*cosphi;
    rm02=sintheta*sinpsi;
    rm12=-sintheta*cospsi;
    rm22=costheta;
    /*printf("coeff: %f %f %f  %f %f %f  %f %f %f /n",rm00,rm10,rm20,rm01,rm11,rm21,rm02,rm12,rm22);*/
}
else
{
    rm00=p_euler_A[0];
    rm10=p_euler_A[1];
    rm20=p_euler_A[2];
    rm01=p_euler_A[3];
    rm11=p_euler_A[4];
    rm21=p_euler_A[5];
    rm02=p_euler_A[6];
    rm12=p_euler_A[7];
    rm22=p_euler_A[8];
}



    sx_1 = sx - 1;
    sy_1 = sy - 1;
    sz_1 = sz - 1;
    sxy = sx * sy;
    index1 = sx;
    index2 = sx + 1;
    index3 = sxy;
    index4 = sxy + 1;
    index5 = sx + sxy;
    index6 = sx + sxy + 1;



    for (k=0; k < sz; k++)
    {

        for (j=0; j < sy; j++)
        {
            for (i=0; i < sx; i++)
            {
                pi = (long)(i-px);
                pj = (long)(j-py);
                pk = (long)(k-pz);
                /* transformation of coordinates */
                r_x = px + rm00 * pi + rm10 * pj + rm20 * pk;
                if (r_x < 0 || r_x > sx_1 )
                {
                    *rotimg++ = 0;   /* this pixel was not inside the image */
                    continue;
                }
                r_y = py + rm01 * pi + rm11 * pj + rm21 * pk;
                if (r_y < 0 || r_y > sy_1 )
                {
                    *rotimg++ = 0;
                    continue;
                }
                r_z = pz + rm02 * pi + rm12 * pj + rm22 * pk;
                if (r_z < 0 || r_z > sz_1 )
                {
                    *rotimg++ = 0;
                    continue;
                }


                *rotimg++ = interpolate3D(image, sx, sy, sz, r_x,r_y,r_z,__TOM_INTERPOL_LINEAR);
            }
        }
    }


}



