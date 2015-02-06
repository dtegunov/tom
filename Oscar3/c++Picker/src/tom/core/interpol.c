#include "tom/core/interpol.h"
/***********************************************************************//**
 * \file interpol.c
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    23.10.2007
 **************************************************************************/
/*Library for interpolation functions for 2D and 3D */
/*implemented interpolation types:*/
/**********************************/
/*Nearest Neighbour*/
/*Linear*/
/*Cubic*/
/*CSpline*/
/**********************************/
/**/
/*When using this library, mind the value access in the data volume at the ends.*/
/**/
/*Linear uses a 2x2(x2) volume for interpolation*/
/*Cubic uses a 3x3(x3) volume for interpolation*/
/*CSpline uses a 4x4(x4) volume for interpolation*/
/**/
/*Documentation for the interpolation methods:*/
/*Nearest - none*/
/*Linear - C Crane - Image Processing*/
/*Cubic - http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html*/
/*CSpline - *http://mathworld.wolfram.com/CubicSpline.html*/











/*macro for cubic interpolation*/
/*calculates the polynom value P_i(x)*y_i according to */
/*http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html*/
#define CUB_INT(x,y,xi,x2,x3,val){\
val = y*(x-x2)/(xi-x2)*(x-x3)/(xi-x3);\
}



/*prestep for cspline interpolation for the calculation of D1 and D2 according to source*/
/*http://mathworld.wolfram.com/CubicSpline.html*/
/*macro for spline interpolation in a row*/
/*3*(f1*(v2-v1)+f2*(v3-v1)+f3*(v4-v2)+f4*(v4-v3))*/
#define CSPL_INT(f1,f2,f3,f4,a,b,c,d,val){\
val = 3*(f1*(b-a)+f2*(c-a)+f3*(d-b)+f4*(d-c));\
}


/*macro to calculate the result of cspline interpolation using 4 features*/
/*c2 and c3 are calculated by CSPL_INT*/
/*value = c2+D1*(r_y-py6)+(3*(c3-c2)-2*D1-D2)*(r_y-py6)*(r_y-py6)+(2*(c2-c3)+D1+D2)*(r_y-py6)*(r_y-py6)*(r_y-py6);*/
/*c2~c2, c3~c3,D1~D1,D2~D2*,off~(r_y-py6) - distance from nearest point*/
#define CSPL_CALC(c2,c3,D1,D2,off,val){\
val = c2+D1*off+(3*(c3-c2)-2*D1-D2)*off*off+(2*(c2-c3)+D1+D2)*off*off*off;\
}




float interpolate2D(const float *image,long sizeX, long sizeY,float x,float y,short method){

long  floorx, floory;		/* rotated coordinates as integer */
float vx1, vx2;		/* interp. parameter */
float vy1, vy2;		/*  ... */
float AA, BB;	    /*  ... */
const float *imgp;
long px1,px2,px3,px4,px5,px6,px7,px8,px9,px10,px11,px12,px13,px14,px15,px16; /*coordinates of 4x4 square round interpolation position for spline & bicubic */
long py1,py2,py3,py4,py5,py6,py7,py8,py9,py10,py11,py12,py13,py14,py15,py16; /*coordinates of 4x4 square round interpolation position for spline & bicubic */
float v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16; /*storage for pixel values*/
float c1,c2,c3,c4; /*storage for cubic interpolation values*/
float xtmp, D1, D2;
/*float a,b,c,d,cub,phi;*/
float f1,f2,f3,f4;
/*long borderx1,borderx2,bordery1,bordery2;*/

if(method == __TOM_INTERPOL_NEAREST){

    return image[(long)x + sizeX*((long)y)];

}else if(method == __TOM_INTERPOL_LINEAR){
			/*Linear Interpolation */
			floorx = (long)x;
			vx2 = x - floorx;
			vx1 = 1 - vx2;
			floory = (long)y;
			vy2 = y - floory;
			vy1 = 1 - vy2;
			imgp = &image[floorx + sizeX*floory];
			if (x<sizeX-1)    			/* not last x pixel, */
				AA = *imgp+(imgp[1]-*imgp)*vx2;	/* interpolation */
			else           	/* last x pixel, no interpolation in x possible */
				AA = *imgp;
			if (y<sizeY-1)
			{				/* not last y pixel, */
				if (x<sizeX-1)			/* not last x pixel, */
					BB = imgp[sizeX]+(imgp[sizeX+1]-imgp[sizeX])*vx2;/* interpolation */
				else
					BB = imgp[sizeX];  /* last x pixel, no interpolation in x possible */

			}
			else
			{
				BB = 0; /* last y pixel, no interpolation in y possible */
			}
			return AA + (BB - AA) * vy2;

}else if(method == __TOM_INTERPOL_CUBIC){
            /*interpolation square
            p1  p2  p3
            p4  p5  p6
                  P
            p7  p8  p9
            */

            px5 = (long)x;
            py5 = (long)y;

            px1 = px5 -1;
            py1 = py5 -1;

            px2 = px5;
            py2 = py5 -1;

            px3 = px5 +1;
            py3 = py5 -1;

            px4 = px5 -1;
            py4 = py5;

            px6 = px5 +1;
            py6 = py5;

            px7 = px5 -1;
            py7 = py5 +1;

            px8 = px5;
            py8 = py5+1;

            px9 = px5+1;
            py9 = py5+1;

            /*load the pixel values*/
            v1 = image[px1 + sizeX*py1];
            v2 = image[px2 + sizeX*py2];
            v3 = image[px3 + sizeX*py3];
            v4 = image[px4 + sizeX*py4];
            v5 = image[px5 + sizeX*py5];
            v6 = image[px6 + sizeX*py6];
            v7 = image[px7 + sizeX*py7];
            v8 = image[px8 + sizeX*py8];
            v9 = image[px9 + sizeX*py9];


            /*interpolate first row 1 2 3*/
            CUB_INT(x,v1,px1,px2,px3,c1);
            CUB_INT(x,v2,px2,px3,px1,c2);
            CUB_INT(x,v3,px3,px1,px2,c3);
            /*store values into v1*/
            /*printf("%f %f %f",c1,c2,c3);*/
            v1 = c1+c2+c3;
            /*same for the next rows*/
            CUB_INT(x,v4,px4,px5,px6,c1);
            CUB_INT(x,v5,px5,px6,px4,c2);
            CUB_INT(x,v6,px6,px5,px4,c3);
            /*printf("%f %f %f",c1,c2,c3);*/
            v2 = c1+c2+c3;

            CUB_INT(x,v7,px7,px8,px9,c1);
            CUB_INT(x,v8,px8,px9,px7,c2);
            CUB_INT(x,v9,px9,px7,px8,c3);

            v3 = c1+c2+c3;

            /*interpolate col 2 5 8 in y direction*/
            CUB_INT(y,v1,py2,py5,py8,c1);
            CUB_INT(y,v2,py5,py8,py2,c2);
            CUB_INT(y,v3,py8,py2,py5,c3);
            /*return value*/

            return c1+c2+c3;
}else if(method == __TOM_INTERPOL_CSPLINE){

            f1= -0.1556;
            f2=0.3111;
            f3=-0.0889;
            f4=0.0444;

            /*interpolation square

             p1   p2   p3 p4
             p5   p6   p7 p8
                     P
             p9   p10  p11 p12
             p13  p14  p15 p16

             */

            px6 = (long)x;
            py6 = (long)y;

            px1 = px6-1;
            py1 = py6-1;

            px2 = px6;
            py2 = py6-1;

            px3 = px6+1;
            py3 = py6+1;

            px4 = px6+2;
            py4 = py6-1;

            px5 = px6-1;
            py5 = py6;

            px7 = px6+1;
            py7 = py6;

            px8 = px6+2;
            py8 = py6;

            px9 = px6-1;
            py9 = py6+1;

            px10 = px6;
            py10 = py6+1;

            px11 = px6+1;
            py11 = py6+1;

            px12 = px6+2;
            py12 = py6+1;

            px13 = px6-1;
            py13 = py6+2;

            px14 = px6;
            py14 = py6+2;

            px15 = px6+1;
            py15 = py6+2;

            px16 = px6+2;
            py16 = py6+2;

            /*load the pixel values*/
            v1 = image[px1 + sizeX*py1];
            v2 = image[px2 + sizeX*py2];
            v3 = image[px3 + sizeX*py3];
            v4 = image[px4 + sizeX*py4];
            v5 = image[px5 + sizeX*py5];
            v6 = image[px6 + sizeX*py6];
            v7 = image[px7 + sizeX*py7];
            v8 = image[px8 + sizeX*py8];
            v9 = image[px9 + sizeX*py9];
            v10 = image[px10 + sizeX*py10];
            v11 = image[px11 + sizeX*py11];
            v12 = image[px12 + sizeX*py12];
            v13 = image[px13 + sizeX*py13];
            v14 = image[px14 + sizeX*py14];
            v15 = image[px15 + sizeX*py15];
            v16 = image[px16 + sizeX*py16];

            /*distance in x direction from pixel of interest to lower pixel in image */
            xtmp = (x-px6);

            /*calculate spline value for line 1 2 3 4 above pixel of interest */
            CSPL_INT(f1,f2,f3,f4,v1,v2,v3,v4,D1);
            CSPL_INT(f4,f3,f2,f1,v1,v2,v3,v4,D2);
            CSPL_CALC(v2,v3,D1,D2,xtmp,c1);


            /*calculate spline value for line 5 6 7 8 above pixel of interest */
            CSPL_INT(f1,f2,f3,f4,v5,v6,v7,v8,D1);
            CSPL_INT(f4,f3,f2,f1,v5,v6,v7,v8,D2);
            CSPL_CALC(v6,v7,D1,D2,xtmp,c2);

            /*calculate spline value for line 9 10 11 12 above pixel of interest */
            CSPL_INT(f1,f2,f3,f4,v9,v10,v11,v12,D1);
            CSPL_INT(f4,f3,f2,f1,v9,v10,v11,v12,D2);
            CSPL_CALC(v10,v11,D1,D2,xtmp,c3);

            /*calculate spline value for line 13 14 15 16 above pixel of interest */
            CSPL_INT(f1,f2,f3,f4,v13,v14,v15,v16,D1);
            CSPL_INT(f4,f3,f2,f1,v13,v14,v15,v16,D2);
            CSPL_CALC(v14,v15,D1,D2,xtmp,c4);

            /*finaly, calculate spline into y direction and return*/
            CSPL_INT(f1,f2,f3,f4,c1,c2,c3,c4,D1);
            CSPL_INT(f4,f3,f2,f1,c1,c2,c3,c4,D2);

            return c2+D1*(y-py6)+(3*(c3-c2)-2*D1-D2)*(y-py6)*(y-py6)+(2*(c2-c3)+D1+D2)*(y-py6)*(y-py6)*(y-py6);
        }
return 0.;
}


float interpolate3D(const float* image, long sizeX, long sizeY, long sizeZ, float x,float y, float z, short method){
    long floorx,floory,floorz;
    float vx1,vx2,vy1,vy2,vz1,vz2;
    long iindex;
    float img0, img1, img2, img3, img4, img5, img6;
    long  index1, index2, index3, index4, index5, index6;
    float AA,BB,CC,DD;
    long runz;
    long px1,px2,px3,px4,px5,px6,px7,px8,px9,px10,px11,px12,px13,px14,px15,px16; /*coordinates of 4x4 square round interpolation position for spline & bicubic */
    long py1,py2,py3,py4,py5,py6,py7,py8,py9,py10,py11,py12,py13,py14,py15,py16; /*coordinates of 4x4 square round interpolation position for spline & bicubic */
    float v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16; /*storage for pixel values*/
    float c1,c2,c3,c4; /*storage for cubic interpolation values*/
    float D1,D2;
    float f1,f2,f3,f4;
    float value[4];
    /*short m;*/


    if(method == __TOM_INTERPOL_NEAREST){

        return image[((long)x) + sizeX*((long)y)+ sizeX*sizeY*((long)z)];

    }else if(method == __TOM_INTERPOL_LINEAR){
                index1 = sizeX;
            	index2 = sizeX + 1;
            	index3 = sizeX*sizeY;
            	index4 = index3 + 1;
            	index5 = sizeX + index3;
            	index6 = index5 + 1;

				/* Interpolation */
				floorx = (long)x;
				vx2 = x - floorx;
				vx1 = 1 - vx2;
				floory = (long)y;
				vy2 = y - floory;
				vy1 = 1 - vy2;
				floorz = (long)z;
				vz2 = z - floorz;
				vz1 = 1 - vz2;

				/* the following section detects border pixels to avoid exceeding dimensions */
				iindex = floorx + floory * sizeX + floorz * sizeX*sizeY;
				img0 = img1 = img2 = img3 = img4 = img5 = img6 = 1;
				if (floorx == sizeX-1)
				{
					img0 = img2 = img4 = img6 = 0;
				}
				if (floory == sizeY-1)
				{
					img1 = img2 = img5 = img6 = 0;
				}
				if (floorz == sizeZ-1)
				{
					img3 = img4 = img5 = img6 = 0;
				}

				if (img0==1)
				{
					img0=image[iindex + 1];
				}
				if (img1==1)
				{
					img1=image[iindex + index1];
				}
				if (img2==1)
				{
					img2=image[iindex + index2];
				}
				if (img3==1)
				{
					img3=image[iindex + index3];
				}
				if (img4==1)
				{
					img4=image[iindex + index4];
				}
				if (img5==1)
				{
					img5=image[iindex + index5];
				}
				if (img6==1)
				{
					img6=image[iindex + index6];
				}

				/* interpolation */
				AA = image[iindex] + (img0 - image[iindex]) * vx2;
				BB = img1 * vx1 + img2 * vx2;
				CC = img3 * vx1 + img4 * vx2;
				DD = img5 * vx1 + img6 * vx2;
				return  (AA * vy1 + BB * vy2) * vz1 + (CC * vy1 + DD * vy2) * vz2;
    }else if(method == __TOM_INTERPOL_CUBIC){

            /*interpolation square - 3 layers for cube
            p1  p2  p3
            p4  p5  p6
                  P
            p7  p8  p9
            */
            floorx = (long)x;
            floory = (long)y;
            floorz = (long)z;

            px5 = (long)x;
            py5 = (long)y;

            px1 = px5 -1;
            py1 = py5 -1;

            px2 = px5;
            py2 = py5 -1;

            px3 = px5 +1;
            py3 = py5 -1;

            px4 = px5 -1;
            py4 = py5;

            px6 = px5 +1;
            py6 = py5;

            px7 = px5 -1;
            py7 = py5 +1;

            px8 = px5;
            py8 = py5+1;

            px9 = px5+1;
            py9 = py5+1;
            for(runz=0;runz<3;runz++){
                value[runz]=0;
                /*load the pixel values*/
                v1 = image[px1 + sizeX*py1+ sizeX*sizeY*((long)z+runz-1)];
                v2 = image[px2 + sizeX*py2+ sizeX*sizeY*((long)z+runz-1)];
                v3 = image[px3 + sizeX*py3+ sizeX*sizeY*((long)z+runz-1)];
                v4 = image[px4 + sizeX*py4+ sizeX*sizeY*((long)z+runz-1)];
                v5 = image[px5 + sizeX*py5+ sizeX*sizeY*((long)z+runz-1)];
                v6 = image[px6 + sizeX*py6+ sizeX*sizeY*((long)z+runz-1)];
                v7 = image[px7 + sizeX*py7+ sizeX*sizeY*((long)z+runz-1)];
                v8 = image[px8 + sizeX*py8+ sizeX*sizeY*((long)z+runz-1)];
                v9 = image[px9 + sizeX*py9+ sizeX*sizeY*((long)z+runz-1)];


                /*interpolate first row 1 2 3*/
                CUB_INT(x,v1,px1,px2,px3,c1);
                CUB_INT(x,v2,px2,px3,px1,c2);
                CUB_INT(x,v3,px3,px1,px2,c3);
                /*store values into v1*/
                v1 = c1+c2+c3;
                /*same for the next rows*/
                CUB_INT(x,v4,px4,px5,px6,c1);
                CUB_INT(x,v5,px5,px6,px4,c2);
                CUB_INT(x,v6,px6,px5,px4,c3);

                v2 = c1+c2+c3;
                CUB_INT(x,v7,px7,px8,px9,c1);
                CUB_INT(x,v8,px8,px9,px7,c2);
                CUB_INT(x,v9,px9,px7,px8,c3);

                v3 = c1+c2+c3;
                /*interpolate col 2 5 8 in y direction*/
                CUB_INT(y,v1,py2,py5,py8,c1);
                CUB_INT(y,v2,py5,py8,py2,c2);
                CUB_INT(y,v3,py8,py2,py5,c3);
                /*return value*/
                value[runz]=c1+c2+c3;
            }
            px1 = floorz-1;
            px2 = floorz;
            px3 = floorz+1;

            CUB_INT(z,value[0],px1,px2,px3,c1);
            CUB_INT(z,value[1],px2,px3,px1,c2);
            CUB_INT(z,value[2],px3,px1,px2,c3);


            return c1+c2+c3;

    }else if(method == __TOM_INTERPOL_CSPLINE){


        f1= -0.1556;
        f2=0.3111;
        f3=-0.0889;
        f4=0.0444;



        floorx = (long)x;
        floory = (long)y;
        floorz = (long)z;


                /*init variables*/


                /*interpolation square

                 p1   p2   p3 p4
                 p5   p6   p7 p8
                         P
                 p9   p10  p11 p12
                 p13  p14  p15 p16

                 */

         px6 = floorx;
         py6 = floory;

         px1 = px6-1;
         py1 = py6-1;

         px2 = px6;
         py2 = py6-1;

         px3 = px6+1;
         py3 = py6+1;

         px4 = px6+2;
         py4 = py6-1;

         px5 = px6-1;
         py5 = py6;

         px7 = px6+1;
         py7 = py6;

         px8 = px6+2;
         py8 = py6;

         px9 = px6-1;
         py9 = py6+1;

         px10 = px6;
         py10 = py6+1;

         px11 = px6+1;
         py11 = py6+1;

         px12 = px6+2;
         py12 = py6+1;

         px13 = px6-1;
         py13 = py6+2;

         px14 = px6;
         py14 = py6+2;

         px15 = px6+1;
         py15 = py6+2;

         px16 = px6+2;
         py16 = py6+2;

         for(runz=0;runz<4;runz++){

                value[runz]=0;
                /*load the pixel values*/
                v1 = image[px1 + sizeX*py1+ sizeX*sizeY*((long)z+runz-2)]; /*-2 because we take two left of z coordinate and 2 right, but z goes 0 - 3*/
                v2 = image[px2 + sizeX*py2+ sizeX*sizeY*((long)z+runz-2)];
                v3 = image[px3 + sizeX*py3+ sizeX*sizeY*((long)z+runz-2)];
                v4 = image[px4 + sizeX*py4+ sizeX*sizeY*((long)z+runz-2)];
                v5 = image[px5 + sizeX*py5+ sizeX*sizeY*((long)z+runz-2)];
                v6 = image[px6 + sizeX*py6+ sizeX*sizeY*((long)z+runz-2)];
                v7 = image[px7 + sizeX*py7+ sizeX*sizeY*((long)z+runz-2)];
                v8 = image[px8 + sizeX*py8+ sizeX*sizeY*((long)z+runz-2)];
                v9 = image[px9 + sizeX*py9+ sizeX*sizeY*((long)z+runz-2)];
                v10 = image[px10 + sizeX*py10+ sizeX*sizeY*((long)z+runz-2)];
                v11 = image[px11 + sizeX*py11+ sizeX*sizeY*((long)z+runz-2)];
                v12 = image[px12 + sizeX*py12+ sizeX*sizeY*((long)z+runz-2)];
                v13 = image[px13 + sizeX*py13+ sizeX*sizeY*((long)z+runz-2)];
                v14 = image[px14 + sizeX*py14+ sizeX*sizeY*((long)z+runz-2)];
                v15 = image[px15 + sizeX*py15+ sizeX*sizeY*((long)z+runz-2)];
                v16 = image[px16 + sizeX*py16+ sizeX*sizeY*((long)z+runz-2)];
                    /*calculate spline value of x row*/
                    /*interpolate into the y col*/
                    /*store this value somewhere*/


                /*calculate spline value for line 1 2 3 4 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v1,v2,v3,v4,D1);
                CSPL_INT(f4,f3,f2,f1,v1,v2,v3,v4,D2);
                CSPL_CALC(v2,v3,D1,D2,(x-floorx),c1);

                /*calculate spline value for line 5 6 7 8 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v5,v6,v7,v8,D1);
                CSPL_INT(f4,f3,f2,f1,v5,v6,v7,v8,D2);
                CSPL_CALC(v6,v7,D1,D2,(x-floorx),c2);

                /*calculate spline value for line 9 10 11 12 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v9,v10,v11,v12,D1);
                CSPL_INT(f4,f3,f2,f1,v9,v10,v11,v12,D2);
                CSPL_CALC(v10,v11,D1,D2,(x-floorx),c3);

                /*calculate spline value for line 13 14 15 16 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v13,v14,v15,v16,D1);
                CSPL_INT(f4,f3,f2,f1,v13,v14,v15,v16,D2);
                CSPL_CALC(v14,v15,D1,D2,(x-floorx),c4);

                /*finaly, calculate spline into y direction and save into value[z]*/
                CSPL_INT(f1,f2,f3,f4,c1,c2,c3,c4,D1);
                CSPL_INT(f4,f3,f2,f1,c1,c2,c3,c4,D2);
                CSPL_CALC(c2,c3,D1,D2,(y-floory),value[runz]);
            }
            /*Interpolate on the 4 results into z direction*/
            CSPL_INT(f1,f2,f3,f4,value[0],value[1],value[2],value[3],D1);
            CSPL_INT(f4,f3,f2,f1,value[0],value[1],value[2],value[3],D2);

            return  value[1]+D1*(z-floorz)+(3*(value[2]-value[1])-2*D1-D2)*(z-floorz)*(z-floorz)+(2*(value[1]-value[2])+D1+D2)*(z-floorz)*(z-floorz)*(z-floorz);
    }
    return 0.;
}
