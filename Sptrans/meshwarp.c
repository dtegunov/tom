/* ======================================================================
 * meshwarp.c -  Mesh warping program.
 * Copyright (C) 1993 by George Wolberg
 *
 * Written by: George Wolberg, 1993
 * ======================================================================
 */

#include "mex.h"
#include "meshwarp.h"


 




/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * saveImage:
 *
 * Save image/mesh I into file.
 */
int saveImage(I, file, type)
imageP	 I;
char	*file;
int	 type;

{
	int	 sz[2];
	FILE	*fp;

	if((fp=fopen(file, "w")) == NULL) {
		fprintf(stderr, "saveImage: Can't create file %s\n", file);
		return (NULL);
	}

	/* write dimensions to file */
	sz[0] = I->width;
	sz[1] = I->height;
	fwrite(sz, sizeof(int), 2, fp);



	
	/* write data to file */
	switch(type) {
	case BW:
		fwrite(I->ch[0], sz[0]*sz[1], 1, fp);
		break;
	case MESH:
		fwrite(I->ch[0], sz[0]*sz[1], 2*sizeof(float), fp);
		break;
	default:
		printf("saveImage: Bad type %d\n", type);
		return(NULL);
	}
	fclose(fp);
}



 


 /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * allocImage:
 *
 * Allocate space for a multi-channel image of width w and height h.
 * If type = BW  : image has 1 channel  of type unsigned char.
 * If type = MESH: image has 2 channels of type float.
 * Return image structure pointer. 
 */
imageP
allocImage(w, h, type)
int w, h, type;
{
	
	imageP	 I;

	/* allocate memory for image data structure */
	I = (imageP) malloc(sizeof(imageS));
	if(I == NULL) {
		fprintf(stderr, "allocImage: Insufficient memory\n");
		return ((imageP) NULL);
	}

	/* init dimensions */
	I->width  = w;
	I->height = h;

	/* init channel pointers */
	switch(type) {
	case BW:
		I->ch[0] = (float *) malloc(w*h*sizeof(float));
		
		break;
	case MESH:
		I->ch[0] = (float *) malloc(2*w*h*sizeof(float));
		I->ch[1] = (float *) I->ch[0] + w*h;
		
		break;
	default:
		printf("allocImage: Bad type %d\n", type);
		return((imageP) NULL);
	}

	return(I);
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * freeImage:
 *
 * Free image/mesh memory.
 */
void freeImage(imageP I)

{
	free((float *) I->ch[0]);
	free((float *) I);
}

 
 

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * catmullRom:
 *
 * Compute a Catmull-Rom spline passing through the len1 points in arrays
 * x1, y1, where y1 = f(x1)
 * len2 positions on the spline are to be computed. Their positions are
 * given in x2. The spline values are stored in y2.
 */
void
catmullRom(x1, y1, len1, x2, y2, len2)
float	*x1, *y1, *x2, *y2;
int	 len1, len2;
{
	int i, j, dir, j1, j2;
	double x,  dx1, dx2;
	double dx, dy, yd1, yd2, p1, p2, p3;
	double a0y, a1y, a2y, a3y;
    FILE *fp_debug;

	/* find direction of monotonic x1; skip ends */
	if(x1[0] < x1[1]) 
	{	/* increasing */
		if(x2[0]<x1[0] || x2[len2-1]>x1[len1-1]) dir=0;
		else dir = 1;
	} 
	else 
	{		/* decreasing */
		if(x2[0]>x1[0] || x2[len2-1]<x1[len1-1]) dir=0;
		else dir = -1;
	}
	
	/* db output*/
	/*fp_debug=fopen("check","a");
	fprintf(fp_debug,"checking (x1[0]=%f x1[1]=%f)  (x2[0]=%f  x1[0]=%f) (x2[len2-1]=%f  x1[len1-1]=%f)  dir=%d  \n",
		    x1[0],x1[1],x2[0],x1[0],x2[len2-1],x1[len1-1],dir); 
	fclose(fp_debug);*/
	/* db output*/
	
	if(dir == 0) {			/* error */
		printf("catmullRom: Output x-coord out of range of input\n");
		return;
	}

	/* p1 is first endpoint of interval
	 * p2 is resampling position
	 * p3 is second endpoint of interval
	 * j  is input index for current interval
	 */

	/* force coefficient initialization */
	if(dir==1)	p3 = x2[0] - 1;
	else		p3 = x2[0] + 1;

	for(i=0; i<len2; i++) {
		/* check if in new interval */
		p2 = x2[i];
		if((dir==1 && p2>p3) || (dir== -1 && p2<p3)) {
			/* find the interval which contains p2 */
			if(dir) {
				for(j=0; j<len1 && p2>x1[j]; j++);
				if(p2 < x1[j]) j--;
			} else {
				for(j=0; j<len1 && p2<x1[j]; j++);
				if(p2 > x1[j]) j--;
			}

			p1 = x1[j];		/* update 1st endpt */
			p3 = x1[j+1];		/* update 2nd endpt */

			/* clamp indices for endpoint interpolation */
			j1 = MAX(j-1, 0);
			j2 = MIN(j+2, len1-1);
				
			/* compute spline coefficients */
			dx  = 1.0 / (p3 - p1);
			dx1 = 1.0 / (p3 - x1[j1]);
			dx2 = 1.0 / (x1[j2] - p1);
			dy  = (y1[j+1] - y1[ j ]) * dx;
			yd1 = (y1[j+1] - y1[ j1]) * dx1;
			yd2 = (y1[j2 ] - y1[ j ]) * dx2;
			a0y =  y1[j];
			a1y =  yd1;
			a2y =  dx *  ( 3*dy - 2*yd1 - yd2);
			a3y =  dx*dx*(-2*dy +   yd1 + yd2);
		}
		/* use Horner's rule to calculate cubic polynomial */
		x = p2 - p1;
		y2[i] = ((a3y*x + a2y)*x + a1y)*x + a0y;
	}
}
 



 /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * meshWarp:
 *
 * Warp I1 with correspondence points given in meshes M1 and M2.
 * Result goes in I2.
 *
 * See "Digital Image Warping" by George Wolberg (IEEE Computer Society
 * Press, 1990) for details.
 * Based on Douglas Smythe's algorithm (in "A Two-Pass Mesh Warping Algorithm
 * for Object Transformation and Image Interpolation", ILM Technical Memo
 * #1030, 1990).
 */
void
meshWarp(I1, M1, M2, I2)
imageP	I1, I2;
imageP	M1, M2;
{
	int	 I_w, I_h, M_w, M_h;
	int	 x, y, u, v, n;
	float	*x1, *y1, *x2, *y2;
	float	*xrow, *yrow, *xcol, *ycol, *coll, *indx, *map;
	float	*src, *dst;
	imageP	 Mx, My, I3;
	FILE *fp_debug;
	
	I_w = I1->width;
	I_h = I1->height;

	M_w = M1->width;
	M_h = M1->height;


	
	/* allocate enough memory for a scanline along the longest dimension */
	n = MAX(I_w, I_h);
	indx = (float *) malloc(n * sizeof(float));
	xrow = (float *) malloc(n * sizeof(float));
	yrow = (float *) malloc(n * sizeof(float));
	map  = (float *) malloc(n * sizeof(float));

	/* db output*/
	/*fp_debug=fopen("check","a");
	fprintf(fp_debug," create table of x-intercepts for source mesh's vertical splines  \n "); 
	fclose(fp_debug);*/
	/* db output*/
	
	
	/* create table of x-intercepts for source mesh's vertical splines */
	Mx = allocImage(M_w, I_h, MESH);
	for(y=0; y < I_h; y++) indx[y] = y;
	for(u=0; u < M_w; u++) {	/* visit each vertical spline   */
		/* store column as row for spline fct */
		xcol = (float *) M1->ch[0] + u;
		ycol = (float *) M1->ch[1] + u;
		coll = (float *) Mx->ch[0] + u;


		/* scan convert vertical splines */
		for(v=0; v < M_h; v++, xcol+=M_w) 
		{
		     xrow[v] = *xcol;
		     
		
		}
		for(v=0; v < M_h; v++, ycol+=M_w) 
		{
			yrow[v] = *ycol;
		    
		}
		
		catmullRom(yrow, xrow, M_h, indx, map, I_h);

		/* store resampled row back into column */
		for(y=0; y < I_h; y++, coll+=M_w) *coll = map[y];
	}

	
	
	/* db output*/
	/*fp_debug=fopen("check","a");
	fprintf(fp_debug,"  create table of x-intercepts for dst mesh's vertical splines \n "); 
	fclose(fp_debug);*/
	/* db output*/
	
	
	/* create table of x-intercepts for dst mesh's vertical splines */
	for(u=0; u < M_w; u++) {	/* visit each  vertical spline  */
		/* store column as row for spline fct */
		xcol = (float *) M2->ch[0] + u;
		ycol = (float *) M2->ch[1] + u;
		coll = (float *) Mx->ch[1] + u;

		/* scan convert vertical splines */
		for(v=0; v < M_h; v++, xcol+=M_w) xrow[v] = *xcol;
		for(v=0; v < M_h; v++, ycol+=M_w) yrow[v] = *ycol;
		catmullRom(yrow, xrow, M_h, indx, map, I_h);

		/* store resampled row back into column */
		for(y=0; y < I_h; y++, coll+=M_w) *coll = map[y];
	}

	
	
		/* db output*/
	/*fp_debug=fopen("check","a");
	fprintf(fp_debug,"  first pass: warp x using tables in Mx \n "); 
	fclose(fp_debug);*/
	/* db output*/
	
	/* first pass: warp x using tables in Mx */
	I3  = allocImage(I_w, I_h, BW);
	x1  = (float *) Mx->ch[0];
	x2  = (float *) Mx->ch[1];
	src = (float *) I1->ch[0];
	dst = (float *) I3->ch[0];
	for(x=0; x < I_w; x++) indx[x] = x;
	for(y=0; y < I_h; y++) {
		/* fit spline to x-intercepts; resample over all cols */
		
	/* db output*/
	/*fp_debug=fopen("check","a");
	fprintf(fp_debug,"outside  x1[0] %f \n ",x1[0]); 
	fclose(fp_debug);*/
	/* db output*/
		
		
		
		catmullRom(x1, x2, M_w, indx, map, I_w);

		/* resample source row based on map */
		resample(src, I_w, 1, map, dst);

		/* advance pointers to next row */
		src += I_w;
		dst += I_w;
		x1  += M_w;
		x2  += M_w;
	}
	freeImage(Mx);

	/* db output*/
	/*fp_debug=fopen("check","a");
	fprintf(fp_debug,"	 create table of y-intercepts for intermediate mesh's hor splines  \n "); 
	fclose(fp_debug);*/
	/* db output*/


    /* create table of y-intercepts for intermediate mesh's hor splines */
	My = allocImage(I_w, M_h, MESH);
	x1 = (float *) M2->ch[0];
	y1 = (float *) M1->ch[1];
	y2 = (float *) My->ch[0];
	for(x=0; x < I_w; x++) indx[x] = x;
	for(v=0; v < M_h; v++) {	/* visit each horizontal spline */
		/* scan convert horizontal splines */
		catmullRom(x1, y1, M_w, indx, y2, I_w);

		/* advance pointers to next row */
		x1 += M_w;
		y1 += M_w;
		y2 += I_w;
	}

	/* create table of y-intercepts for dst mesh's horizontal splines */
	x1 = (float *) M2->ch[0];
	y1 = (float *) M2->ch[1];
	y2 = (float *) My->ch[1];
	for(v=0; v < M_h; v++) {	/* visit each horizontal spline   */
		/* scan convert horizontal splines */
		catmullRom(x1, y1, M_w, indx, y2, I_w);

		/* advance pointers to next row */
		x1 += M_w;
		y1 += M_w;
		y2 += I_w;
	}

	/* second pass: warp y */
	src = (float *) I3->ch[0];
	dst = (float *) I2->ch[0];
	for(y=0; y < I_h; y++) indx[y] = y;
	for(x=0; x < I_w; x++) {
		/* store column as row for spline fct */
		xcol = (float *) My->ch[0] + x;
		ycol = (float *) My->ch[1] + x;
		for(v=0; v < M_h; v++, xcol+=I_w) xrow[v] = *xcol;
		for(v=0; v < M_h; v++, ycol+=I_w) yrow[v] = *ycol;

		/* fit spline to y-intercepts; resample over all rows */
		catmullRom(xrow, yrow, M_h, indx, map, I_h);

		/* resample source column based on map */
		resample(src, I_h, I_w, map, dst);

		/* advance pointers to next column */
		src++;
		dst++;
	}
	freeImage(My);
	freeImage(I3);
	free((float *) indx);
	free((float *) xrow);
	free((float *) yrow);
	free((float *)  map);
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * resample:
 *
 * Resample the len elements of src (with stride offst) into dst according
 * to the spatial mapping given in xmap.
 * Perform linear interpolation for magnification and box filtering
 * (unweighted averaging) for minification.
 * Based on Fant's algorithm (IEEE Computer Graphics & Applications, 1/86).
 */
void
resample(src, len, offst, xmap, dst)
float	*src, *dst;
float	*xmap;
int	 len, offst;
{
	int u, x, v0, v1;
	double val, sizfac, inseg, outseg, acc, inpos[1024];

	/* precompute input index for each output pixel */
	for(u=x=0; x<len; x++) {
		while(xmap[u+1]<x) u++;
		inpos[x] = u + (double) (x-xmap[u]) / (xmap[u+1]-xmap[u]);
	}

	inseg  = 1.0;
	outseg = inpos[1];
	sizfac = outseg;
	acc = 0.;
	v0 = *src;	src += offst;
	v1 = *src;	src += offst;
	for(u=1; u<len; ) {
		val = inseg*v0 + (1-inseg)*v1;
		if(inseg < outseg) {
			acc += (val * inseg);
			outseg -= inseg;
			inseg = 1.0;
			v0 = v1;
			v1 = *src;
			src += offst;
		} else {
			acc += (val * outseg);
			acc /= sizfac;
			*dst = (int) MIN(acc, 0xff);
			dst += offst;
			acc = 0.;
			inseg -= outseg;
			outseg = inpos[u+1] - inpos[u];
			sizfac = outseg;
			u++;
		}
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
const int *dims_im1, *dims_im2, *dims_mash1, *dims_mash2;
int w,h,lauf;
float *src;
FILE *fp_debug;

imageP	I1, I2;
imageP	M1, M2;

/* get dimensions */
dims_im1=mxGetDimensions(prhs[0]);
dims_mash1=mxGetDimensions(prhs[1]);
dims_mash2=mxGetDimensions(prhs[2]);
dims_im2=mxGetDimensions(prhs[3]);

/*
printf(" %d %d %d %d",dims_mash1[0],dims_mash1[1],dims_mash2[0],dims_mash2[1]);
*/


/* allocate memory and init structure for output image */
I1 = allocImage(dims_im1[0],dims_im1[1], BW);
I2 = allocImage(dims_im2[0],dims_im2[1], BW);


M1 = allocImage(dims_mash1[0],dims_mash1[1], MESH);
M2 = allocImage(dims_mash2[0],dims_mash2[1], MESH);



/*transfer data*/
I1->ch[0]=  mxGetData(prhs[0]);
M1->ch[0]=  mxGetData(prhs[1]); 
M2->ch[0]=  mxGetData(prhs[2]);
w=M1->width;
h=M1->height;
M1->ch[1] = (float *) M1->ch[0] + w*h;
w=M2->width;
h=M2->height;
M2->ch[1] = (float *) M2->ch[0] + w*h;
I2->ch[0]=  mxGetData(prhs[3]);



/*db output*/
src=M2->ch[0];

/*fp_debug=fopen("check","a");
fprintf(fp_debug,"pre start %f \n ",src[0]); 
fclose(fp_debug);*/
/* db output*/


/* do actual computation */
meshWarp(I1, M1, M2, I2); 


freeImage(I1);
freeImage(I2);
freeImage(M1);    
freeImage(M2);
/*saveImage(I2,"test", BW); */


}  
