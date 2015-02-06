#include "mex.h"
#include "math.h"
/* undef needed for LCC compiler */
#undef EXTERN_C
#ifdef _WIN32
	#include <windows.h>
	#include <process.h>
#else
	#include <pthread.h>
#endif

/*
 Affine transformation function (Rotation, Translation, Resize)
 This function transforms a volume with a 3x3 transformation matrix 

 Iout=affine_transform_2d_double(Iin,Minv,ImageSize,check_bil_intp)

 inputs,
   Iin: The greyscale input image
   Minv: The (inverse) 3x3 transformation matrix
   ImageSize: Size of output imgae
   check_bil_intp: If 1 bilinear interpolation otherwise nearest neigb.
 output,
   Iout: The transformed image

 example,
   % Read image
   I=im2double(imread('lenag2.png'))
   % Make a transformation matrix
   M=make_transformation_matrix([2 3],[1.0 1.1],2);
   % Transform the image
   Iout=rigid_transform_2d_double(I,M,size(I))
   % Show the image
   figure, imshow(Iout);

% Function is written by D.Kroon University of Twente (February 2009)
*/

/* Convert 2D/3D matrix index to 1D index */
int mindex2(int x, int y, int sizx, int sizy) 
{ 
    if(x<0) { x=0; }
    if(y<0) { y=0; }
    if(x>(sizx-1)) { x=sizx-1; }
    if(y>(sizy-1)) { y=sizy-1; }
    
    return y*sizx+x; 
}

/* Convert 2D/3D matrix index to 1D index */
int mindex2c(int x, int y, int sizx, int sizy, double *perc, int index_perc) 
{ 
    if(x<0) { perc[index_perc]=0; return 0;}
    if(y<0) { perc[index_perc]=0; return 0;}
    if(x>(sizx-1)) { perc[index_perc]=0; return 0;}
    if(y>(sizy-1)) { perc[index_perc]=0; return 0;}
    
    return y*sizx+x; 
}


int mindex3(int x, int y, int z, int sizx, int sizy,int sizz) 
{ 
    if(x<0) { x=0; }
    if(y<0) { y=0; }
    if(z<0) { z=0; }
    if(x>(sizx-1)) { x=sizx-1; }
    if(y>(sizy-1)) { y=sizy-1; }
    if(z>(sizz-1)) { z=sizz-1; }
    return z*sizx*sizy+y*sizx+x;
}


#ifdef _WIN32
  unsigned __stdcall transformvolume(double **Args) {
#else
  void transformvolume(double **Args) {
#endif
    double *Isize_d, *mean_in, *A, *Iin, *Iout, *ThreadID;
	double *ImageSize_d, *check_bil_intp;
	double mean_out[2]={0,0};
	int bilintp=0;
    int Isize[2]={0,0};
	int ImageSize[2]={0,0};
	
    int x,y;
	double *Nthreadsd;
    int Nthreads;

    /* Location of pixel which will be come the current pixel */
    double Tlocalx;
    double Tlocaly;
    
    /* Linear interpolation variables */
    int xBas[4], yBas[4];
    double perc[4]={0,0,0,0};
    double xCom, yCom;
    double color[4]={0,0,0,0};
    
    /* Color offsets */
    int offset_out[3]={0,0,0};
    int offset_in[3]={0,0,0};
    
    /* X,Y,Z coordinates of current pixel */
    double xd,yd;

    /* Variables to store 1D index */
    int indexI;
    int indexI1, indexI2, indexI3, indexI4;
    
    /* Multiple threads, one does the odd the other even indexes */
    int offset;
    
    /* Loop through all colors r,g,b or only gray */
    int c;
    
    
    Isize_d=Args[0];
    mean_in=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iout=Args[4];
    ThreadID=Args[5];
	ImageSize_d=Args[6];
	check_bil_intp=Args[7]; bilintp=(int)check_bil_intp[0];
	Nthreadsd=Args[8];  Nthreads=(int)Nthreadsd[0];
   
      /* Center of the output image */
	mean_out[0]=ImageSize_d[0]/2;  mean_out[1]=ImageSize_d[1]/2;  
  
    ImageSize[0] = (int)ImageSize_d[0]; 
    ImageSize[1] = (int)ImageSize_d[1]; 
    
	Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    Isize[2] = (int)Isize_d[2]; 
    
    offset=(int) ThreadID[0];
	
    offset_out[0]=0;
    offset_in[0]=0;
    offset_out[1]=ImageSize[0]*ImageSize[1];
    offset_in[1]=Isize[0]*Isize[1];
    offset_out[2]=2*ImageSize[0]*ImageSize[1];
    offset_in[2]=2*Isize[0]*Isize[1];
             
    /* Loop through all image pixel coordinates */
    for (y=offset; y<ImageSize[1]; y=y+Nthreads)
    {
        for (x=0; x<ImageSize[0]; x++)
        {
            xd=(double)x-mean_out[0]; yd=(double)y-mean_out[1];

            Tlocalx = mean_in[0] + A[0] * xd + A[1] *yd + A[2] * 1;
            Tlocaly = mean_in[1] + A[3] * xd + A[4] *yd + A[5] * 1;

			if(bilintp>0)
			{
				/* Determine the coordinates of the pixel(s) which will be come the current pixel */
				/* (using linear interpolation)  */
				xBas[0]=(int) floor(Tlocalx); yBas[0]=(int) floor(Tlocaly);
				xBas[1]=xBas[0]+0;      yBas[1]=yBas[0]+1;
				xBas[2]=xBas[0]+1;      yBas[2]=yBas[0]+0;
				xBas[3]=xBas[0]+1;      yBas[3]=yBas[0]+1;

				
                /* Linear interpolation constants (percentages) */
				xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);
				perc[0]=(1-xCom) * (1-yCom);
				perc[1]=(1-xCom) * yCom;
				perc[2]=xCom * (1-yCom);
				perc[3]=xCom * yCom;

                indexI1=mindex2c(xBas[0],yBas[0],Isize[0],Isize[1],perc,0);
                indexI2=mindex2c(xBas[1],yBas[1],Isize[0],Isize[1],perc,1);
                indexI3=mindex2c(xBas[2],yBas[2],Isize[0],Isize[1],perc,2);
                indexI4=mindex2c(xBas[3],yBas[3],Isize[0],Isize[1],perc,3);
    
                /* Get index current pixel value */
				indexI=mindex2(x,y,ImageSize[0],ImageSize[1]);

                for(c=0; c<Isize[2]; c++)
                {
                    color[0]=Iin[indexI1+offset_in[c]]; 
                    color[1]=Iin[indexI2+offset_in[c]];
                    color[2]=Iin[indexI3+offset_in[c]];
                    color[3]=Iin[indexI4+offset_in[c]];
                    Iout[indexI+offset_out[c]]=color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
                }
			}
			else
            {  
                /* Get index current pixel value */
				indexI=mindex2(x,y,ImageSize[0],ImageSize[1]);
                perc[0]=1;
                indexI1=mindex2c((int) floor(Tlocalx+0.5),(int) floor(Tlocaly+0.5),Isize[0],Isize[1],perc,0);
                for(c=0; c<Isize[2]; c++)
                {
                    Iout[indexI+offset_out[c]]=Iin[indexI1+offset_in[c]]*perc[0];
                }
			}
        }
    }
   
    /* explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    #ifdef _WIN32
	_endthreadex( 0 );
    return 0;
	#else
	pthread_exit(NULL);
	#endif
    
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Ox and Oy are the grid points
       Zo is the input image
       Zi is the transformed image
       nx and ny are the number of grid points (inside the image) */
    double *Iin, *Iout, *M, *ImageSize, *check_bil_intp; 
    mxArray *matlabCallOut[1]={0};
    mxArray *matlabCallIn[1]={0};
    double *Nthreadsd;
    int Nthreads;
	
    /* double pointer array to store all needed function variables */
    double ***ThreadArgs;
    double **ThreadArgs1;
    
    /* Handles to the worker threads */
	#ifdef _WIN32
		HANDLE *ThreadList; 
    #else
		pthread_t *ThreadList;
	#endif
    
    /* ID of Threads */
    double **ThreadID;              
    double *ThreadID1;
    
    /* Transformation matrix */
    double A[9]={0,0,0,0,0,0,0,0,0};

    /* Loop variable */
    int i;
    
    /* Size of input image */
    double Isize_d[3]={0,0,0};
    const mwSize *dims;
	mwSize dims_out3[3]={0,0,0};
    
    mwSize dimnum;
    
    double mean_in[2]={0,0};

  /* Check for proper number of arguments. */
  if(nrhs!=4) {
    mexErrMsgTxt("Four inputs are required.");
  } else if(nlhs!=1) {
    mexErrMsgTxt("One output required");
  }
     
  /* Get the sizes of the image */
  dimnum=mxGetNumberOfDimensions(prhs[0]);
    
  dims = mxGetDimensions(prhs[0]);   
  Isize_d[0] = (double)dims[0]; Isize_d[1] = (double)dims[1]; 
  if(dimnum>2) 
  { 
      Isize_d[2] = (double)dims[2];
  } 
  else 
  { 
      Isize_d[2]=1; 
  }
          
  /* Assign pointers to each input. */
  Iin=mxGetPr(prhs[0]);
  M=mxGetPr(prhs[1]);
  ImageSize=mxGetPr(prhs[2]);
  check_bil_intp=mxGetPr(prhs[3]);

  dims_out3[0]=(mwSize)ImageSize[0]; dims_out3[1]=(mwSize)ImageSize[1]; dims_out3[2]=(mwSize)Isize_d[2];
  plhs[0] = mxCreateNumericArray(3, dims_out3, mxDOUBLE_CLASS, mxREAL);
 
  A[0] = M[mindex2(0,0,3,3)]; A[1] = M[mindex2(0,1,3,3)]; A[2] = M[mindex2(0,2,3,3)]; 
  A[3] = M[mindex2(1,0,3,3)]; A[4] = M[mindex2(1,1,3,3)]; A[5] = M[mindex2(1,2,3,3)]; 
  A[6] = M[mindex2(2,0,3,3)]; A[7] = M[mindex2(2,1,3,3)]; A[8] = M[mindex2(2,2,3,3)]; 
  
	mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
	Nthreadsd=mxGetPr(matlabCallOut[0]);
	Nthreads=(int)Nthreadsd[0];
	/* Reserve room for handles of threads in ThreadList */
    #ifdef _WIN32
		ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
    #else
		ThreadList = (pthread_t*)malloc(Nthreads* sizeof( pthread_t ));
	#endif
	
    
	ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
	ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );
	
  
  /* Assign pointer to output. */
  Iout = mxGetPr(plhs[0]);
  
  /* Center of the volume */
  mean_in[0]=Isize_d[0]/2;  mean_in[1]=Isize_d[1]/2;  
   
  for (i=0; i<Nthreads; i++)
  {
    /* Make Thread ID */
    ThreadID1= (double *)malloc( 1* sizeof(double) );
    ThreadID1[0]=i;
    ThreadID[i]=ThreadID1;  
    
	/* Make Thread Structure */
    ThreadArgs1 = (double **)malloc( 8* sizeof( double * ) );  
	ThreadArgs1[0]=Isize_d;
	ThreadArgs1[1]=mean_in;
	ThreadArgs1[2]=A;
	ThreadArgs1[3]=Iin;
	ThreadArgs1[4]=Iout;
	ThreadArgs1[5]=ThreadID[i];
	ThreadArgs1[6]=ImageSize;
	ThreadArgs1[7]=check_bil_intp;
	ThreadArgs1[8]=Nthreadsd;
    /* Start a Thread  */
	ThreadArgs[i]=ThreadArgs1;
    	

    #ifdef _WIN32
        ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs[i] , 0, NULL );
    #else
        pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume, ThreadArgs[i]);
    #endif
  }

  #ifdef _WIN32
	for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
 	for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
  #else
	for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
  #endif

  for (i=0; i<Nthreads; i++) 
  { 
    free(ThreadArgs[i]);
    free(ThreadID[i]);
  }

  free(ThreadArgs);
  free(ThreadID );
  free(ThreadList);

}
        

