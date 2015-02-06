

/* switch betwenn mex and C */
#define MATLAB  
/*#define ANSI_C*/

#define	 PI ((double)3.14159265358979323846264338327950288419716939937510)

/* Input Arguments */
#ifdef MATLAB  
	#define    DIM_WFUNK  prhs[0]
    #define    ALL_ANGLES  prhs[1]
	#define	   THICKNESS  prhs[2]
    #define    ANGLE_PROJ prhs[3]
    #define    OUT plhs[0]   
     
#endif


/* void rot2d (float *,float *,long,long,*float,char,float,float,int); */
