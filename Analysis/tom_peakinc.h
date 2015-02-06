

/* switch betwenn mex and C */
#define MATLAB  
/*#define ANSI_C*/

#define	 PI ((double)3.14159265358979323846264338327950288419716939937510)

/* Input Arguments */
#ifdef MATLAB  
	#define    STREAM    prhs[0]
    #define    DIMS    prhs[1]
    #define     MAX  prhs[2]
	#define     POS_MAX    prhs[3]
#endif


void find_max (float *,int *,float *,int *);