

/* switch betwenn mex and C */
#define MATLAB  
/*#define ANSI_C*/

#define	 PI ((double)3.14159265358979323846264338327950288419716939937510)

/* Input Arguments */
#ifdef MATLAB  
	#define    THETA   prhs[0]
	#define    ACT_TH  prhs[1]
	#define    PSI     prhs[2]	
	#define    D       prhs[3]
	#define    XST     prhs[4]
	#define    YST     prhs[5]
	#define    ZST     prhs[6] 
	#define    W	   prhs[7]
#endif

void dist (float *,int,float *,float *,float *,float *,float *,float *,int ,int,float *);  

