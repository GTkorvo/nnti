
/*************************************************************************************       
  HAIM: GLOBAL EIGENVALUE CALCULATIONS FOR MULTIGRID           
        GENERALIZED GLOBAL BASIS (GGB) METHOD USING ARPACK                      
**************************************************************************************/     

#define SHIFTS   0 
#define MAX_ITRS 2 
#define MODE     6 


#include <stdio.h>
#include <stdlib.h>
//#include "ml_struct.h" 
#include "ml_mat_formats.h" 
struct ML_Eigenvalue_Struct  {
  int     Max_Iter;                  /* User input from input file */ 
  int     Num_Eigenvalues;
  int     Arnoldi;
  double  Residual_Tol;
  int     Fattening;


  int     Nflag;          /* Flag to indicate the first Newton iteration */
  int     Pnconv;         /* Previous number of converged eigenvalues */
  double *Evec, *Eval;    /* eigenvectors and eigenvalues to be reused 
			     with MGGB */

};


void  ML_ARPACK_driver(char which[],
			 char bmat[], int iparam[], int mode,
			 int nev, int ncv, double tol,  ML *ml,
		       struct ML_CSR_MSRdata  *mydata, int Fattening,
		       struct ML_Eigenvalue_Struct *eigen_struct,
		       int Debug_Flag, int GGB_alp_flag);


void ML_GGB_2_CSR (double **eigvec, int nconv, int MatSize, 
		   struct ML_CSR_MSRdata  *mydata, int Debug_Flag );
void ML_GGB2CSR (double *v, int nconv, int MatSize, 
		   struct ML_CSR_MSRdata  *mydata, int Debug_Flag );


void  ML_mgs (double *NewVec, int nconv, int nloc2, struct ML_Eigenvalue_Struct 
	      *eigen_struct,  int Debug_Flag);
	      
void  ML_GGBalp (double *NewVec, int nconv, int nloc2, struct ML_Eigenvalue_Struct 
		   *eigen_struct);

extern double  ML_subspace (int nrows, double *inp1, int ncols1, double *inp2, int ncols2);



extern void ML_ARPACK_GGB( 
		    struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
		    struct ML_CSR_MSRdata *mydata, int Debug_Flag, 
		    int GGB_alp_flag);

extern int  ML_MGGB_angle(struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
		    struct ML_CSR_MSRdata *mydata);

extern double ML_normc(double *real, double *imag,  int leng );  




