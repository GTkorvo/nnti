#ifndef __MLGGB_
#define __MLGGB_



#define SHIFTS   0 
#define MAX_ITRS 2 
#define MODE     6 

#include <stdio.h>
#include <stdlib.h>
#include "ml_common.h"
#include "ml_mat_formats.h" 
#include "ml_lapack.h"
#include "ml_eigf2c.h"


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

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif



void  ML_ARPACK_driver(char which[],
			 char bmat[], int iparam[], int mode,
			 int nev, int ncv, double tol,  ML *ml,
		       struct ML_CSR_MSRdata  *mydata, int Fattening,
		       struct ML_Eigenvalue_Struct *eigen_struct,
		       int Debug_Flag, int GGB_alp_flag);


  void ML_GGB2CSR (double *v, int nconv, int MatSize, int proc_id,
		   struct ML_CSR_MSRdata  *mydata, int Debug_Flag );


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



void dnaupd_(int *, char *, int *, char *, int *, double *, double *,
		 int *, double *, int *, int *, int *, double *, double *,
		 int *, int *);

  void pdnaupd_(int *, int *, char *, int *, char *, int *, double *, double *,
		int *, double *, int *, int *, int *, double *, double *,
		int *, int *);

 
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
