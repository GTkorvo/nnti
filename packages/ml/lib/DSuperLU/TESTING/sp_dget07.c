

/*
 * -- SuperLU routine (version 1.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <math.h>
#include "dsp_defs.h"
#include "util.h"

int sp_dget07(char *trans, int n, int nrhs, SuperMatrix *A, double *b, 
	      int ldb, double *x, int ldx, double *xact, 
              int ldxact, double *ferr, double *berr, double *reslts)
{
/*
    Purpose   
    =======   

    SP_DGET07 tests the error bounds from iterative refinement for the   
    computed solution to a system of equations op(A)*X = B, where A is a 
    general n by n matrix and op(A) = A or A**T, depending on TRANS.
    
    RESLTS(1) = test of the error bound   
              = norm(X - XACT) / ( norm(X) * FERR )   
    A large value is returned if this ratio is not less than one.   

    RESLTS(2) = residual from the iterative refinement routine   
              = the maximum of BERR / ( (n+1)*EPS + (*) ), where   
                (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i ) 

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations.   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

    N       (input) INT
            The number of rows of the matrices X and XACT.  N >= 0.   

    NRHS    (input) INT   
            The number of columns of the matrices X and XACT.  NRHS >= 0. 
  

    A       (input) SuperMatrix *, dimension (A->nrow, A->ncol)
            The original n by n matrix A.   

    B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side vectors for the system of linear   
            equations.   

    LDB     (input) INT   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)   
            The computed solution vectors.  Each vector is stored as a   
            column of the matrix X.   

    LDX     (input) INT   
            The leading dimension of the array X.  LDX >= max(1,N).   

    XACT    (input) DOUBLE PRECISION array, dimension (LDX,NRHS)   
            The exact solution vectors.  Each vector is stored as a   
            column of the matrix XACT.   

    LDXACT  (input) INT   
            The leading dimension of the array XACT.  LDXACT >= max(1,N). 
  

    FERR    (input) DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bounds for each solution vector   
            X.  If XTRUE is the true solution, FERR bounds the magnitude 
            of the largest entry in (X - XTRUE) divided by the magnitude 
            of the largest entry in X.   

    BERR    (input) DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector (i.e., the smallest relative change in any entry of A 
  
            or B that makes X an exact solution).   

    RESLTS  (output) DOUBLE PRECISION array, dimension (2)   
            The maximum over the NRHS solution vectors of the ratios:   
            RESLTS(1) = norm(X - XACT) / ( norm(X) * FERR )   
            RESLTS(2) = BERR / ( (n+1)*EPS + (*) )   

    ===================================================================== 
*/
    
    /* Table of constant values */
    int c__1 = 1;

    /* System generated locals */
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    double diff, axbi;
    int    imax, irow, n__1;
    int    i, j, k;
    double unfl, ovfl;
    double xnorm;
    double errbnd;
    int    notran;
    double eps, tmp;
    double *rwork;
    double *Aval;
    NCformat *Astore;

    /* Function prototypes */
    extern int    lsame_(char *, char *);
    extern int    idamax_(int *, double *, int *);
    extern double dlamch_(char *);

    /* Quick exit if N = 0 or NRHS = 0. */
    if ( n <= 0 || nrhs <= 0 ) {
	reslts[0] = 0.;
	reslts[1] = 0.;
	return 0;
    }

    eps = dlamch_("Epsilon");
    unfl = dlamch_("Safe minimum");
    ovfl   = 1. / unfl;
    notran = lsame_(trans, "N");

    rwork  = (double *) SUPERLU_MALLOC(n*sizeof(double));
    if ( !rwork ) ABORT("SUPERLU_MALLOC fails for rwork");
    Astore = A->Store;
    Aval   = (double *) Astore->nzval;
    
    /* Test 1:  Compute the maximum of   
       norm(X - XACT) / ( norm(X) * FERR )   
       over all the vectors X and XACT using the infinity-norm. */

    errbnd = 0.;
    for (j = 0; j < nrhs; ++j) {
	n__1 = n;
	imax = idamax_(&n__1, &x[j*ldx], &c__1);
	d__1 = fabs(x[imax-1 + j*ldx]);
	xnorm = MAX(d__1,unfl);
	diff = 0.;
	for (i = 0; i < n; ++i) {
	    d__1 = fabs(x[i+j*ldx] - xact[i+j*ldxact]);
	    diff = MAX(diff, d__1);
	}

	if (xnorm > 1.) {
	    goto L20;
	} else if (diff <= ovfl * xnorm) {
	    goto L20;
	} else {
	    errbnd = 1. / eps;
	    goto L30;
	}

L20:
#if 0	
	if (diff / xnorm <= ferr[j]) {
	    d__1 = diff / xnorm / ferr[j];
	    errbnd = MAX(errbnd,d__1);
	} else {
	    errbnd = 1. / eps;
	}
#endif
	d__1 = diff / xnorm / ferr[j];
	errbnd = MAX(errbnd,d__1);
	/*printf("Ferr: %f\n", errbnd);*/
L30:
	;
    }
    reslts[0] = errbnd;

    /* Test 2: Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where 
       (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) + abs(b))_i ) */

    for (k = 0; k < nrhs; ++k) {
	for (i = 0; i < n; ++i) 
            rwork[i] = fabs( b[i + k*ldb] );
	if ( notran ) {
	    for (j = 0; j < n; ++j) {
		tmp = fabs( x[j + k*ldx] );
		for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
		    rwork[Astore->rowind[i]] += fabs(Aval[i]) * tmp;
                }
	    }
	} else {
	    for (j = 0; j < n; ++j) {
		tmp = 0.;
		for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
		    irow = Astore->rowind[i];
		    d__1 = fabs( x[irow + k*ldx] );
		    tmp += fabs(Aval[i]) * d__1;
		}
		rwork[j] += tmp;
	    }
	}

	axbi = rwork[0];
	for (i = 1; i < n; ++i) axbi = MIN(axbi, rwork[i]);
	
	/* Computing MAX */
	d__1 = axbi, d__2 = (n + 1) * unfl;
	tmp = berr[k] / ((n + 1) * eps + (n + 1) * unfl / MAX(d__1,d__2));
	
	if (k == 0) {
	    reslts[1] = tmp;
	} else {
	    reslts[1] = MAX(reslts[1],tmp);
	}
    }

    SUPERLU_FREE(rwork);
    return 0;

} /* sp_dget07 */
