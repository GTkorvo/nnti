

/*
 * -- SuperLU routine (version 1.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */

#include "dsp_defs.h"
#include "util.h"

/*
 * Convert a full matrix into a sparse matrix format. 
 */
int
sp_dconvert(int m, int n, double *A, int lda, int kl, int ku,
	   double *a, int *asub, int *xa, int *nnz)
{
    int     lasta = 0;
    int     i, j, ilow, ihigh;
    int     *row;
    double  *val;

    for (j = 0; j < n; ++j) {
	xa[j] = lasta;
	val = &a[xa[j]];
	row = &asub[xa[j]];

	ilow = MAX(0, j - ku);
	ihigh = MIN(n-1, j + kl);
	for (i = ilow; i <= ihigh; ++i) {
	    val[i-ilow] = A[i + j*lda];
	    row[i-ilow] = i;
	}
	lasta += ihigh - ilow + 1;
    }

    xa[n] = *nnz = lasta;
    return 0;
}


