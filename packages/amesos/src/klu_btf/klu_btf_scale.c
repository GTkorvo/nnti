/* ========================================================================== */
/* === klu_btf_scale ======================================================== */
/* ========================================================================== */

/* Scale a matrix.  Can be called by the user, but not needed.
 * This is called by klu_btf_factor and klu_btf_refactor.
 * Returns KLU_OK if the input matrix is valid, KLU_INVALID otherwise.
 * If the W input argument is non-NULL, then the input matrix is checked for
 * duplicate entries.
 */

#include "klu_btf_internal.h"

int klu_btf_scale
(
    /* inputs, not modified */
    int n,
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],
    /* workspace, not defined on input or output */
    int W [ ]		/* size n, can be NULL */
)
{
    int row, col, p, pend, check_duplicates ;

    if (n <= 0 || !Ap || !Ai || !Ax)
    {
	/* Ap, Ai, and Ax must be present, and n must be > 0 */
	return (KLU_INVALID) ;
    }
    if (Ap [0] != 0 || Ap [n] < 0)
    {
	/* nz = Ap [n] must be >= 0 and Ap [0] must equal zero */
	return (KLU_INVALID) ;
    }
    for (col = 0 ; col < n ; col++)
    {
	if (Ap [col] > Ap [col+1])
	{
	    /* column pointers must be non-decreasing */
	    return (KLU_INVALID) ;
	}
    }

    /* initialize row sum and workspace */
    for (row = 0 ; row < n ; row++)
    {
	Rs [row] = 0.0 ;
    }

    /* check for duplicates only if W is present */
    check_duplicates = (W != (int *) NULL) ;
    if (check_duplicates)
    {
	for (row = 0 ; row < n ; row++)
	{
	    W [row] = EMPTY ;
	}
    }

    /* the scale factor for row i is sum (abs (A (i,:))) */
    for (col = 0 ; col < n ; col++)
    {
	pend = Ap [col+1] ;
	for (p = Ap [col] ; p < pend ; p++)
	{
	    row = Ai [p] ;
	    if (row < 0 || row >= n)
	    {
		/* row index out of range, or duplicate entry */
		return (KLU_INVALID) ;
	    }
	    if (check_duplicates)
	    {
		if (W [row] == col)
		{
		    /* duplicate entry */
		    return (KLU_INVALID) ;
		}
		/* flag row i as appearing in column col */
		W [row] = col ;
	    }
	    /* accumulate the row sum */
	    Rs [row] += ABS (Ax [p]) ;
	}
    }

    /* do not scale empty rows */
    for (row = 0 ; row < n ; row++)
    {
	/* matrix is singular */
	PRINTF (("Rs [%d] = %g\n", row, Rs [row])) ;
	if (Rs [row] == 0.0)
	{
	    PRINTF (("Row %d of A is all zero\n", row)) ;
	    Rs [row] = 1.0 ;
	}
    }
    return (KLU_OK) ;
}
