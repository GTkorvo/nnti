#ifndef _IFP_DENSEMAT_H_
#define _IFP_DENSEMAT_H_

#include <math.h>
#include <assert.h>
#include "ifp_LocalMat.h"
#include "ifp_BlockVec.h"
#include "ifp_LocalPrecon.h"
#include "ifp_ifpack.h"
#include "ifp_blas3.h"
#include "ifp_lapackd.h"

// ifp_DenseMat objects typically allocate their own storage for matrix data.
// The default constructor is the only constructor that does not allocate
// storage.  Thus this constructor may be used to control memory
// allocation, for example by allocating an array of ifp_DenseMat using new.
// In this case, the data pointer must be set to NULL before the destructor
// is called.

class ifp_DenseMat : public ifp_LocalMat
{
// Temporary change protected:
public:
    double *a;
    int    nrow;
    int    ncol;

public:
    ifp_DenseMat() {a = NULL; nrow = ncol = 0;}
   ~ifp_DenseMat() {delete [] a;}

    ifp_DenseMat(const int r, const int c)
        {nrow = r; ncol = c; a = new double[r*c];}

    ifp_DenseMat(const ifp_DenseMat& A);
// cannot inline on Cray
//      {nrow = A.nrow; 
//       ncol = A.ncol;
//       register double *p = a = new double[nrow*ncol];
//       register double *q = A.a;
//       for (int i=0; i<nrow*ncol; i++) *p++ = *q++;}

    void set(const int r, const int c, double *d)
        {nrow = r; ncol = c; a = d;}

    int numrow() const {return nrow;}
    int numcol() const {return ncol;}

    double operator=(double s)
        {register double *p = a;
         for (int i=0; i<nrow*ncol; i++) *p++ = s; return s;}

    // 0-based indexing
    const double& operator()(const unsigned int i, const unsigned int j) const
        {return a[j*nrow+i];}
    double& operator()(const unsigned int i, const unsigned int j)
        {return a[j*nrow+i];}

    // virtual functions

    double *& Data() {return a;}
    const double *Data() const {return a;}
    ifp_LocalMat *CreateEmpty() const {return new ifp_DenseMat();}
    inline ifp_LocalMat *CreateInv(ifp_LocalPrecon&) const;
    inline void SetToZero(int, int);
    void MatCopy(const ifp_LocalMat& A)  // UNDONE: realloc if necessary
        {register double *p = a;
         register double *q = ((ifp_DenseMat *)&A)->a;
         for (int i=0; i<nrow*ncol; i++) *p++ = *q++;}
    void Print(ostream&) const;

    // virtual mathematical functions

    inline void Mat_Trans(ifp_LocalMat *B) const;
    inline void Mat_Mat_Add(const ifp_LocalMat *B, ifp_LocalMat *C, double alpha) const;
    inline void Mat_Mat_Mult(const ifp_LocalMat *B, ifp_LocalMat *C, 
        double alpha, double beta) const;
    inline void Mat_Vec_Mult(const ifp_BlockVec& B, ifp_BlockVec& C, 
        double alpha, double beta) const;
    inline void Mat_Trans_Vec_Mult(const ifp_BlockVec& B, ifp_BlockVec& C,
        double alpha, double beta) const;
    inline void Mat_Vec_Solve(const ifp_BlockVec& B, ifp_BlockVec& C) const;
    inline void Mat_Trans_Vec_Solve(const ifp_BlockVec& B, ifp_BlockVec& C) const;
};

ostream& operator << (ostream& os, const ifp_DenseMat& mat);

class ifp_DenseMat_LU : public ifp_DenseMat
{
private:
    double  *lu;
    integer *ipiv;  // FORTRAN integer

public:
    inline ifp_DenseMat_LU(const ifp_DenseMat&);
    inline void Mat_Vec_Solve(const ifp_BlockVec&, ifp_BlockVec&) const;
    inline ~ifp_DenseMat_LU();
};

class ifp_DenseMat_INVERSE : public ifp_DenseMat
{
private:
public:
    inline ifp_DenseMat_INVERSE(const ifp_DenseMat&);
    void Mat_Vec_Solve(const ifp_BlockVec& B, ifp_BlockVec& X) const
        {solve_is_mult(B, X);}
};

class ifp_DenseMat_SVD : public ifp_DenseMat
{
private:
public:
    inline ifp_DenseMat_SVD(const ifp_DenseMat&, double, double);
    void Mat_Vec_Solve(const ifp_BlockVec& B, ifp_BlockVec& X) const
        {solve_is_mult(B, X);}
};

class ifp_DenseMat_DIAGDOM : public ifp_DenseMat
{
private:
public:
    inline ifp_DenseMat_DIAGDOM(const ifp_DenseMat&, double alpha);
    void Mat_Vec_Solve(const ifp_BlockVec& B, ifp_BlockVec& X) const
        {solve_is_mult(B, X);}
};

class ifp_DenseMat_GERSH : public ifp_DenseMat
{
private:
public:
    inline ifp_DenseMat_GERSH(const ifp_DenseMat&, double);
    void Mat_Vec_Solve(const ifp_BlockVec& B, ifp_BlockVec& X) const
        {solve_is_mult(B, X);}
};

// B = A'
// must copy data
inline void ifp_DenseMat::Mat_Trans(ifp_LocalMat *B) const
{
    ifp_DenseMat& b = *(ifp_DenseMat *) B;

    assert (this != &b);       // in-place not allowed

    if (b.a == NULL && b.nrow == 0 && b.ncol == 0)
    {
	b.nrow = ncol;
	b.ncol = nrow;
	b.a = new double[b.nrow*b.ncol];
    }

    assert (a != NULL);
    assert (b.a != NULL);
    assert (b.nrow == ncol);
    assert (b.ncol == nrow);

    int i, j;
    double *p, *q;
    
    // Traverse B and fill from A.

    p = b.a;
    for (i=0; i<nrow; i++)
    {
	q = a+i;
        for (j=0; j<ncol; j++)
	{
	    *p++ = *q;
	    q += nrow;
	}
    }
}

// C = A + alpha B
inline void ifp_DenseMat::Mat_Mat_Add(const ifp_LocalMat *B, ifp_LocalMat *C,
    double alpha) const
{
    ifp_DenseMat& b = *(ifp_DenseMat *) B;
    ifp_DenseMat& c = *(ifp_DenseMat *) C;

    if (c.a == NULL && c.nrow == 0 && c.ncol == 0)
    {
	c.nrow = nrow;
	c.ncol = ncol;
	c.a = new double[c.nrow*c.ncol];
    }

    assert (a != NULL);
    assert (b.a != NULL);
    assert (c.a != NULL);
    assert (nrow == b.nrow);
    assert (ncol == b.ncol);
    assert (nrow == c.nrow);
    assert (ncol == c.ncol);

    int i;
    double *ap = a;
    double *bp = b.a;
    double *cp = c.a;

    if (alpha == 1.0)
        for (i=0; i<nrow*ncol; i++)
            *cp++ = *ap++ + *bp++;

    else if (alpha == -1.0)
        for (i=0; i<nrow*ncol; i++)
	    *cp++ = *ap++ - *bp++;

    else
        for (i=0; i<nrow*ncol; i++)
	    *cp++ = *ap++ + alpha * *bp++;
}

// C = alpha A B + beta C
inline void ifp_DenseMat::Mat_Mat_Mult(const ifp_LocalMat *B, ifp_LocalMat *C, 
    double alpha, double beta) const
{
    assert (B != C);
    ifp_DenseMat& b = *(ifp_DenseMat *) B;
    ifp_DenseMat& c = *(ifp_DenseMat *) C;

    if (beta == 0.0 && c.a == NULL && c.nrow == 0 && c.ncol == 0)
    {
	c.nrow = nrow;
	c.ncol = b.ncol;
	c.a = new double[c.nrow*c.ncol];
    }

    if (beta == 0.0 && (c.nrow != nrow || c.ncol != b.ncol))
    {
        assert (c.a != NULL);
	delete c.a;

	c.nrow = nrow;
	c.ncol = b.ncol;
	c.a = new double[c.nrow*c.ncol];
    }

    assert (a != NULL);
    assert (b.a != NULL);
    assert (c.a != NULL);

    char transa = 'N';
    char transb = 'N';
    integer M = nrow;
    integer N = c.ncol;
    integer K = ncol;
    integer LDA = nrow;
    integer LDB = b.nrow;
    integer LDC = c.nrow;

    F77NAME(dgemm)(&transa, &transb, &M, &N, &K, &alpha, a, &LDA, 
        b.a, &LDB, &beta, c.a, &LDC);
}

// C = alpha A B + beta C
inline void ifp_DenseMat::Mat_Vec_Mult(const ifp_BlockVec& B, ifp_BlockVec& C,
    double alpha, double beta) const
{
    char trans = 'N';
    integer M = nrow;
    integer N = C.dim1;
    integer K = ncol;
    integer LDA = nrow;
    integer LDB = B.dim0;
    integer LDC = C.dim0;

    assert (a != NULL);
    assert (&B != &C);

    F77NAME(dgemm)(&trans, &trans, &M, &N, &K, &alpha, a, &LDA, 
        B.v, &LDB, &beta, C.v, &LDC);
}

inline void ifp_DenseMat::Mat_Trans_Vec_Mult(const ifp_BlockVec& B, 
    ifp_BlockVec& C, double alpha, double beta) const
{
    char transa = 'T';
    char transb = 'N';
    integer M = ncol;
    integer N = C.dim1;
    integer K = nrow;
    integer LDA = nrow;
    integer LDB = B.dim0;
    integer LDC = C.dim0;

    assert (a != NULL);
    assert (&B != &C);

    F77NAME(dgemm)(&transa, &transb, &M, &N, &K, &alpha, a, &LDA,
        B.v, &LDB, &beta, C.v, &LDC);
}

inline void ifp_DenseMat::Mat_Vec_Solve(const ifp_BlockVec&, ifp_BlockVec&) const
{
    ifp_error("ifp_DenseMat::Mat_Vec_Solve: to solve with a local matrix,\n"
            "a local preconditioner should be used.\n", 0);
}

inline void ifp_DenseMat::Mat_Trans_Vec_Solve(const ifp_BlockVec&,
ifp_BlockVec&) const
{
    ifp_error("ifp_DenseMat::Mat_Trans_Vec_Solve: to solve with a local matrix,\n"
            "a local preconditioner should be used.\n", 0);
}

inline void ifp_DenseMat::SetToZero(int r, int c)
{
    if (a == NULL && nrow == 0 && ncol == 0)
    {
	nrow = r;
	ncol = c;
	a = new double[nrow*ncol];
    }

    assert (a != NULL);
    assert (nrow == r);
    assert (ncol == c);

    double *p = a;
    for (int i=0; i<nrow*ncol; i++)
	*p++ = 0.0;
}

inline ifp_LocalMat *ifp_DenseMat::CreateInv(ifp_LocalPrecon& local_precon) const
{
    assert (a != NULL);

    switch (local_precon.name)
    {
    case LP_LU:
        return new ifp_DenseMat_LU(*this);
    case LP_INVERSE:
        return new ifp_DenseMat_INVERSE(*this);
    case LP_SVD:
        return new ifp_DenseMat_SVD(*this, local_precon.darg1, local_precon.darg2);
    case LP_DIAGDOM:
        return new ifp_DenseMat_DIAGDOM(*this, local_precon.darg1);
    case LP_GERSH:
        return new ifp_DenseMat_GERSH(*this, local_precon.darg1);
    default:
        ifp_error("The local preconditioner you have chosen is not available\n"
                "for the type of blocks (ifp_DenseMat) in the block matrix.", 0);
    }
    return(0); // Zero pointer (never returned because ifp_error aborts, but this satisfies the compiler)
}

inline ifp_DenseMat_LU::ifp_DenseMat_LU(const ifp_DenseMat& A)
{
    a = NULL; // indicate this is implied inverse
    nrow = A.numrow();
    ncol = A.numcol();

    assert (nrow == ncol);

    lu = new double[nrow*ncol];
    // copy to lu

    double *p = lu;
    const double *q = &A(0,0);
    for (int i=0; i<nrow*ncol; i++)
        *p++ = *q++;

    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    ipiv = new integer[nrow];

    F77NAME(dgetrf)(&M, &N, lu, &LDA, ipiv, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_LU: dgetrf error", INFO);
}

inline void ifp_DenseMat_LU::Mat_Vec_Solve(const ifp_BlockVec& b, ifp_BlockVec& x) const
{
    char trans = 'N';
    integer N = ncol;
    integer NRHS = b.dim1;
    integer LDA = nrow;
    integer LDB = b.dim0;
    integer INFO;

    if (b.v != x.v)
        x.BlockCopy(b);

    F77NAME(dgetrs)(&trans, &N, &NRHS, lu, &LDA, ipiv, 
        x.v, &LDB, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_LU: dgetrs error", INFO);
}

inline ifp_DenseMat_LU::~ifp_DenseMat_LU()
{
    delete [] ipiv;
    delete [] lu;
}

inline ifp_DenseMat_INVERSE::ifp_DenseMat_INVERSE(const ifp_DenseMat& A)
    : ifp_DenseMat(A)
{
    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    assert (nrow == ncol);

    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_INVERSE: dgetrf error", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_INVERSE: dgetri error", INFO);

    delete [] ipiv;
    delete [] work;
}

inline ifp_DenseMat_SVD::ifp_DenseMat_SVD(const ifp_DenseMat& A, double rthresh, 
  double athresh)
    : ifp_DenseMat(A)
{
    assert (nrow == ncol);

    int i, j;
    double thresh;
    double *u  = new double[nrow*ncol];
    double *s  = new double[nrow];
    double *vt = new double[nrow*ncol];

    char job = 'A';
    integer N = nrow;
    integer LWORK = 5*N;
    double *work = new double[LWORK];
    integer INFO;
    integer nzero_diag = 0;
    for (i=0;i<nrow; i++) if (a[i*N+i] == 0.0) nzero_diag++;

    F77NAME(dgesvd) (&job, &job, &N, &N, a, &N, s, u, &N, vt, &N,
	work, &LWORK, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_SVD: dgesvd error", INFO);

    delete [] work;

    double s0 = s[0];
    double sn = s[nrow-1];

    // apply threshold
    thresh = s[0]*rthresh + athresh;
    int n_replaced = 0;
    for (i=0; i<nrow; i++)
	if (s[i] < thresh) 
	  {
	    s[i] = thresh;
	    n_replaced++;
	  }
    //   cout <<nzero_diag<<" "<<s0<<" "<<sn<<" "<<s0/sn<<" "<<thresh<<" "<<n_replaced << endl;


    if (s[0] == 0.0)
	ifp_error("ifp_DenseMat_SVD: block is zero after thresholding", 0);

    // scale the columns of u with reciprocal singular values
    double *p = u;
    for (i=0; i<ncol; i++)
    {
	double scal = s[i];
	for (j=0; j<nrow; j++)
	    *p++ /= scal;
    }

    // multiply through and store the result
    char trans = 'T';
    double alpha = 1.0;
    double beta = 0.0;
    F77NAME(dgemm)(&trans, &trans, &N, &N, &N, &alpha, vt, &N,
        u, &N, &beta, a, &N);

    delete [] u;
    delete [] s;
    delete [] vt;
}

inline ifp_DenseMat_DIAGDOM::ifp_DenseMat_DIAGDOM(const ifp_DenseMat& A, double alpha)
    : ifp_DenseMat(A)
{
    int i, j;
    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    assert (nrow == ncol);

    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // compute sum of abs of off-diagonals and put in work array
    double *p = a;
    for (i=0; i<nrow; i++)
        work[i] = 0.0;

    for (j=0; j<ncol; j++)
    {
        for (i=0; i<nrow; i++)
        {
            if (i != j)
            {
                work[i] += ABS(*p);
            }
            p++;
        }
    }

    for (i=0; i<nrow; i++)
    {
	if (ABS(a[i*nrow+i]) < alpha*work[i])
	    a[i*nrow+i] = SGN(a[i*nrow+i])*alpha*work[i];
    }

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_DIAGDOM: dgetrf error", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_DIAGDOM: dgetri error", INFO);

    delete [] ipiv;
    delete [] work;
}

inline ifp_DenseMat_GERSH::ifp_DenseMat_GERSH(const ifp_DenseMat& A, double alpha)
    : ifp_DenseMat(A)
{
    int i, j;
    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    assert (nrow == ncol);

    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // compute sum of abs of off-diagonals and put in work array
    double *p = a;
    for (i=0; i<nrow; i++)
        work[i] = 0.0;

    for (j=0; j<ncol; j++)
    {
        for (i=0; i<nrow; i++)
        {
            if (i != j)
            {
                work[i] += ABS(*p);
            }
            p++;
        }
    }

    double aii;
    for (i=0; i<nrow; i++)
    {
	aii = a[i*nrow+i];

	if (aii >= 0.0)
        {
	    if (aii - work[i] < alpha)
		a[i*nrow+i] += alpha - aii + work[i];
	}
	else
	{
	    if (aii + work[i] > -alpha)
		a[i*nrow+i] -= alpha + aii + work[i];
	}
    }

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_GERSH: dgetrf error", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
	ifp_error("ifp_DenseMat_GERSH: dgetri error", INFO);

    delete [] ipiv;
    delete [] work;
}

#endif // _IFP_DENSEMAT_H_
