/*

    Banded matrices are stored row-wise. Thus, element A(i,j)
   is stored at location j-i+j0 + M*i, where
   
      j0 = bl
      M  = bl+bu+1

    This means that element A(i,i-bl) is stored at 0 + M*i, and
    element A(i,i+bu) is stored at M-1 + M*i.

*/

// General purpose routines
r8 *c8band_new(int n, int bl, int bu, void *buf);
void c8band_print(int n, int bl, int bu, r8 *A);

// Access a particular element
r8 *c8band_element(int n, int bl, int bu, r8 *A, int i, int j);
c8 *cband_element(int n, int bl, int bu, c8 *A, int i, int j);

// Compute the matrix-vector product y = A*x
void c8band_matvec(int n, int bu, int bl, r8 *A, r8 *x, r8 *y);
void c8band_matvec1(int n, int bu, int bl, r8 *A, r8 *x, r8 *y);

// Compute the eigenvalues and eigenvectors of a banded Hermitian matrix
void c8band_eigen(int n, int b, r8 *A, int kstart, int kcount, r8 *e, r8 *v, int vstride);

// Return stats for last call to c8band_eigen
r8 c8band_eigen_mflops();
r8 c8band_eigen_time();

// Support routines for the eigensolver
int  c8band_oracle(int n, int b, r8 *A, r8 x, r8 *tmp);
void c8band_gershgorin(int n, int b, r8 *A, r8 *emin, r8 *emax);
void c8band_trace(int n, int b, r8 *A, r8 *tr, r8 *e2);
void c8band_invit(int n, int bl, int bu, r8 *r, r8 *y, r8 *x);
void c8band_qr_shifted(int n, int bl, int bu, r8 *A, c8 z, r8 *Q, r8 *R);

void c8band_q_apply(int n, int bl, r8 *Q, r8 *x);
void c8band_qh_apply(int n, int bl, r8 *Q, r8 *x);

// Householder routines
void c8_householder(int n, r8 *a, r8 *v);
void c8_apply_householder(int n, r8 *v, r8 *x);

// Generate a test matrix.
void c8band_test1(int n, int b, r8 *A);

