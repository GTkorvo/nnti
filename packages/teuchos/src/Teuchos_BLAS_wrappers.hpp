// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_BLAS_WRAPPERS_HPP_
#define _TEUCHOS_BLAS_WRAPPERS_HPP_

#include "Teuchos_ConfigDefs.hpp"

/*! \file Teuchos_BLAS_wrappers.hpp  
    \brief The Templated BLAS wrappers.
*/

/* Define fcd (Fortran Teuchos_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define Teuchos_fcd fcd 

#define DROTG_F77   F77_BLAS_MANGLE(srotg,SROTG)
#define DROT_F77    F77_BLAS_MANGLE(srot,SROT)
#define DASUM_F77   F77_BLAS_MANGLE(sasum,SASUM) 
#define DAXPY_F77   F77_BLAS_MANGLE(saxpy,SAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(scopy,SCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(sdot,SDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(snrm2,SNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(sscal,SSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(isamax,ISAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(sgemv,SGEMV)
#define DGER_F77    F77_BLAS_MANGLE(sger,SGER)
#define DTRMV_F77   F77_BLAS_MANGLE(strmv,STRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(sgemm,SGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(ssymm,SSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(strmm,STRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(strsm,STRSM)
 
#ifdef HAVE_TEUCHOS_COMPLEX

#define ZROTG_F77   F77_BLAS_MANGLE(crotg,CROTG)
#define ZROT_F77    F77_BLAS_MANGLE(crot,CROT)
#define ZASUM_F77   F77_BLAS_MANGLE(scasum,SCASUM) 
#define ZAXPY_F77   F77_BLAS_MANGLE(caxpy,CAXPY)
#define ZCOPY_F77   F77_BLAS_MANGLE(ccopy,CCOPY)
#define ZDOT_F77    F77_BLAS_MANGLE(cdotc,CDOTC)
#define ZNRM2_F77   F77_BLAS_MANGLE(scnrm2,SCNRM2)
#define ZSCAL_F77   F77_BLAS_MANGLE(cscal,CSCAL)
#define IZAMAX_F77  F77_BLAS_MANGLE(icamax,ICAMAX)
#define ZGEMV_F77   F77_BLAS_MANGLE(cgemv,CGEMV)
#define ZGER_F77    F77_BLAS_MANGLE(cgerc,CGERC)
#define ZTRMV_F77   F77_BLAS_MANGLE(ctrmv,CTRMV)
#define ZGEMM_F77   F77_BLAS_MANGLE(cgemm,CGEMM)
#define ZSYMM_F77   F77_BLAS_MANGLE(csymm,CSYMM)
#define ZTRMM_F77   F77_BLAS_MANGLE(ctrmm,CTRMM)
#define ZTRSM_F77   F77_BLAS_MANGLE(ctrsm,CTRSM)

#endif /* HAVE_TEUCHOS_COMPLEX */

#elif defined(INTEL_CXML)

#define PREFIX __stdcall 
#define Teuchos_fcd const char *, unsigned int 

#define DROTG_F77   F77_BLAS_MANGLE(drotg,DROTG)
#define DROT_F77    F77_BLAS_MANGLE(drot,DROT)
#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)  
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZROTG_F77   F77_BLAS_MANGLE(zrotg,ZROTG)
#define ZROT_F77    F77_BLAS_MANGLE(zrot,ZROT)
#define ZASUM_F77   F77_BLAS_MANGLE(dzasum,DZASUM)
#define ZAXPY_F77   F77_BLAS_MANGLE(zaxpy,ZAXPY)
#define ZCOPY_F77   F77_BLAS_MANGLE(zcopy,ZCOPY)
#define ZDOT_F77    F77_BLAS_MANGLE(zdotc,ZDOTC)  
#define ZNRM2_F77   F77_BLAS_MANGLE(dznrm2,DZNRM2)
#define ZSCAL_F77   F77_BLAS_MANGLE(zscal,ZSCAL)
#define IZAMAX_F77  F77_BLAS_MANGLE(izamax,IZAMAX)
#define ZGEMV_F77   F77_BLAS_MANGLE(zgemv,ZGEMV)
#define ZGER_F77    F77_BLAS_MANGLE(zgerc,ZGERC)
#define ZTRMV_F77   F77_BLAS_MANGLE(ztrmv,ZTRMV)
#define ZGEMM_F77   F77_BLAS_MANGLE(zgemm,ZGEMM)
#define ZSYMM_F77   F77_BLAS_MANGLE(zsymm,ZSYMM)
#define ZTRMM_F77   F77_BLAS_MANGLE(ztrmm,ZTRMM)
#define ZTRSM_F77   F77_BLAS_MANGLE(ztrsm,ZTRSM)

#endif /* HAVE_TEUCHOS_COMPLEX */

#elif defined(INTEL_MKL)

#define PREFIX
#define Teuchos_fcd const char *

#define DROTG_F77   F77_BLAS_MANGLE(drotg,DROTG)
#define DROT_F77    F77_BLAS_MANGLE(drot,DROT)
#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)  
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZROTG_F77   F77_BLAS_MANGLE(zrotg,ZROTG)
#define ZROT_F77    F77_BLAS_MANGLE(zrot,ZROT)
#define ZASUM_F77   F77_BLAS_MANGLE(dzasum,DZASUM)
#define ZAXPY_F77   F77_BLAS_MANGLE(zaxpy,ZAXPY)
#define ZCOPY_F77   F77_BLAS_MANGLE(zcopy,ZCOPY)
#define ZDOT_F77    F77_BLAS_MANGLE(zdotc,ZDOTC)  
#define ZNRM2_F77   F77_BLAS_MANGLE(dznrm2,DZNRM2)
#define ZSCAL_F77   F77_BLAS_MANGLE(zscal,ZSCAL)
#define IZAMAX_F77  F77_BLAS_MANGLE(izamax,IZAMAX)
#define ZGEMV_F77   F77_BLAS_MANGLE(zgemv,ZGEMV)
#define ZGER_F77    F77_BLAS_MANGLE(zgerc,ZGERC)
#define ZTRMV_F77   F77_BLAS_MANGLE(ztrmv,ZTRMV)
#define ZGEMM_F77   F77_BLAS_MANGLE(zgemm,ZGEMM)
#define ZSYMM_F77   F77_BLAS_MANGLE(zsymm,ZSYMM)
#define ZTRMM_F77   F77_BLAS_MANGLE(ztrmm,ZTRMM)
#define ZTRSM_F77   F77_BLAS_MANGLE(ztrsm,ZTRSM)

#endif /* HAVE_TEUCHOS_COMPLEX */

#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_BLAS_MANGLE is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif

#define F77_FUNC(lcase,UCASE) PREFIX UCASE
#define F77_BLAS_MANGLE(lcase,UCASE) PREFIX UCASE

#else /* Define Teuchos_fcd for all other machines */

#define PREFIX
#define Teuchos_fcd const char * 

#ifdef TRILINOS_NO_CONFIG_H

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE 
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) lcase
#define F77_BLAS_MANGLE(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) lcase ## _
#define F77_BLAS_MANGLE(lcase,UCASE) lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */

#endif /* TRILINOS_NO_CONFIG_H */

#define DROTG_F77   F77_BLAS_MANGLE(drotg,DROTG)
#define DROT_F77    F77_BLAS_MANGLE(drot,DROT)
#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM) 
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZROTG_F77   F77_BLAS_MANGLE(zrotg,ZROTG)
#define ZROT_F77    F77_BLAS_MANGLE(zrot,ZROT)
#define ZASUM_F77   F77_BLAS_MANGLE(dzasum,DZASUM)
#define ZAXPY_F77   F77_BLAS_MANGLE(zaxpy,ZAXPY)
#define ZCOPY_F77   F77_BLAS_MANGLE(zcopy,ZCOPY)
#define ZDOT_F77    F77_BLAS_MANGLE(zdotc,ZDOTC)  
#define ZNRM2_F77   F77_BLAS_MANGLE(dznrm2,DZNRM2)
#define ZSCAL_F77   F77_BLAS_MANGLE(zscal,ZSCAL)
#define IZAMAX_F77  F77_BLAS_MANGLE(izamax,IZAMAX)
#define ZGEMV_F77   F77_BLAS_MANGLE(zgemv,ZGEMV)
#define ZGER_F77    F77_BLAS_MANGLE(zgerc,ZGERC)
#define ZTRMV_F77   F77_BLAS_MANGLE(ztrmv,ZTRMV)
#define ZGEMM_F77   F77_BLAS_MANGLE(zgemm,ZGEMM)
#define ZSYMM_F77   F77_BLAS_MANGLE(zsymm,ZSYMM)
#define ZTRMM_F77   F77_BLAS_MANGLE(ztrmm,ZTRMM)
#define ZTRSM_F77   F77_BLAS_MANGLE(ztrsm,ZTRSM)

#endif /* HAVE_TEUCHOS_COMPLEX */

#endif


/* Explicitly define each F77 name for all BLAS kernels */

#define SROTG_F77   F77_BLAS_MANGLE(srotg,SROTG)
#define SROT_F77    F77_BLAS_MANGLE(srot,SROT)
#define SSCAL_F77   F77_BLAS_MANGLE(sscal,SSCAL) 
#define SCOPY_F77   F77_BLAS_MANGLE(scopy,SCOPY)
#define SAXPY_F77   F77_BLAS_MANGLE(saxpy,SAXPY)
#define SDOT_F77    F77_BLAS_MANGLE(sdot,SDOT)
#define SNRM2_F77   F77_BLAS_MANGLE(snrm2,SNRM2)
#define SASUM_F77   F77_BLAS_MANGLE(sasum,SASUM)
#define ISAMAX_F77  F77_BLAS_MANGLE(isamax,ISAMAX)

#define SGEMV_F77   F77_BLAS_MANGLE(sgemv,SGEMV)
#define SGER_F77    F77_BLAS_MANGLE(sger,SGER)
#define STRMV_F77   F77_BLAS_MANGLE(strmv,STRMV)
#define SGEMM_F77   F77_BLAS_MANGLE(sgemm,SGEMM)
#define SSYMM_F77   F77_BLAS_MANGLE(ssymm,SSYMM)
#define STRMM_F77   F77_BLAS_MANGLE(strmm,STRMM)
#define STRSM_F77   F77_BLAS_MANGLE(strsm,STRSM)

#ifdef HAVE_TEUCHOS_COMPLEX

#define CROTG_F77   F77_BLAS_MANGLE(crotg,CROTG)
#define CROT_F77    F77_BLAS_MANGLE(crot,CROT)
#define CASUM_F77   F77_BLAS_MANGLE(scasum,SCASUM) 
#define CAXPY_F77   F77_BLAS_MANGLE(caxpy,CAXPY)
#define CCOPY_F77   F77_BLAS_MANGLE(ccopy,CCOPY)
#define CDOT_F77    F77_BLAS_MANGLE(cdotc,CDOTC)
#define CNRM2_F77   F77_BLAS_MANGLE(scnrm2,SCNRM2)
#define CSCAL_F77   F77_BLAS_MANGLE(cscal,CSCAL)
#define ICAMAX_F77  F77_BLAS_MANGLE(icamax,ICAMAX)
#define CGEMV_F77   F77_BLAS_MANGLE(cgemv,CGEMV)
#define CGER_F77    F77_BLAS_MANGLE(cgerc,CGERC)
#define CTRMV_F77   F77_BLAS_MANGLE(ctrmv,CTRMV)
#define CGEMM_F77   F77_BLAS_MANGLE(cgemm,CGEMM)
#define CSYMM_F77   F77_BLAS_MANGLE(csymm,CSYMM)
#define CTRMM_F77   F77_BLAS_MANGLE(ctrmm,CTRMM)
#define CTRSM_F77   F77_BLAS_MANGLE(ctrsm,CTRSM)

#endif /* HAVE_TEUCHOS_COMPLEX */

#ifdef __cplusplus
extern "C" {
#endif


/* Double precision BLAS 1 */
void PREFIX DROTG_F77(double* da, double* db, double* c, double* s);
void PREFIX DROT_F77(const int* n, double* dx, const int* incx, double* dy, const int* incy, double* c, double* s);
double PREFIX DASUM_F77(const int* n, const double x[], const int* incx);
void PREFIX DAXPY_F77(const int* n, const double* alpha, const double x[], const int* incx, double y[], const int* incy);
void PREFIX DCOPY_F77(const int* n, const double *x, const int* incx, double *y, const int* incy);
double PREFIX DDOT_F77(const int* n, const double x[], const int* incx, const double y[], const int* incy);
double PREFIX DNRM2_F77(const int* n, const double x[], const int* incx); 
void PREFIX DSCAL_F77(const int* n, const double* alpha, double *x, const int* incx);
int PREFIX IDAMAX_F77(const int* n, const double *x, const int* incx);

/* Double std::complex precision BLAS 1 */
#ifdef HAVE_TEUCHOS_COMPLEX

void PREFIX ZROTG_F77(std::complex<double>* da, std::complex<double>* db, double* c, std::complex<double>* s);
void PREFIX ZROT_F77(const int* n, std::complex<double>* dx, const int* incx, std::complex<double>* dy, const int* incy, double* c, std::complex<double>* s);
double PREFIX ZASUM_F77(const int* n, const std::complex<double> x[], const int* incx);
void PREFIX ZAXPY_F77(const int* n, const std::complex<double>* alpha, const std::complex<double> x[], const int* incx, std::complex<double> y[], const int* incy);
void PREFIX ZCOPY_F77(const int* n, const std::complex<double> *x, const int* incx, std::complex<double> *y, const int* incy);
std::complex<double>
PREFIX ZDOT_F77(const int* n, const std::complex<double> x[], const int* incx, const std::complex<double> y[], const int* incy);
double PREFIX ZNRM2_F77(const int* n, const std::complex<double> x[], const int* incx); 
void PREFIX ZSCAL_F77(const int* n, const std::complex<double>* alpha, std::complex<double> *x, const int* incx);
int PREFIX IZAMAX_F77(const int* n, const std::complex<double> *x, const int* incx);

#endif // HAVE_TEUCHOS_COMPLEX

#ifdef HAVE_TEUCHOS_BLASFLOAT

/* Single precision BLAS 1 */ 
void PREFIX SROTG_F77(float* da, float* db, float* c, float* s);
void PREFIX SROT_F77(const int* n, float* dx, const int* incx, float* dy, const int* incy, float* c, float* s);
float PREFIX SASUM_F77(const int* n, const float x[], const int* incx);
void PREFIX SAXPY_F77(const int* n, const float* alpha, const float x[], const int* incx, float y[], const int* incy);
void PREFIX SCOPY_F77(const int* n, const float *x, const int* incx, float *y, const int* incy);
float PREFIX SDOT_F77(const int* n, const float x[], const int* incx, const float y[], const int* incy);
float PREFIX SNRM2_F77(const int* n, const float x[], const int* incx); 
void PREFIX SSCAL_F77(const int* n, const float* alpha, float *x, const int* incx);
int PREFIX ISAMAX_F77(const int* n, const float *x, const int* incx);

#endif // HAVE_TEUCHOS_BLASFLOAT

/* Single std::complex precision BLAS 1 */ 
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_BLASFLOAT)

void PREFIX CROTG_F77(std::complex<float>* da, std::complex<float>* db, float* c, std::complex<float>* s);
void PREFIX CROT_F77(const int* n, std::complex<float>* dx, const int* incx, std::complex<float>* dy, const int* incy, float* c, std::complex<float>* s);
float PREFIX CASUM_F77(const int* n, const std::complex<float> x[], const int* incx);
void PREFIX CAXPY_F77(const int* n, const std::complex<float>* alpha, const std::complex<float> x[], const int* incx, std::complex<float> y[], const int* incy);
void PREFIX CCOPY_F77(const int* n, const std::complex<float> *x, const int* incx, std::complex<float> *y, const int* incy);
std::complex<float> PREFIX CDOT_F77(const int* n, const std::complex<float> x[], const int* incx, const std::complex<float> y[], const int* incy);
float PREFIX CNRM2_F77(const int* n, const std::complex<float> x[], const int* incx); 
void PREFIX CSCAL_F77(const int* n, const std::complex<float>* alpha, std::complex<float> *x, const int* incx);
int PREFIX ICAMAX_F77(const int* n, const std::complex<float> *x, const int* incx);

#endif //  defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_BLASFLOAT)

/* Double precision BLAS 2 */
void PREFIX DGEMV_F77(Teuchos_fcd, const int* m, const int* n, const double* alpha, const double A[], const int* lda,
                 const double x[], const int* incx, const double* beta, double y[], const int* incy);
void PREFIX DTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n, 
                const double *a, const int *lda, double *x, const int *incx); 
void PREFIX DGER_F77(const int *m, const int *n, const double *alpha, const double *x, const int *incx, const double *y,
               const int *incy, double *a, const int *lda);

/* Double precision BLAS 2 */
#ifdef HAVE_TEUCHOS_COMPLEX

void PREFIX ZGEMV_F77(Teuchos_fcd, const int* m, const int* n, const std::complex<double>* alpha, const std::complex<double> A[], const int* lda,
                 const std::complex<double> x[], const int* incx, const std::complex<double>* beta, std::complex<double> y[], const int* incy);
void PREFIX ZTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n, 
                const std::complex<double> *a, const int *lda, std::complex<double> *x, const int *incx); 
void PREFIX ZGER_F77(const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *x, const int *incx, const std::complex<double> *y,
               const int *incy, std::complex<double> *a, const int *lda);

#endif /* HAVE_TEUCHOS_COMPLEX */

#ifdef HAVE_TEUCHOS_BLASFLOAT

/* Single precision BLAS 2 */
void PREFIX SGEMV_F77(Teuchos_fcd, const int* m, const int* n, const float* alpha, const float A[], const int* lda,
                 const float x[], const int* incx, const float* beta, float y[], const int* incy);
void PREFIX STRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,
                const float *a, const int *lda, float *x, const int *incx); 
void PREFIX SGER_F77(const int *m, const int *n, const float *alpha, const float *x, const int *incx, const float *y,
               const int *incy, float *a, const int *lda);

#endif // HAVE_TEUCHOS_BLASFLOAT

/* Single std::complex precision BLAS 2 */
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_BLASFLOAT)

void PREFIX CGEMV_F77(Teuchos_fcd, const int* m, const int* n, const std::complex<float>* alpha, const std::complex<float> A[], const int* lda,
                 const std::complex<float> x[], const int* incx, const std::complex<float>* beta, std::complex<float> y[], const int* incy);
void PREFIX CTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,
                const std::complex<float> *a, const int *lda, std::complex<float> *x, const int *incx); 
void PREFIX CGER_F77(const int *m, const int *n, const std::complex<float> *alpha, const std::complex<float> *x, const int *incx, const std::complex<float> *y,
               const int *incy, std::complex<float> *a, const int *lda);

#endif // defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_BLASFLOAT)

/* Double precision BLAS 3 */
void PREFIX DGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * 
                n, const int *k, const double *alpha, const double *a, const int *lda, 
                const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
void PREFIX DSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const double *alpha, const double *a, const int *lda, 
                const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
void PREFIX DTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,  
                const int *m, const int *n, const double *alpha, const double *a, const int * lda, double *b, const int *ldb);
void PREFIX DTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, 
                const int *m, const int *n, const double *alpha, const double *a, const int *
                lda, double *b, const int *ldb);

/* Double std::complex precision BLAS 3 */
#ifdef HAVE_TEUCHOS_COMPLEX

void PREFIX ZGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * 
                n, const int *k, const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, 
                const std::complex<double> *b, const int *ldb, const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
void PREFIX ZSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, 
                const std::complex<double> *b, const int *ldb, const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
void PREFIX ZTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,  
                const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *a, const int * lda, std::complex<double> *b, const int *ldb);
void PREFIX ZTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, 
                const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *a, const int *
                lda, std::complex<double> *b, const int *ldb);

#endif /* HAVE_TEUCHOS_COMPLEX */

#ifdef HAVE_TEUCHOS_BLASFLOAT

/* Single precision BLAS 3 */
void PREFIX SGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *
                n, const int *k, const float *alpha, const float *a, const int *lda, 
                const float *b, const int *ldb, const float *beta, float *c, const int *ldc);
void PREFIX SSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const float *alpha, const float *a, const int *lda, 
                const float *b, const int *ldb, const float *beta, float *c, const int *ldc);
void PREFIX STRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, 
                const int *m, const int *n, const float *alpha, const float *a, const int * lda, float *b, const int *ldb);
void PREFIX STRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const float *alpha, const float *a, const int *
                lda, float *b, const int *ldb);

#endif // HAVE_TEUCHOS_BLASFLOAT

/* Single std::complex precision BLAS 3 */

#ifdef HAVE_TEUCHOS_COMPLEX

void PREFIX CGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *
                n, const int *k, const std::complex<float> *alpha, const std::complex<float> *a, const int *lda, 
                const std::complex<float> *b, const int *ldb, const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void PREFIX CSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const std::complex<float> *alpha, const std::complex<float> *a, const int *lda, 
                const std::complex<float> *b, const int *ldb, const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void PREFIX CTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, 
                const int *m, const int *n, const std::complex<float> *alpha, const std::complex<float> *a, const int * lda, std::complex<float> *b, const int *ldb);
void PREFIX CTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const std::complex<float> *alpha, const std::complex<float> *a, const int *
                lda, std::complex<float> *b, const int *ldb);

#endif /* HAVE_TEUCHOS_COMPLEX */

#ifdef __cplusplus
}
#endif

#endif // end of TEUCHOS_BLAS_WRAPPERS_HPP_
