/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Include file that should be included in every header file.           */
/* ******************************************************************** */
/* Author        : Jonathan Hu                                          */
/* Date          : July, 2003                                           */
/* ******************************************************************** */

#ifndef __MLCOMMON__
#define __MLCOMMON__

/* this avoids the classic build system from picking up a spurious
 * macro from other packages that may have been autotooled */
#ifdef ML_CLASSIC_BUILD
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif
#endif

/* The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
   be undef'd here to avoid warnings when this file is included from another package.
   Based on what Kevin Long did for Epetra. */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#ifdef HAVE_CONFIG_H 
#include "ml_config.h"

/* aztecoo depends on epetra ...*/
#if defined(HAVE_ML_AZTEC) || defined(HAVE_ML_EPETRA)
#define AZTEC
#define ML_WITH_EPETRA
#endif

/* ... but not vice versa */
#ifdef HAVE_ML_EPETRA
#ifndef ML_WITH_EPETRA
#define ML_WITH_EPETRA
#endif
#endif

#ifdef HAVE_ML_SUPERLU
#define SUPERLU
#endif

#ifdef HAVE_MPI
#define ML_MPI
#endif

#if defined(HAVE_ML_EXTERNAL_MPI_FUNCTIONS) && defined(HAVE_MPI)
#define ML_USING_MPI_FUNCTIONS
#endif

#ifdef HAVE_BLAS
#define USE_VENDOR_BLAS
#endif

#ifdef HAVE_LAPACK
#define USE_VENDOR_LAPACK
#endif

#ifdef HAVE_ML_ENRICH
#define ML_ENRICH
#endif

#ifdef HAVE_ML_NEW_T_PE
#define ML_NEW_T_PE
#endif

#ifdef HAVE_ML_COMPLEX_MAXWELL
#define GREG
#endif

#ifdef HAVE_ML_TIMING
#define ML_TIMING
#endif

#ifdef HAVE_ML_MULTIPLE_RHS
#define WKC
#endif

#ifdef HAVE_ML_FLOPS
#define ML_FLOPS
#endif

#ifdef HAVE_ML_BENCHMARKING
#define ML_BENCHMARK
#endif

#endif /*ifdef HAVE_CONFIG_H*/

#endif
