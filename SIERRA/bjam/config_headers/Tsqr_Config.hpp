//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __TSQR_Tsqr_Config_hpp
#define __TSQR_Tsqr_Config_hpp

/* #undef HAVE_LAPACK_DLARFP */
/* #undef HAVE_LAPACK_SLARFP */
/* #undef HAVE_LAPACK_ZLARFP */
/* #undef HAVE_LAPACK_CLARFP */

/* #undef HAVE_LAPACK_DLARFGP */
/* #undef HAVE_LAPACK_SLARFGP */
/* #undef HAVE_LAPACK_ZLARFGP */
/* #undef HAVE_LAPACK_CLARFGP */

/* #undef HAVE_LAPACK_DGEQRFP */
/* #undef HAVE_LAPACK_SGEQRFP */
/* #undef HAVE_LAPACK_ZGEQRFP */
/* #undef HAVE_LAPACK_CGEQRFP */

/* #undef HAVE_LAPACK_DGEQR2P */
/* #undef HAVE_LAPACK_SGEQR2P */
/* #undef HAVE_LAPACK_ZGEQR2P */
/* #undef HAVE_LAPACK_CGEQR2P */

#define F77_BLAS_MANGLE(name,NAME) name ## _

/* #undef HAVE_TSQR_NEW_MPI_COMMUNICATORS */
/* #undef HAVE_MPI_COMM_NETWORK */
/* #undef HAVE_MPI_COMM_NODE */
/* #undef HAVE_MPI_COMM_SOCKET */
/* #undef HAVE_MPI_COMM_CACHE */

#endif // __TSQR_Tsqr_Config_hpp
