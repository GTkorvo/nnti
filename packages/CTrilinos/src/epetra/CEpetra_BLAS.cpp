
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"

#include "CTrilinos_enums.h"
#include "CEpetra_BLAS.h"
#include "CEpetra_BLAS_Cpp.hpp"
#include "Epetra_BLAS.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_BLAS */
Table<Epetra_BLAS>& tableOfBLASs()
{
    static Table<Epetra_BLAS> loc_tableOfBLASs(CT_Epetra_BLAS_ID);
    return loc_tableOfBLASs;
}


} // namespace


//
// Definitions from CEpetra_BLAS.h
//


extern "C" {


CT_Epetra_BLAS_ID_t Epetra_BLAS_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_BLAS_Generalize ( 
  CT_Epetra_BLAS_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id);
}

CT_Epetra_BLAS_ID_t Epetra_BLAS_Create (  )
{
    return CEpetra::storeNewBLAS(new Epetra_BLAS());
}

CT_Epetra_BLAS_ID_t Epetra_BLAS_Duplicate ( 
  CT_Epetra_BLAS_ID_t BLASID )
{
    const Teuchos::RCP<const Epetra_BLAS> BLAS = CEpetra::getConstBLAS(BLASID);
    return CEpetra::storeNewBLAS(new Epetra_BLAS(*BLAS));
}

void Epetra_BLAS_Destroy ( CT_Epetra_BLAS_ID_t * selfID )
{
    CEpetra::removeBLAS(selfID);
}

float Epetra_BLAS_ASUM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->ASUM(N, X, INCX);
}

double Epetra_BLAS_ASUM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->ASUM(N, X, INCX);
}

float Epetra_BLAS_DOT_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const float * Y, const int INCX, const int INCY )
{
    return CEpetra::getConstBLAS(selfID)->DOT(N, X, Y, INCX, INCY);
}

double Epetra_BLAS_DOT_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const double * Y, const int INCX, const int INCY )
{
    return CEpetra::getConstBLAS(selfID)->DOT(N, X, Y, INCX, INCY);
}

float Epetra_BLAS_NRM2_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->NRM2(N, X, INCX);
}

double Epetra_BLAS_NRM2_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->NRM2(N, X, INCX);
}

void Epetra_BLAS_SCAL_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  float * X, const int INCX )
{
    CEpetra::getConstBLAS(selfID)->SCAL(N, ALPHA, X, INCX);
}

void Epetra_BLAS_SCAL_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  double * X, const int INCX )
{
    CEpetra::getConstBLAS(selfID)->SCAL(N, ALPHA, X, INCX);
}

void Epetra_BLAS_COPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  float * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->COPY(N, X, Y, INCX, INCY);
}

void Epetra_BLAS_COPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  double * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->COPY(N, X, Y, INCX, INCY);
}

int Epetra_BLAS_IAMAX_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->IAMAX(N, X, INCX);
}

int Epetra_BLAS_IAMAX_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->IAMAX(N, X, INCX);
}

void Epetra_BLAS_AXPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  const float * X, float * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->AXPY(N, ALPHA, X, Y, INCX, INCY);
}

void Epetra_BLAS_AXPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  const double * X, double * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->AXPY(N, ALPHA, X, Y, INCX, INCY);
}

void Epetra_BLAS_GEMV_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const float ALPHA, const float * A, const int LDA, 
  const float * X, const float BETA, float * Y, const int INCX, 
  const int INCY )
{
    CEpetra::getConstBLAS(selfID)->GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, 
        INCX, INCY);
}

void Epetra_BLAS_GEMV_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const double ALPHA, const double * A, const int LDA, 
  const double * X, const double BETA, double * Y, const int INCX, 
  const int INCY )
{
    CEpetra::getConstBLAS(selfID)->GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, 
        INCX, INCY);
}

void Epetra_BLAS_GEMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const float ALPHA, 
  const float * A, const int LDA, const float * B, const int LDB, 
  const float BETA, float * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->GEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, 
        B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_GEMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const double ALPHA, 
  const double * A, const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->GEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, 
        B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_SYMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const float ALPHA, const float * A, 
  const int LDA, const float * B, const int LDB, const float BETA, 
  float * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->SYMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, 
        LDB, BETA, C, LDC);
}

void Epetra_BLAS_SYMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const double ALPHA, const double * A, 
  const int LDA, const double * B, const int LDB, const double BETA, 
  double * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->SYMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, 
        LDB, BETA, C, LDC);
}

void Epetra_BLAS_TRMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const float ALPHA, const float * A, const int LDA, float * B, 
  const int LDB )
{
    CEpetra::getConstBLAS(selfID)->TRMM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, 
        A, LDA, B, LDB);
}

void Epetra_BLAS_TRMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const double ALPHA, const double * A, const int LDA, double * B, 
  const int LDB )
{
    CEpetra::getConstBLAS(selfID)->TRMM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, 
        A, LDA, B, LDB);
}


} // extern "C"


//
// Definitions from CEpetra_BLAS_Cpp.hpp
//


/* get Epetra_BLAS from non-const table using CT_Epetra_BLAS_ID */
const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CT_Epetra_BLAS_ID_t id )
{
    if (tableOfBLASs().isType(id.table))
        return tableOfBLASs().get<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

/* get Epetra_BLAS from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CTrilinos_Universal_ID_t id )
{
    if (tableOfBLASs().isType(id.table))
        return tableOfBLASs().get<Epetra_BLAS>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_BLAS>(id);
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CT_Epetra_BLAS_ID */
const Teuchos::RCP<const Epetra_BLAS>
CEpetra::getConstBLAS( CT_Epetra_BLAS_ID_t id )
{
    if (tableOfBLASs().isType(id.table))
        return tableOfBLASs().getConst<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BLAS>
CEpetra::getConstBLAS( CTrilinos_Universal_ID_t id )
{
    if (tableOfBLASs().isType(id.table))
        return tableOfBLASs().getConst<Epetra_BLAS>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_BLAS>(id);
}

/* store Epetra_BLAS (owned) in non-const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeNewBLAS( Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().store<Epetra_BLAS>(pobj, true));
}

/* store Epetra_BLAS in non-const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeBLAS( Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().store<Epetra_BLAS>(pobj, false));
}

/* store const Epetra_BLAS in const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeConstBLAS( const Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().store<Epetra_BLAS>(pobj, false));
}

/* remove Epetra_BLAS from table using CT_Epetra_BLAS_ID */
void
CEpetra::removeBLAS( CT_Epetra_BLAS_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(*id);
    if (tableOfBLASs().isType(aid.table))
        tableOfBLASs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(aid);
}

/* remove Epetra_BLAS from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeBLAS( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfBLASs().isType(aid->table))
        tableOfBLASs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_BLAS table */
void
CEpetra::purgeBLAS(  )
{
    tableOfBLASs().purge();
}

/* store Epetra_BLAS in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasBLAS( const Teuchos::RCP< Epetra_BLAS > & robj )
{
    return tableOfBLASs().alias(robj);
}

/* store const Epetra_BLAS in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstBLAS( const Teuchos::RCP< const Epetra_BLAS > & robj )
{
    return tableOfBLASs().alias(robj);
}



