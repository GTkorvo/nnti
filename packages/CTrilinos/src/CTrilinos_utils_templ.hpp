
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


#ifndef CTRILINOS_UTILS_TEMPL_HPP
#define CTRILINOS_UTILS_TEMPL_HPP


#include "CTrilinos_config.h"


#include <string>

#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CTeuchos_CommandLineProcessor_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_any_Cpp.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CAmesos_Cpp.hpp"
#include "CEpetra_FECrsMatrix_Cpp.hpp"
#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#include "CAztecOO_Cpp.hpp"
#include "CAztecOO_StatusTest_Cpp.hpp"
#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#include "CIfpack_Cpp.hpp"
#include "CIfpack_Preconditioner_Cpp.hpp"

#include "CTrilinos_enums.h"
#include "CTrilinos_Table.hpp"
#include "CTrilinos_utils.hpp"


namespace CTrilinos {


template <class T>
CTrilinos_Universal_ID_t
cast( Table<const T> &destTable, CTrilinos_Universal_ID_t id )
{
    CTrilinos_Universal_ID_t newid;
    newid.type = CT_Invalid_ID;
    newid.index = -1;
    newid.is_const = id.is_const;

    if (id.is_const) {
        switch (id.type) {
        case CT_Epetra_Distributor_ID:
            newid = destTable.cast(CEpetra::getConstDistributor(id));
            break;
        case CT_Epetra_SerialComm_ID:
            newid = destTable.cast(CEpetra::getConstSerialComm(id));
            break;
        case CT_Epetra_BLAS_ID:
            newid = destTable.cast(CEpetra::getConstBLAS(id));
            break;
        case CT_Epetra_Comm_ID:
            newid = destTable.cast(CEpetra::getConstComm(id));
            break;
        case CT_Epetra_Operator_ID:
            newid = destTable.cast(CEpetra::getConstOperator(id));
            break;
        case CT_Epetra_MultiVector_ID:
            newid = destTable.cast(CEpetra::getConstMultiVector(id));
            break;
        case CT_Epetra_OffsetIndex_ID:
            newid = destTable.cast(CEpetra::getConstOffsetIndex(id));
            break;
        case CT_Epetra_Object_ID:
            newid = destTable.cast(CEpetra::getConstObject(id));
            break;
        case CT_Epetra_RowMatrix_ID:
            newid = destTable.cast(CEpetra::getConstRowMatrix(id));
            break;
        case CT_Epetra_CompObject_ID:
            newid = destTable.cast(CEpetra::getConstCompObject(id));
            break;
        case CT_Epetra_Directory_ID:
            newid = destTable.cast(CEpetra::getConstDirectory(id));
            break;
        case CT_Epetra_Flops_ID:
            newid = destTable.cast(CEpetra::getConstFlops(id));
            break;
        case CT_Epetra_SrcDistObject_ID:
            newid = destTable.cast(CEpetra::getConstSrcDistObject(id));
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            newid = destTable.cast(CEpetra::getConstMpiComm(id));
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            newid = destTable.cast(CEpetra::getConstCrsMatrix(id));
            break;
        case CT_Epetra_CrsGraph_ID:
            newid = destTable.cast(CEpetra::getConstCrsGraph(id));
            break;
        case CT_Epetra_DistObject_ID:
            newid = destTable.cast(CEpetra::getConstDistObject(id));
            break;
        case CT_Epetra_Vector_ID:
            newid = destTable.cast(CEpetra::getConstVector(id));
            break;
        case CT_Epetra_Export_ID:
            newid = destTable.cast(CEpetra::getConstExport(id));
            break;
        case CT_Epetra_Map_ID:
            newid = destTable.cast(CEpetra::getConstMap(id));
            break;
        case CT_Epetra_BlockMap_ID:
            newid = destTable.cast(CEpetra::getConstBlockMap(id));
            break;
        case CT_Epetra_Import_ID:
            newid = destTable.cast(CEpetra::getConstImport(id));
            break;
        case CT_Epetra_Time_ID:
            newid = destTable.cast(CEpetra::getConstTime(id));
            break;
        case CT_Epetra_JadMatrix_ID:
            newid = destTable.cast(CEpetra::getConstJadMatrix(id));
            break;
        case CT_Epetra_LinearProblem_ID:
            newid = destTable.cast(CEpetra::getConstLinearProblem(id));
            break;
        case CT_Epetra_LAPACK_ID:
            newid = destTable.cast(CEpetra::getConstLAPACK(id));
            break;
        case CT_Teuchos_ParameterList_ID:
            newid = destTable.cast(CTeuchos::getConstParameterList(id));
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            newid = destTable.cast(CAmesos::getConstBaseSolver(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            newid = destTable.cast(CEpetra::getConstFECrsMatrix(id));
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            newid = destTable.cast(CEpetra::getConstIntSerialDenseVector(id));
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            newid = destTable.cast(CEpetra::getConstSerialDenseMatrix(id));
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            newid = destTable.cast(CAztecOO::getConstStatusTest(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            newid = destTable.cast(CAztecOO::getConstStatusTestCombo(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            newid = destTable.cast(CAztecOO::getConstStatusTestMaxIters(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            newid = destTable.cast(CAztecOO::getConstStatusTestResNorm(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            newid = destTable.cast(CIfpack::getConstPreconditioner(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            break;
        }
    } else {
        switch (id.type) {
        case CT_Epetra_Distributor_ID:
            newid = destTable.cast(CEpetra::getDistributor(id));
            break;
        case CT_Epetra_SerialComm_ID:
            newid = destTable.cast(CEpetra::getSerialComm(id));
            break;
        case CT_Epetra_BLAS_ID:
            newid = destTable.cast(CEpetra::getBLAS(id));
            break;
        case CT_Epetra_Comm_ID:
            newid = destTable.cast(CEpetra::getComm(id));
            break;
        case CT_Epetra_Operator_ID:
            newid = destTable.cast(CEpetra::getOperator(id));
            break;
        case CT_Epetra_MultiVector_ID:
            newid = destTable.cast(CEpetra::getMultiVector(id));
            break;
        case CT_Epetra_OffsetIndex_ID:
            newid = destTable.cast(CEpetra::getOffsetIndex(id));
            break;
        case CT_Epetra_Object_ID:
            newid = destTable.cast(CEpetra::getObject(id));
            break;
        case CT_Epetra_RowMatrix_ID:
            newid = destTable.cast(CEpetra::getRowMatrix(id));
            break;
        case CT_Epetra_CompObject_ID:
            newid = destTable.cast(CEpetra::getCompObject(id));
            break;
        case CT_Epetra_Directory_ID:
            newid = destTable.cast(CEpetra::getDirectory(id));
            break;
        case CT_Epetra_Flops_ID:
            newid = destTable.cast(CEpetra::getFlops(id));
            break;
        case CT_Epetra_SrcDistObject_ID:
            newid = destTable.cast(CEpetra::getSrcDistObject(id));
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            newid = destTable.cast(CEpetra::getMpiComm(id));
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            newid = destTable.cast(CEpetra::getCrsMatrix(id));
            break;
        case CT_Epetra_CrsGraph_ID:
            newid = destTable.cast(CEpetra::getCrsGraph(id));
            break;
        case CT_Epetra_DistObject_ID:
            newid = destTable.cast(CEpetra::getDistObject(id));
            break;
        case CT_Epetra_Vector_ID:
            newid = destTable.cast(CEpetra::getVector(id));
            break;
        case CT_Epetra_Export_ID:
            newid = destTable.cast(CEpetra::getExport(id));
            break;
        case CT_Epetra_Map_ID:
            newid = destTable.cast(CEpetra::getMap(id));
            break;
        case CT_Epetra_BlockMap_ID:
            newid = destTable.cast(CEpetra::getBlockMap(id));
            break;
        case CT_Epetra_Import_ID:
            newid = destTable.cast(CEpetra::getImport(id));
            break;
        case CT_Epetra_Time_ID:
            newid = destTable.cast(CEpetra::getTime(id));
            break;
        case CT_Epetra_JadMatrix_ID:
            newid = destTable.cast(CEpetra::getJadMatrix(id));
            break;
        case CT_Epetra_LinearProblem_ID:
            newid = destTable.cast(CEpetra::getLinearProblem(id));
            break;
        case CT_Epetra_LAPACK_ID:
            newid = destTable.cast(CEpetra::getLAPACK(id));
            break;
        case CT_Teuchos_ParameterList_ID:
            newid = destTable.cast(CTeuchos::getParameterList(id));
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            newid = destTable.cast(CAmesos::getBaseSolver(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            newid = destTable.cast(CEpetra::getFECrsMatrix(id));
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            newid = destTable.cast(CEpetra::getIntSerialDenseVector(id));
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            newid = destTable.cast(CEpetra::getSerialDenseMatrix(id));
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            newid = destTable.cast(CAztecOO::getStatusTest(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            newid = destTable.cast(CAztecOO::getStatusTestCombo(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            newid = destTable.cast(CAztecOO::getStatusTestMaxIters(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            newid = destTable.cast(CAztecOO::getStatusTestResNorm(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            newid = destTable.cast(CIfpack::getPreconditioner(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            break;
        }
    }

    return newid;
}

template <class T>
CTrilinos_Universal_ID_t
cast( Table<T> &destTable, CTrilinos_Universal_ID_t id )
{
    CTrilinos_Universal_ID_t newid;
    newid.type = CT_Invalid_ID;
    newid.index = -1;
    newid.is_const = id.is_const;

    if (id.is_const) {
        throw CTrilinosTypeMismatchError("[CTrilinos::cast]: casting away constness");
    } else {
        switch (id.type) {
        case CT_Epetra_Distributor_ID:
            newid = destTable.cast(CEpetra::getDistributor(id));
            break;
        case CT_Epetra_SerialComm_ID:
            newid = destTable.cast(CEpetra::getSerialComm(id));
            break;
        case CT_Epetra_BLAS_ID:
            newid = destTable.cast(CEpetra::getBLAS(id));
            break;
        case CT_Epetra_Comm_ID:
            newid = destTable.cast(CEpetra::getComm(id));
            break;
        case CT_Epetra_Operator_ID:
            newid = destTable.cast(CEpetra::getOperator(id));
            break;
        case CT_Epetra_MultiVector_ID:
            newid = destTable.cast(CEpetra::getMultiVector(id));
            break;
        case CT_Epetra_OffsetIndex_ID:
            newid = destTable.cast(CEpetra::getOffsetIndex(id));
            break;
        case CT_Epetra_Object_ID:
            newid = destTable.cast(CEpetra::getObject(id));
            break;
        case CT_Epetra_RowMatrix_ID:
            newid = destTable.cast(CEpetra::getRowMatrix(id));
            break;
        case CT_Epetra_CompObject_ID:
            newid = destTable.cast(CEpetra::getCompObject(id));
            break;
        case CT_Epetra_Directory_ID:
            newid = destTable.cast(CEpetra::getDirectory(id));
            break;
        case CT_Epetra_Flops_ID:
            newid = destTable.cast(CEpetra::getFlops(id));
            break;
        case CT_Epetra_SrcDistObject_ID:
            newid = destTable.cast(CEpetra::getSrcDistObject(id));
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            newid = destTable.cast(CEpetra::getMpiComm(id));
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            newid = destTable.cast(CEpetra::getCrsMatrix(id));
            break;
        case CT_Epetra_CrsGraph_ID:
            newid = destTable.cast(CEpetra::getCrsGraph(id));
            break;
        case CT_Epetra_DistObject_ID:
            newid = destTable.cast(CEpetra::getDistObject(id));
            break;
        case CT_Epetra_Vector_ID:
            newid = destTable.cast(CEpetra::getVector(id));
            break;
        case CT_Epetra_Export_ID:
            newid = destTable.cast(CEpetra::getExport(id));
            break;
        case CT_Epetra_Map_ID:
            newid = destTable.cast(CEpetra::getMap(id));
            break;
        case CT_Epetra_BlockMap_ID:
            newid = destTable.cast(CEpetra::getBlockMap(id));
            break;
        case CT_Epetra_Import_ID:
            newid = destTable.cast(CEpetra::getImport(id));
            break;
        case CT_Epetra_Time_ID:
            newid = destTable.cast(CEpetra::getTime(id));
            break;
        case CT_Epetra_JadMatrix_ID:
            newid = destTable.cast(CEpetra::getJadMatrix(id));
            break;
        case CT_Epetra_LinearProblem_ID:
            newid = destTable.cast(CEpetra::getLinearProblem(id));
            break;
        case CT_Epetra_LAPACK_ID:
            newid = destTable.cast(CEpetra::getLAPACK(id));
            break;
        case CT_Teuchos_ParameterList_ID:
            newid = destTable.cast(CTeuchos::getParameterList(id));
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            newid = destTable.cast(CAmesos::getBaseSolver(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            newid = destTable.cast(CEpetra::getFECrsMatrix(id));
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            newid = destTable.cast(CEpetra::getIntSerialDenseVector(id));
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            newid = destTable.cast(CEpetra::getSerialDenseMatrix(id));
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            newid = destTable.cast(CAztecOO::getStatusTest(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            newid = destTable.cast(CAztecOO::getStatusTestCombo(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            newid = destTable.cast(CAztecOO::getStatusTestMaxIters(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            newid = destTable.cast(CAztecOO::getStatusTestResNorm(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            newid = destTable.cast(CIfpack::getPreconditioner(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            break;
        }
    }

    return newid;
}

/* convert struct from specific type to generic CTrilinos_Universal_ID_t
 * but keep the content in tact */
template <typename T>
CTrilinos_Universal_ID_t
abstractType( T id )
{
    CTrilinos_Universal_ID_t newid;

    newid.type = id.type;
    newid.index = id.index;
    newid.is_const = id.is_const;

    return newid;
}

/* convert struct from generic CTrilinos_Universal_ID_t to specific type
 * but keep the content in tact */
template <typename T>
T
concreteType( CTrilinos_Universal_ID_t id )
{
    T newid;

    newid.type = id.type;
    newid.index = id.index;
    newid.is_const = id.is_const;

    return newid;
}

/* isSameObject(RCP, RCP) */
template <class T1, class T2>
bool
isSameObject( const Teuchos::RCP<T1> &rcp1, const Teuchos::RCP<T2> &rcp2 )
{
    return (rcp1.shares_resource(rcp2));
}

/* isSameObject(RCP, generic_id) */
template <class T>
bool
isSameObject( const Teuchos::RCP<T> &rcp, CTrilinos_Universal_ID_t id )
{
    bool shares = false;

    if (id.is_const) {
        switch (id.type) {
        case CT_Epetra_Distributor_ID:
            shares = rcp.shares_resource(CEpetra::getConstDistributor(id));
            break;
        case CT_Epetra_SerialComm_ID:
            shares = rcp.shares_resource(CEpetra::getConstSerialComm(id));
            break;
        case CT_Epetra_BLAS_ID:
            shares = rcp.shares_resource(CEpetra::getConstBLAS(id));
            break;
        case CT_Epetra_Comm_ID:
            shares = rcp.shares_resource(CEpetra::getConstComm(id));
            break;
        case CT_Epetra_Operator_ID:
            shares = rcp.shares_resource(CEpetra::getConstOperator(id));
            break;
        case CT_Epetra_MultiVector_ID:
            shares = rcp.shares_resource(CEpetra::getConstMultiVector(id));
            break;
        case CT_Epetra_OffsetIndex_ID:
            shares = rcp.shares_resource(CEpetra::getConstOffsetIndex(id));
            break;
        case CT_Epetra_Object_ID:
            shares = rcp.shares_resource(CEpetra::getConstObject(id));
            break;
        case CT_Epetra_RowMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getConstRowMatrix(id));
            break;
        case CT_Epetra_CompObject_ID:
            shares = rcp.shares_resource(CEpetra::getConstCompObject(id));
            break;
        case CT_Epetra_Directory_ID:
            shares = rcp.shares_resource(CEpetra::getConstDirectory(id));
            break;
        case CT_Epetra_Flops_ID:
            shares = rcp.shares_resource(CEpetra::getConstFlops(id));
            break;
        case CT_Epetra_SrcDistObject_ID:
            shares = rcp.shares_resource(CEpetra::getConstSrcDistObject(id));
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            shares = rcp.shares_resource(CEpetra::getConstMpiComm(id));
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getConstCrsMatrix(id));
            break;
        case CT_Epetra_CrsGraph_ID:
            shares = rcp.shares_resource(CEpetra::getConstCrsGraph(id));
            break;
        case CT_Epetra_DistObject_ID:
            shares = rcp.shares_resource(CEpetra::getConstDistObject(id));
            break;
        case CT_Epetra_Vector_ID:
            shares = rcp.shares_resource(CEpetra::getConstVector(id));
            break;
        case CT_Epetra_Export_ID:
            shares = rcp.shares_resource(CEpetra::getConstExport(id));
            break;
        case CT_Epetra_Map_ID:
            shares = rcp.shares_resource(CEpetra::getConstMap(id));
            break;
        case CT_Epetra_BlockMap_ID:
            shares = rcp.shares_resource(CEpetra::getConstBlockMap(id));
            break;
        case CT_Epetra_Import_ID:
            shares = rcp.shares_resource(CEpetra::getConstImport(id));
            break;
        case CT_Epetra_Time_ID:
            shares = rcp.shares_resource(CEpetra::getConstTime(id));
            break;
        case CT_Epetra_JadMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getConstJadMatrix(id));
            break;
        case CT_Epetra_LinearProblem_ID:
            shares = rcp.shares_resource(CEpetra::getConstLinearProblem(id));
            break;
        case CT_Epetra_LAPACK_ID:
            shares = rcp.shares_resource(CEpetra::getConstLAPACK(id));
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            shares = rcp.shares_resource(CTeuchos::getConstCommandLineProcessor(id));
            break;
        case CT_Teuchos_ParameterList_ID:
            shares = rcp.shares_resource(CTeuchos::getConstParameterList(id));
            break;
        case CT_Teuchos_ParameterEntry_ID:
            shares = rcp.shares_resource(CTeuchos::getConstParameterEntry(id));
            break;
        case CT_Teuchos_any_ID:
            shares = rcp.shares_resource(CTeuchos::getConstany(id));
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            shares = rcp.shares_resource(CAmesos::getConstBaseSolver(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            shares = rcp.shares_resource(CAmesos::getConstAmesos(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getConstFECrsMatrix(id));
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            shares = rcp.shares_resource(CEpetra::getConstIntSerialDenseVector(id));
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getConstSerialDenseMatrix(id));
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_ID:
            shares = rcp.shares_resource(CAztecOO::getConstAztecOO(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            shares = rcp.shares_resource(CAztecOO::getConstStatusTest(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            shares = rcp.shares_resource(CAztecOO::getConstStatusTestCombo(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            shares = rcp.shares_resource(CAztecOO::getConstStatusTestMaxIters(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            shares = rcp.shares_resource(CAztecOO::getConstStatusTestResNorm(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_ID:
            shares = rcp.shares_resource(CIfpack::getConstIfpack(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            shares = rcp.shares_resource(CIfpack::getConstPreconditioner(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            break;
        }
    } else {
        switch (id.type) {
        case CT_Epetra_Distributor_ID:
            shares = rcp.shares_resource(CEpetra::getDistributor(id));
            break;
        case CT_Epetra_SerialComm_ID:
            shares = rcp.shares_resource(CEpetra::getSerialComm(id));
            break;
        case CT_Epetra_BLAS_ID:
            shares = rcp.shares_resource(CEpetra::getBLAS(id));
            break;
        case CT_Epetra_Comm_ID:
            shares = rcp.shares_resource(CEpetra::getComm(id));
            break;
        case CT_Epetra_Operator_ID:
            shares = rcp.shares_resource(CEpetra::getOperator(id));
            break;
        case CT_Epetra_MultiVector_ID:
            shares = rcp.shares_resource(CEpetra::getMultiVector(id));
            break;
        case CT_Epetra_OffsetIndex_ID:
            shares = rcp.shares_resource(CEpetra::getOffsetIndex(id));
            break;
        case CT_Epetra_Object_ID:
            shares = rcp.shares_resource(CEpetra::getObject(id));
            break;
        case CT_Epetra_RowMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getRowMatrix(id));
            break;
        case CT_Epetra_CompObject_ID:
            shares = rcp.shares_resource(CEpetra::getCompObject(id));
            break;
        case CT_Epetra_Directory_ID:
            shares = rcp.shares_resource(CEpetra::getDirectory(id));
            break;
        case CT_Epetra_Flops_ID:
            shares = rcp.shares_resource(CEpetra::getFlops(id));
            break;
        case CT_Epetra_SrcDistObject_ID:
            shares = rcp.shares_resource(CEpetra::getSrcDistObject(id));
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            shares = rcp.shares_resource(CEpetra::getMpiComm(id));
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getCrsMatrix(id));
            break;
        case CT_Epetra_CrsGraph_ID:
            shares = rcp.shares_resource(CEpetra::getCrsGraph(id));
            break;
        case CT_Epetra_DistObject_ID:
            shares = rcp.shares_resource(CEpetra::getDistObject(id));
            break;
        case CT_Epetra_Vector_ID:
            shares = rcp.shares_resource(CEpetra::getVector(id));
            break;
        case CT_Epetra_Export_ID:
            shares = rcp.shares_resource(CEpetra::getExport(id));
            break;
        case CT_Epetra_Map_ID:
            shares = rcp.shares_resource(CEpetra::getMap(id));
            break;
        case CT_Epetra_BlockMap_ID:
            shares = rcp.shares_resource(CEpetra::getBlockMap(id));
            break;
        case CT_Epetra_Import_ID:
            shares = rcp.shares_resource(CEpetra::getImport(id));
            break;
        case CT_Epetra_Time_ID:
            shares = rcp.shares_resource(CEpetra::getTime(id));
            break;
        case CT_Epetra_JadMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getJadMatrix(id));
            break;
        case CT_Epetra_LinearProblem_ID:
            shares = rcp.shares_resource(CEpetra::getLinearProblem(id));
            break;
        case CT_Epetra_LAPACK_ID:
            shares = rcp.shares_resource(CEpetra::getLAPACK(id));
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            shares = rcp.shares_resource(CTeuchos::getCommandLineProcessor(id));
            break;
        case CT_Teuchos_ParameterList_ID:
            shares = rcp.shares_resource(CTeuchos::getParameterList(id));
            break;
        case CT_Teuchos_ParameterEntry_ID:
            shares = rcp.shares_resource(CTeuchos::getParameterEntry(id));
            break;
        case CT_Teuchos_any_ID:
            shares = rcp.shares_resource(CTeuchos::getany(id));
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            shares = rcp.shares_resource(CAmesos::getBaseSolver(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            shares = rcp.shares_resource(CAmesos::getAmesos(id));
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getFECrsMatrix(id));
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            shares = rcp.shares_resource(CEpetra::getIntSerialDenseVector(id));
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            shares = rcp.shares_resource(CEpetra::getSerialDenseMatrix(id));
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_ID:
            shares = rcp.shares_resource(CAztecOO::getAztecOO(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            shares = rcp.shares_resource(CAztecOO::getStatusTest(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            shares = rcp.shares_resource(CAztecOO::getStatusTestCombo(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            shares = rcp.shares_resource(CAztecOO::getStatusTestMaxIters(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            shares = rcp.shares_resource(CAztecOO::getStatusTestResNorm(id));
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_ID:
            shares = rcp.shares_resource(CIfpack::getIfpack(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            shares = rcp.shares_resource(CIfpack::getPreconditioner(id));
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        default:
            break;
        }
    }

    return shares;
}

/* isSameObject(RCP, specific_id) */
template <class T, typename S>
bool
isSameObject( const Teuchos::RCP<T> &rcp, S id )
{
    return isSameObject<T>(rcp, abstractType<S>(id));
}

/* isSameObject(specific_id, RCP) */
template <typename S, class T>
bool
isSameObject( S id, const Teuchos::RCP<T> &rcp )
{
    return isSameObject<T>(rcp, abstractType<S>(id));
}

/* isSameObject(specific_id, specific_id) */
template <typename S1, typename S2>
bool
isSameObject( S1 id1, S2 id2 )
{
    CTrilinos_Universal_ID_t id1a = abstractType<S1>(id1);
    CTrilinos_Universal_ID_t id2a = abstractType<S2>(id2);

    return isSameObject(id1a, id2a);
}


} // namespace CTrilinos


#endif

