
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

#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_Operator.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_Operator.h
//


extern "C" {


CT_Epetra_Operator_ID_t Epetra_Operator_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Operator_Generalize ( 
  CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id);
}

void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose )
{
    return CTrilinos::tableRepos().get<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->SetUseTranspose(((UseTranspose) != 
        FALSE ? true : false));
}

int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->Apply(*X, *Y);
}

int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CTrilinos::tableRepos().getConst<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = 
        CTrilinos::tableRepos().get<Epetra_MultiVector, 
        CT_Epetra_MultiVector_ID_t>(YID);
    return CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->ApplyInverse(*X, *Y);
}

double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->NormInf();
}

const char * Epetra_Operator_Label ( CT_Epetra_Operator_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->Label();
}

boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->UseTranspose()) ? TRUE : FALSE);
}

boolean Epetra_Operator_HasNormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return ((CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->HasNormInf()) ? TRUE : FALSE);
}

CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Comm, CT_Epetra_Comm_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->Comm() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Map, CT_Epetra_Map_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_Operator, 
        CT_Epetra_Operator_ID_t>(selfID)->OperatorRangeMap() ));
}


} // extern "C"


//
// Definitions from CEpetra_Operator_Cpp.hpp
//


/* get Epetra_Operator from non-const table using CT_Epetra_Operator_ID */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Operator, CT_Epetra_Operator_ID_t>(id);
}

/* get Epetra_Operator from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Operator, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_Operator from either the const or non-const table
 * using CT_Epetra_Operator_ID */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Operator, CT_Epetra_Operator_ID_t>(id);
}

/* get const Epetra_Operator from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Operator, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_Operator in non-const table */
CT_Epetra_Operator_ID_t
CEpetra::storeOperator( Epetra_Operator *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Operator, CT_Epetra_Operator_ID_t>(pobj, false);
}

/* store const Epetra_Operator in const table */
CT_Epetra_Operator_ID_t
CEpetra::storeConstOperator( const Epetra_Operator *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Operator, CT_Epetra_Operator_ID_t>(pobj, false);
}



