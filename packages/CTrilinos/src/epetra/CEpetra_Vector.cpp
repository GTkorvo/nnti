
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

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Vector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_Vector.h
//


extern "C" {


CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut )
{
    return CTrilinos::tableRepos().store<Epetra_Vector, 
        CT_Epetra_Vector_ID_t>(new Epetra_Vector(
        *CEpetra::getConstBlockMap(MapID), ((
        zeroOut) != FALSE ? true : false)));
}

CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID )
{
    return CTrilinos::tableRepos().store<Epetra_Vector, 
        CT_Epetra_Vector_ID_t>(new Epetra_Vector(
        *CEpetra::getConstVector(SourceID)));
}

CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * V )
{
    return CTrilinos::tableRepos().store<Epetra_Vector, 
        CT_Epetra_Vector_ID_t>(new Epetra_Vector(
        (Epetra_DataAccess) CV, *CEpetra::getConstBlockMap(MapID), 
        V));
}

CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index )
{
    return CTrilinos::tableRepos().store<Epetra_Vector, 
        CT_Epetra_Vector_ID_t>(new Epetra_Vector(
        (Epetra_DataAccess) CV, *CEpetra::getConstMultiVector(
        SourceID), Index));
}

void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceGlobalValues(
        NumEntries, Values, Indices);
}

int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceMyValues(NumEntries, 
        Values, Indices);
}

int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoGlobalValues(
        NumEntries, Values, Indices);
}

int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoMyValues(NumEntries, 
        Values, Indices);
}

int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceGlobalValues(
        NumEntries, BlockOffset, Values, Indices);
}

int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceMyValues(NumEntries, 
        BlockOffset, Values, Indices);
}

int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoGlobalValues(
        NumEntries, BlockOffset, Values, Indices);
}

int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoMyValues(NumEntries, 
        BlockOffset, Values, Indices);
}

int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V )
{
    return CEpetra::getConstVector(selfID)->ExtractCopy(V);
}

int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V )
{
    return CEpetra::getConstVector(selfID)->ExtractView(V);
}

double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index )
{
    const Epetra_Vector& self = *( CEpetra::getConstVector(selfID) );

    return self[index];
}


} // extern "C"


//
// Definitions from CEpetra_Vector_Cpp.hpp
//


/* get Epetra_Vector from non-const table using CT_Epetra_Vector_ID */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Vector, CT_Epetra_Vector_ID_t>(id);
}

/* get Epetra_Vector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Vector, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_Vector from either the const or non-const table
 * using CT_Epetra_Vector_ID */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Vector, CT_Epetra_Vector_ID_t>(id);
}

/* get const Epetra_Vector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Vector, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_Vector in non-const table */
CT_Epetra_Vector_ID_t
CEpetra::storeVector( Epetra_Vector *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Vector, CT_Epetra_Vector_ID_t>(pobj, false);
}

/* store const Epetra_Vector in const table */
CT_Epetra_Vector_ID_t
CEpetra::storeConstVector( const Epetra_Vector *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Vector, CT_Epetra_Vector_ID_t>(pobj, false);
}



