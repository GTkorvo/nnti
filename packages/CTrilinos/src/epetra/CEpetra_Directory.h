
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


/*! @file CEpetra_Directory.h
 * @brief Wrappers for Epetra_Directory */

/* True C header file! */


#ifndef CEPETRA_DIRECTORY_H
#define CEPETRA_DIRECTORY_H


#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Functions Epetra_Directory_Cast() and Epetra_Directory_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Directory table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/*! @name Explicit casting methods */
/*@{*/

CT_Epetra_Directory_ID_t Epetra_Directory_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Directory_Abstract ( 
  CT_Epetra_Directory_ID_t id );

/*@}*/

/*! @name Epetra_Directory destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Directory::~Epetra_Directory()
*/
void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID );

/*@}*/

/*! @name Epetra_Directory member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual int Epetra_Directory::GetDirectoryEntries( const Epetra_BlockMap& Map, const int NumEntries, const int * GlobalEntries, int * Procs, int * LocalEntries, int * EntrySizes, bool high_rank_sharing_procs=false) const = 0
*/
int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs );

/*! @brief Wrapper for 
   virtual bool Epetra_Directory::GIDsAllUniquelyOwned() const = 0
*/
boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DIRECTORY_H */

