
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


#ifdef HAVE_MPI



/*! @file CEpetra_MpiComm.h
 * @brief Wrappers for Epetra_MpiComm */

/* True C header file! */


#ifndef CEPETRA_MPICOMM_H
#define CEPETRA_MPICOMM_H


#include "mpi.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Distributor.h"
#include "CEpetra_Directory.h"
#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Functions Epetra_MpiComm_Cast() and Epetra_MpiComm_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_MpiComm table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/*! @name Explicit casting methods */
/*@{*/

CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_MpiComm_Abstract ( 
  CT_Epetra_MpiComm_ID_t id );

/*@}*/

/*! @name Epetra_MpiComm constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MpiComm::Epetra_MpiComm(MPI_Comm comm)
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

/*! @brief Wrapper for 
   Epetra_MpiComm::Epetra_MpiComm(const Epetra_MpiComm & Comm)
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID );

/*@}*/

/*! @name Epetra_MpiComm destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_MpiComm::~Epetra_MpiComm()
*/
void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );

/*@}*/

/*! @name Epetra_MpiComm member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Comm * Epetra_MpiComm::Clone() const
*/
CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( 
  CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_MpiComm::Barrier() const
*/
void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(double * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(int * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(long * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(char * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GatherAll(double * MyVals, double * AllVals, int Count) const
*/
int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GatherAll(int * MyVals, int * AllVals, int Count) const
*/
int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GatherAll(long * MyVals, long * AllVals, int Count) const
*/
int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::SumAll(double * PartialSums, double * GlobalSums, int Count) const
*/
int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::SumAll(int * PartialSums, int * GlobalSums, int Count) const
*/
int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::SumAll(long * PartialSums, long * GlobalSums, int Count) const
*/
int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const
*/
int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const
*/
int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const
*/
int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MinAll(double * PartialMins, double * GlobalMins, int Count) const
*/
int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MinAll(int * PartialMins, int * GlobalMins, int Count) const
*/
int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MinAll(long * PartialMins, long * GlobalMins, int Count) const
*/
int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::ScanSum(double * MyVals, double * ScanSums, int Count) const
*/
int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::ScanSum(int * MyVals, int * ScanSums, int Count) const
*/
int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::ScanSum(long * MyVals, long * ScanSums, int Count) const
*/
int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   MPI_Comm Epetra_MpiComm::Comm() const
*/
MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MyPID() const
*/
int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::NumProc() const
*/
int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Distributor * Epetra_MpiComm::CreateDistributor() const
*/
CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( 
  CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Directory * Epetra_MpiComm::CreateDirectory(const Epetra_BlockMap & Map) const
*/
CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GetMpiTag() const
*/
int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   MPI_Comm Epetra_MpiComm::GetMpiComm() const
*/
MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID );

/*@}*/

/*! @name Epetra_MpiComm operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MpiComm & Epetra_MpiComm::operator=(const Epetra_MpiComm & Comm)
*/
void Epetra_MpiComm_Assign ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MPICOMM_H */

#endif /* HAVE_MPI */


