/*! \@HEADER */
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
Questions? Contact M. Nicole Lemaster (mnlemas\@sandia.gov)

************************************************************************
*/
/*! \@HEADER */


#include "CTrilinos_config.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Vector.h"
#include "CEpetra_Export.h"
#include "CEpetra_Import.h"
#include "CEpetra_OffsetIndex.h"
#include "CEpetra_SrcDistObject.h"
#include "Epetra_CombineMode.h"
#include "Epetra_DistObject.h"
#include "CEpetra_DistObject.h"
#include "CEpetra_DistObject_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_DistObject_ID_t Epetra_DistObject_Cast ( 
  CTrilinos_Universal_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_DistObject , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID))));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Create(MapID, zeroOut));

  ECHO(CT_Epetra_DistObject_ID_t selfID = Epetra_DistObject_Cast(Epetra_Vector_Abstract(vecID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(vecID, selfID), true);
  TEST_EQUALITY(selfID.type, CT_Epetra_DistObject_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);

  /* These casts should be valid */
  ECHO(CT_Epetra_DistObject_ID_t dupID = Epetra_DistObject_Cast(Epetra_DistObject_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), true);
  ECHO(CT_Epetra_SrcDistObject_ID_t sdoID = Epetra_SrcDistObject_Cast(Epetra_DistObject_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, sdoID), true);
  ECHO(CT_Epetra_DistObject_ID_t dupID2 = Epetra_DistObject_Cast(Epetra_SrcDistObject_Abstract(sdoID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(dupID2, sdoID), true);
}

/**********************************************************************
void Epetra_DistObject_Destroy ( 
  CT_Epetra_DistObject_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_DistObject , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID))));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Create(MapID, zeroOut));
  ECHO(CT_Epetra_DistObject_ID_t selfID = Epetra_DistObject_Cast(Epetra_Vector_Abstract(vecID)))

  ECHO(Epetra_DistObject_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );
 **********************************************************************/

/**********************************************************************
int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );
 **********************************************************************/

/**********************************************************************
int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );
 **********************************************************************/

/**********************************************************************
int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_DistObject , Map )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID))));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Create(MapID, zeroOut));
  ECHO(CT_Epetra_DistObject_ID_t selfID = Epetra_DistObject_Cast(Epetra_Vector_Abstract(vecID)));

  ECHO(CT_Epetra_BlockMap_ID_t mapID = Epetra_DistObject_Map(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int els = Epetra_BlockMap_NumGlobalElements(mapID));
  TEST_EQUALITY(els, NumGlobalElements);
}

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_DistObject , Comm )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID))));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Create(MapID, zeroOut));
  ECHO(CT_Epetra_DistObject_ID_t selfID = Epetra_DistObject_Cast(Epetra_Vector_Abstract(vecID)));

  ECHO(CT_Epetra_Comm_ID_t cID = Epetra_DistObject_Comm(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int np = Epetra_Comm_NumProc(cID));
  TEST_EQUALITY(NumProc, np);
}

/**********************************************************************
boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_DistObject , DistributedGlobal )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(
       Epetra_Map_Create(NumGlobalElements, IndexBase, CommID))));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t vecID = Epetra_Vector_Create(MapID, zeroOut));
  ECHO(CT_Epetra_DistObject_ID_t selfID = Epetra_DistObject_Cast(Epetra_Vector_Abstract(vecID)));

  ECHO(boolean dg = Epetra_DistObject_DistributedGlobal(selfID));

  /* Now check the result of the call to the wrapper function */
  boolean isdg = (NumProc > 1 ? TRUE : FALSE);
  TEST_EQUALITY_CONST(dg, isdg);
}


/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

