//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_MigrationFunctions.C$
//
// Purpose        : Static methods which are directly registered with
//		    Zoltan.  They us the static container to access
//		    the dynamic object methods.
//
// Special Notes  : 
//
// Creator        : Robert J. Hoekstra
//
// Creation Date  : 08/04/2000
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-------------------------------------------------------------------------

#include <Zoltan_MigrationFunctions.h>
#include <Zoltan_MigrationContainer.h>
#include <Zoltan_MigrationObject.h>

int Zoltan_MigrationFunctions::Object_Size    (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						int * ierr )
{
  Zoltan_MigrationObject * obj_ptr =
	Zoltan_MigrationContainer::getMigrationObject(
	Zoltan_MigrationContainer::getMigrationID() );

  return obj_ptr->Object_Size( data, num_gid_entries, num_lid_entries,
		global_id, local_id, ierr );
}

void Zoltan_MigrationFunctions::Pre_Migrate   (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						int num_import,
						LB_ID_PTR import_global_ids,
						LB_ID_PTR import_local_ids,
						int * import_procs,
						int num_export,
						LB_ID_PTR export_global_ids,
						LB_ID_PTR export_local_ids,
						int * export_procs,
						int * ierr )
{
  Zoltan_MigrationObject * obj_ptr =
	Zoltan_MigrationContainer::getMigrationObject(
	Zoltan_MigrationContainer::getMigrationID() );

  obj_ptr->Pre_Migrate( data, num_gid_entries, num_lid_entries,
	num_import, import_global_ids, import_local_ids,
	import_procs, num_export, export_global_ids, export_local_ids,
	export_procs, ierr );
}

void Zoltan_MigrationFunctions::Mid_Migrate   (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						int num_import,
						LB_ID_PTR import_global_ids,
						LB_ID_PTR import_local_ids,
						int * import_procs,
						int num_export,
						LB_ID_PTR export_global_ids,
						LB_ID_PTR export_local_ids,
						int * export_procs,
						int * ierr )
{
  Zoltan_MigrationObject * obj_ptr =
	Zoltan_MigrationContainer::getMigrationObject(
	Zoltan_MigrationContainer::getMigrationID() );

  obj_ptr->Mid_Migrate( data, num_gid_entries, num_lid_entries,
	num_import, import_global_ids, import_local_ids,
	import_procs, num_export, export_global_ids, export_local_ids,
	export_procs, ierr );
}

void Zoltan_MigrationFunctions::Post_Migrate  (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						int num_import,
						LB_ID_PTR import_global_ids,
						LB_ID_PTR import_local_ids,
						int * import_procs,
						int num_export,
						LB_ID_PTR export_global_ids,
						LB_ID_PTR export_local_ids,
						int * export_procs,
						int * ierr )
{
  Zoltan_MigrationObject * obj_ptr =
	Zoltan_MigrationContainer::getMigrationObject(
	Zoltan_MigrationContainer::getMigrationID() );

  obj_ptr->Post_Migrate( data, num_gid_entries, num_lid_entries,
	num_import, import_global_ids, import_local_ids,
	import_procs, num_export, export_global_ids, export_local_ids,
	export_procs, ierr );
}

void Zoltan_MigrationFunctions::Pack_Object   (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						int destination_processor,
						int size,
						char * buffer,
						int * ierr )
{
  Zoltan_MigrationObject * obj_ptr =
	Zoltan_MigrationContainer::getMigrationObject(
	Zoltan_MigrationContainer::getMigrationID() );

  obj_ptr->Pack_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, destination_processor, size, buffer, ierr );
}

void Zoltan_MigrationFunctions::Unpack_Object (	void * data,
						int num_gid_entries,
						LB_ID_PTR global_id,
						int size,
						char * buffer,
						int * ierr )
{
  Zoltan_MigrationObject * obj_ptr =
	Zoltan_MigrationContainer::getMigrationObject(
	Zoltan_MigrationContainer::getMigrationID() );

  obj_ptr->Unpack_Object( data, num_gid_entries, global_id, size,
	buffer, ierr );
}
