
#ifndef EPETRA_ZOLTANQUERY_H
#define EPETRA_ZOLTANQUERY_H

#ifdef HAVE_CONFIG_H
#include <EpetraExt_config.h>
#endif

#include <Zoltan_QueryObject.h>

#include <vector>

class Epetra_CrsGraph;

namespace EpetraExt {

class Epetra_ZoltanQuery : public Zoltan_QueryObject
{

  const Epetra_CrsGraph & graph_;
  const Epetra_CrsGraph * tgraph_;

  std::vector< std::vector<int> > LBProc_;
  std::vector< std::vector<int> > LBProc_Trans_;

  const bool localEdgesOnly_;

 public:

  Epetra_ZoltanQuery( const Epetra_CrsGraph & graph,
                      const Epetra_CrsGraph * tgraph = 0,
                      bool localEdgesOnly = false );

  //General Functions
  int Number_Objects  ( void * data,
                        int * ierr );

  void Object_List    ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_ids,
                        ZOLTAN_ID_PTR local_ids,
                        int weight_dim,
                        float * object_weights,
                        int * ierr );

  int First_Object    ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR first_global_id,
                        ZOLTAN_ID_PTR first_local_id,
                        int weight_dim,
                        float * first_weight,
                        int * ierr );

  int Next_Object     ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        ZOLTAN_ID_PTR next_global_id,
                        ZOLTAN_ID_PTR next_local_id,
                        int weight_dim,
                        float * next_weight,
                        int * ierr );

  int Number_Border_Objects   ( void * data,
                                int number_neighbor_procs,
                                int * ierr );

  void Border_Object_List     ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                int number_neighbor_procs,
                                ZOLTAN_ID_PTR global_ids,
                                ZOLTAN_ID_PTR local_ids,
                                int weight_dim,
                                float * object_weights,
                                int * ierr );
  
  int First_Border_Object     ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                int number_neighbor_procs,
                                ZOLTAN_ID_PTR first_global_id,
                                ZOLTAN_ID_PTR first_local_id,
                                int weight_dim,
                                float * first_weight,
                                int * ierr );
  
  int Next_Border_Object      ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                ZOLTAN_ID_PTR global_id,
                                ZOLTAN_ID_PTR local_id,
                                int number_neighbor_procs,
                                ZOLTAN_ID_PTR next_global_id,
                                ZOLTAN_ID_PTR next_local_id,
                                int weight_dim,
                                float * next_weight,
                                int * ierr );
  
  //Geometry Based Functions
  int Number_Geometry_Objects ( void * data,
                                int * ierr );
  
  void Geometry_Values        ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                ZOLTAN_ID_PTR global_id,
                                ZOLTAN_ID_PTR local_id,
                                double * geometry_vector,
                                int * ierr );
  
  //Graph Based Functions
  int Number_Edges    ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        int * ierr );
  
  void Edge_List      ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        ZOLTAN_ID_PTR neighbor_global_ids,
                        int * neighbor_procs,
                        int weight_dim,
                        float * edge_weights,
                        int * ierr );
  
  //Tree Based Functions
  int Number_Coarse_Objects   ( void * data,
                                int * ierr );

  void Coarse_Object_List     ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                ZOLTAN_ID_PTR global_ids,
                                ZOLTAN_ID_PTR local_ids,
                                int * assigned,
                                int * number_vertices,
                                ZOLTAN_ID_PTR vertices,
                                int * in_order,
                                ZOLTAN_ID_PTR in_vertex,
                                ZOLTAN_ID_PTR out_vertex,
                                int * ierr );

  int First_Coarse_Object     ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                ZOLTAN_ID_PTR first_global_id,
                                ZOLTAN_ID_PTR first_local_id,
                                int * assigned,
                                int * number_vertices,
                                ZOLTAN_ID_PTR vertices,
                                int * in_order,
                                ZOLTAN_ID_PTR in_vertex,
                                ZOLTAN_ID_PTR out_vertex,
                                int * ierr );

  int Next_Coarse_Object      ( void * data,
                                int num_gid_entries,
                                int num_lid_entries,
                                ZOLTAN_ID_PTR global_id,
                                ZOLTAN_ID_PTR local_id,
                                ZOLTAN_ID_PTR next_global_id,
                                ZOLTAN_ID_PTR next_local_id,
                                int * assigned,
                                int * number_vertices,
                                ZOLTAN_ID_PTR vertices,
                                ZOLTAN_ID_PTR in_vertex,
                                ZOLTAN_ID_PTR out_vertex,
                                int * ierr );

  int Number_Children ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        int * ierr );

  void Child_List     ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR parent_global_id,
                        ZOLTAN_ID_PTR parent_local_id,
                        ZOLTAN_ID_PTR child_global_id,
                        ZOLTAN_ID_PTR child_local_id,
                        int * assigned,
                        int * number_vertices,
                        ZOLTAN_ID_PTR vertices,
                        ZOLTAN_REF_TYPE * reference_type,
                        ZOLTAN_ID_PTR in_vertex,
                        ZOLTAN_ID_PTR out_vertex,
                        int * ierr  );

  void Child_Weight   ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        int weight_dim,
                        float * object_weight,
                        int * ierr ); 

  int Number_HG_Edges ( void * data,
    		        int * ierr );

  int Number_HG_Pins ( void * data,
		        int * ierr );

  int HG_Edge_List   ( void * data,
                       int num_gid_entries,
                       int ewgt_dim,
                       int nedge,
                       int maxsize,
                       int * edge_sizes,
                       ZOLTAN_ID_PTR edge_verts,
                       int * edge_procs,
                       float * edge_weights );

};

} //namespace EpetraExt

#endif //EPETRA_ZOLTANQUERY_H

