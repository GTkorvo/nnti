
#include <stk_rebalance/ZoltanPartition.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

/* The Zoltan Include file has an odd name: lbi_ */
#include <lbi_const.h>

#include <Teuchos_ParameterList.hpp>

/* sstream to convert Zoltan double version number
   to string.
*/
#include <sstream>

using namespace std;
using namespace stk;
using namespace stk::rebalance;

#define STK_GEOMDECOMP_DEBUG 0

namespace {


double static_zoltan_version(const double v=0) {
  static double version=0;
  if (v) { version=v;}
  return version;
}

inline unsigned num_gid_entries() {
  static const unsigned n=2;
  return n;
}

inline unsigned num_lid_entries() {
  static const unsigned n=2;
  return n;
}
inline unsigned wdim() {
  static const unsigned n=1;
  return n;
}

//inline void convert_param_to_string(const Parameters &from,
//				    vector < pair<std::string, std::string> > &to)
//{
//  Parameters::const_iterator
//    from_iter  = from.begin(),
//    from_end   = from.end();
//
//  for (; from_iter != from_end; ++from_iter) {
//    to.push_back(std::pair<std::string,std::string>((*from_iter).name(),
//					  (*from_iter).get()));
//  }
//}
//
//inline void fill_parameters (const char *T[][2],
//			     const int  i,
//			     Parameters &Entry)
//{
//  for (int j=0; j<i; ++j) Entry.set(T[j][0], T[j][1]);
//}
//
//void fill_name_conversion (Parameters &Name_Conversion)
//{
//  const char *General[][2] =
//    {
//      { "LOAD BALANCING METHOD"      , "LB_METHOD"  },
//      { "ZOLTAN DEBUG LEVEL"         , "DEBUG_LEVEL" },
//      { "DEBUG PROCESSOR NUMBER"     , "DEBUG_PROCESSOR" },
//      { "TIMER"                      , "TIMER" },
//      { "DETERMINISTIC DECOMPOSITION", "DETERMINISTIC" },
//      { "DEBUG MEMORY"               , "DEBUG_MEMORY" },
//      { "IMBALANCE TOLERANCE"        , "IMBALANCE_TOL" },
//      { "RENUMBER PARTITIONS"        , "REMAP" },
//      { "KEEP CUTS"                  , "KEEP_CUTS" },
//      { "REUSE CUTS"                 , "RCB_REUSE" },
//      { "RCB RECOMPUTE BOX"          , "RCB_RECOMPUTE_BOX" },
//      { "CHECK GEOMETRY"             , "CHECK_GEOM" },
//      { "LOCK RCB DIRECTIONS"        , "RCB_LOCK_DIRECTIONS" },
//      { "SET RCB DIRECTIONS"         , "RCB_SET_DIRECTIONS" },
//      { "RCB MAX ASPECT RATIO"       , "RCB_MAX_ASPECT_RATIO" },
//      { "RECTILINEAR RCB BLOCKS"     , "RCB_RECTILINEAR_BLOCKS" },
//      { "OCTREE DIMENSION"           , "OCT_DIM" },
//      { "OCTREE METHOD"              , "OCT_METHOD" },
//      { "OCTREE MIN OBJECTS"         , "OCT_MINOBJECTS" },
//      { "OCTREE MAX OBJECTS"         , "OCT_MAXOBJECTS" },
//      // These values are never changed, but must
//      // be set so default values work correctly.
//      { "NUMBER GLOBAL ID ENTRIES"   , "NUM_GID_ENTRIES" },
//      { "NUMBER LOCAL ID ENTRIES"    , "NUM_LID_ENTRIES" },
//      { "OBJECT WEIGHT DIMENSION"    , "OBJ_WEIGHT_DIM" },
//      { "RETURN LISTS"               , "RETURN_LISTS" },
//      { "AUTOMATIC MIGRATION"        , "AUTO_MIGRATE" },
//      { "DISTANCE"                   , "DISTANCE" }
//    };
//  const char *RCB[][2] =
//    {
//      { "OVER ALLOCATE MEMORY"       , "RCB_OVERALLOC" },
//      { "ALGORITHM DEBUG LEVEL"      , "RCB_OUTPUT_LEVEL" }
//    };
//  const char *RIB[][2] =
//    {
//      { "OVER ALLOCATE MEMORY"       , "RIB_OVERALLOC" },
//      { "ALGORITHM DEBUG LEVEL"      , "RIB_OUTPUT_LEVEL" }
//    };
//  const char *HSFC[][2] =
//    {
//      { "OVER ALLOCATE MEMORY"       , "" },
//      { "ALGORITHM DEBUG LEVEL"      , "" }
//    };
//  const char *OCT[][2] =
//    {
//      { "OVER ALLOCATE MEMORY"       , "" },
//      { "ALGORITHM DEBUG LEVEL"      , "OCT_OUTPUT_LEVEL" }
//    };
//
//  const int  Table_lens[] = {sizeof(General)/(2*sizeof(char *)),
//			     sizeof(RCB)    /(2*sizeof(char *)),
//			     sizeof(RIB)    /(2*sizeof(char *)),
//			     sizeof(HSFC)   /(2*sizeof(char *)),
//			     sizeof(OCT)    /(2*sizeof(char *))};
//
//  const char *Table_names[] = {"General",
//			       "0",
//			       "1",
//			       "2",
//			       "3"};
//
//  ThrowAssert (sizeof(Table_names)/sizeof(char *) ==
//	  sizeof(Table_lens)/sizeof(int));
//
//  fill_parameters (General,
//		   Table_lens[0],
//		   Name_Conversion.set_nested(Table_names[0]));
//  fill_parameters (RCB,
//		   Table_lens[1],
//		   Name_Conversion.set_nested(Table_names[1]));
//  fill_parameters (RIB,
//		   Table_lens[2],
//		   Name_Conversion.set_nested(Table_names[2]));
//  fill_parameters (HSFC,
//		   Table_lens[3],
//		   Name_Conversion.set_nested(Table_names[3]));
//  fill_parameters (OCT,
//		   Table_lens[4],
//		   Name_Conversion.set_nested(Table_names[4]));
//}
//
//
//void fill_value_conversion (Parameters &Value_Conversion)
//{
//  const char *LB_METHOD[][2] =
//    {
//      { "0"   , "RCB"  },
//      { "1"   , "RIB" },
//      { "2"   , "HSFC" },
//      { "3"   , "OCTPART" },
//    };
//  const char *TIMER[][2] =
//    {
//      { "0"   , "WALL"  },
//      { "1"   , "CPU" },
//    };
//
//  const int   Table_lens[]  = {sizeof(LB_METHOD)/(2*sizeof(char *)),
//			       sizeof(TIMER)    /(2*sizeof(char *))};
//  const char *Table_names[] = {"LOAD BALANCING METHOD",
//			       "TIMER"};
//
//  ThrowAssert (sizeof(Table_names)/sizeof(char *) ==
//	  sizeof(Table_lens)/sizeof(int));
//
//  fill_parameters (LB_METHOD,
//		   Table_lens[0],
//		   Value_Conversion.set_nested(Table_names[0]));
//  fill_parameters (TIMER,
//		   Table_lens[1],
//		   Value_Conversion.set_nested(Table_names[1]));
//}
//
//void fill_default_value (Parameters &Default_Value)
//{
//  const char *General[][2] =
//    {
//      { "LOAD BALANCING METHOD"      , "0"  },
//      // NOTE: "LOAD BALANCING METHOD" default
//      // Is also hard coded in convert_names_and_values().
//      { "RENUMBER PARTITIONS"        , "1" },
//      { "ZOLTAN DEBUG LEVEL"         , "0" },
//      { "TIMER"                      , "0" },
//      { "DETERMINISTIC DECOMPOSITION", "1" },
//      { "DEBUG MEMORY"               , "1" },
//      { "IMBALANCE TOLERANCE"        , "1.1" },
//      { "KEEP CUTS"                  , "1" },
//      { "REUSE CUTS"                 , "1" },
//      //      { "RCB RECOMPUTE BOX"          , "0" },
//      { "OVER ALLOCATE MEMORY"       , "1.1" },
//      { "ALGORITHM DEBUG LEVEL"      , "0" },
//      { "OCTREE MIN OBJECTS"         , "1" },
//      { "OCTREE MAX OBJECTS"         , "1" },
//      // These values are never changed, but must
//      // be set so default values work correctly.
//      { "NUMBER GLOBAL ID ENTRIES"   , "2" },
//      { "NUMBER LOCAL ID ENTRIES"    , "2" },
//      { "OBJECT WEIGHT DIMENSION"    , "1" },
//      { "RETURN LISTS"               , "EXPORT" }
//    };
//
//  const int   Table_lens[]  = {sizeof(General)/(2*sizeof(char *))};
//  const char *Table_names[] = {"General"};
//
//  ThrowAssert (sizeof(Table_names)/sizeof(char *) ==
//	  sizeof(Table_lens)/sizeof(int));
//
//  fill_parameters (General,
//		   Table_lens[0],
//		   Default_Value.set_nested(Table_names[0]));
//}
//
//
//
//#if STK_GEOMDECOMP_DEBUG>=2
//void debug_print_decomp_export(Zoltan   *zoltan,
//			       Zoltan *zoltan_id )
//{
//  int i;
//  int           num_export   = zoltan->Num_Exported();
//  ZOLTAN_ID_PTR export_lids  = zoltan->Export_Local_IDs();
//  ZOLTAN_ID_PTR export_gids  = zoltan->Export_Global_IDs();
//  int*          export_procs = zoltan->Export_Proc_IDs();
//  int           Z_LID_SIZE   = zoltan->Num_Lid_Entries();
//  int           Z_GID_SIZE   = zoltan->Num_Gid_Entries();
//
//  if ( export_gids!=NULL && export_lids!=NULL && export_procs!=NULL ) {
//    Env::output() << ": Zoltan RCB EXPORTS" << std::endl;
//    for ( i = 0; i < num_export; i++ ) {
//      Env::output()
//	<< "  " << i
//	<< ":  GID = "
//	<< "T" << zoltan_id->Type( &export_gids[i*Z_GID_SIZE]) << "  "
//	<< "I" << zoltan_id->Index(&export_gids[i*Z_GID_SIZE]) << "  "
//	<< "P" << zoltan_id->Proc( &export_gids[i*Z_GID_SIZE]) << "  "
//	<< "    LID = "
//	<< "T" << zoltan_id->Type( &export_lids[i*Z_LID_SIZE]) << "  "
//	<< "I" << zoltan_id->Index(&export_lids[i*Z_LID_SIZE]) << "  "
//	<< "  EXPORT_PROC_ID = "
//	<< export_procs[i]
//	<< std::endl;
//    }
//    Env::output_flush();
//  }
//}
//#endif



extern "C" {
  int Callback_Num_Elements( void *data, int *ierr );
  void Callback_Element_List( void *data,
			     int Num_gid_entries,
			     int Num_lid_entries,
			     ZOLTAN_ID_PTR global_ids,
			     ZOLTAN_ID_PTR local_ids,     //{Iterator or index}
			     int wdim,
			     float *weights,
			     int *ierr );
  int Callback_First_Element( void *data,
			      int Num_gid_entries,
			      int Num_lid_entries,
			      ZOLTAN_ID_PTR global_id,
			      ZOLTAN_ID_PTR local_id,   //{Iterator or index}
			      int wdim,
			      float *weight,
			      int *ierr );
  int Callback_Next_Element( void *data,
			     int Num_gid_entries,
			     int Num_lid_entries,
			     ZOLTAN_ID_PTR global_id,
			     ZOLTAN_ID_PTR local_id,     //{Iterator or index}
			     ZOLTAN_ID_PTR next_global_id,
			     ZOLTAN_ID_PTR next_local_id,
			     int wdim,
			     float *next_weight,
			     int *ierr );

  int Callback_Num_Dimensions( void *data, int *ierr );
  void Callback_Centroid_Coord( void *data,
				int Num_gid_entries,
				int Num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,  //{Iterator or index}
				double *geom,
				int *ierr );
  int Callback_Num_Edges( void *data,
                          int Num_gid_entries,
                          int Num_lid_entries,
                          ZOLTAN_ID_PTR global_id,
                          ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                          int *ierr );
  void Callback_Edge_List( void *data,
                           int Num_gid_entries,
                           int Num_lid_entries,
                           ZOLTAN_ID_PTR global_id,
                           ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                           ZOLTAN_ID_PTR nbor_global_id,
                           int *nbor_procs,
                           int wgt_dim,
                           float *ewgts,
                           int *ierr );
}


int Callback_Num_Elements( void *data, int *ierr )
{
  if ( data == NULL ) {
    *ierr = ZOLTAN_FATAL;  // Set FATAL Zoltan error flag
    return 0;
  }
  stk::rebalance::GeomDecomp   *gdata =  static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan       *zdata = dynamic_cast<stk::rebalance::Zoltan*>(gdata);

  if (zdata == NULL ) {
    *ierr = ZOLTAN_FATAL;  // Set FATAL Zoltan error flag
    return 0;
  }

  int ne = zdata->num_elems();

  //ThrowAssert( ne > -1 );

  *ierr = ZOLTAN_OK;
  return ne;
}

void Callback_Element_List( void *data,
			    int Num_gid_entries,
			    int Num_lid_entries,
			    ZOLTAN_ID_PTR global_ids,
			    ZOLTAN_ID_PTR local_ids,     //{Iterator or index}
			    int weightdim,
			    float *weights,
			    int *ierr )
{
  //ThrowRequire( num_gid_entries() == static_cast<unsigned>(Num_gid_entries) );
  //ThrowRequire( num_lid_entries() == static_cast<unsigned>(Num_lid_entries) );

  if (!data) {
    *ierr = ZOLTAN_FATAL;           // Set FATAL Zoltan error flag
    return;
  }

  stk::rebalance::GeomDecomp   *gdata =  static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan       *zdata = dynamic_cast<stk::rebalance::Zoltan*>     (gdata);

  if (!zdata) {
    *ierr = ZOLTAN_FATAL;           // Set FATAL Zoltan error flag
    return;
  }

  unsigned k=0, l=0;
  const unsigned num_local_ids = gdata->num_moid();
  for (unsigned j=0; j<num_local_ids; ++j) {
    local_ids [k] = j;
    global_ids[k] = 0; // region_id
    ++k;
    local_ids [k] = 0;
    global_ids[k] = static_cast<ZOLTAN_ID_TYPE>(gdata->globalID(j));
    ++k;
    if (weightdim) {
      weights   [l++] = gdata->object_weight(j);
    } else {
      ++l;
    }
  }
  //ThrowRequire(2*zdata->num_elems() == k && 2*l == k);

  *ierr = ZOLTAN_OK;
  return;

}


int Callback_First_Element( void *data,
			int Num_gid_entries,
			int Num_lid_entries,
			ZOLTAN_ID_PTR global_id,
			ZOLTAN_ID_PTR local_id,   //{Iterator or index}
			int wdim,
			float *weight,
			int *ierr )
{

  //ThrowAssert( 2 == Num_gid_entries );
  //ThrowAssert( 2 == Num_lid_entries );

  if ( data == NULL ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  stk::rebalance::GeomDecomp   *gdata =  static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan *zdata       = dynamic_cast<stk::rebalance::Zoltan*>     (gdata);

  if ( zdata == NULL ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  zdata->iter_init();

  //: Set first element
  local_id [ 0 ] = static_cast<ZOLTAN_ID_TYPE>(zdata->iter_current_key());
  global_id[ 0 ] = 0;
  global_id[ 1 ] =
    static_cast<ZOLTAN_ID_TYPE>(zdata->iter_mesh_object()->identifier());
  //: Set weight for first element
  weight[ 0 ] = zdata->iter_object_weight();
  for (int j=1; j < wdim ; j++) {
    weight[ j ] = weight[0];
  }

  *ierr = ZOLTAN_OK;
  return 1;

}

int Callback_Next_Element( void *data,
				  int Num_gid_entries,
				  int Num_lid_entries,
				  ZOLTAN_ID_PTR global_id,
				  ZOLTAN_ID_PTR local_id,     //{Iterator or index}
				  ZOLTAN_ID_PTR next_global_id,
				  ZOLTAN_ID_PTR next_local_id,
				  int wdim,
				  float *next_weight,
				  int *ierr )
{
  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  //ThrowAssert( 2 == Num_gid_entries ); // global_id()'s Global meshObj ID index
  //ThrowAssert( 2 == Num_lid_entries ); // local_id()'s  Local element ID index

  if (!data) {
    *ierr = ZOLTAN_FATAL;           // Set FATAL Zoltan error flag
    return 0;
  }

  stk::rebalance::GeomDecomp   *gdata =  static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan *zdata       = dynamic_cast<stk::rebalance::Zoltan*>     (gdata);

  if (!zdata) {
    *ierr = ZOLTAN_FATAL;           // Set FATAL Zoltan error flag
    return 0;
  }

  // Check that we are in sync with Zoltan.
  unsigned key = zdata->iter_current_key();
  //ThrowAssert (local_id [0] == key);
  //ThrowAssert (global_id[0] == 0);
  //ThrowAssert (global_id[1] == (unsigned) zdata->iter_mesh_object()->global_id());

  // Increment local id in current region
  ++(*zdata);

  // Store region, local, and global ids in the "next" arrays
  key = zdata->iter_current_key();
  next_local_id [ 0 ] = key;
  next_global_id[ 0 ] = 0;
  next_global_id[ 1 ] =
    static_cast<ZOLTAN_ID_TYPE>(zdata->iter_mesh_object()->identifier());
  //: Set weight for next element
  next_weight[ 0 ] = zdata->iter_object_weight();
  for (int j=1; j < wdim ; j++) {
    next_weight[ j ] = next_weight[0];
  }

  *ierr = ZOLTAN_OK;
  return 1;

}

int Callback_Num_Dimensions( void *data, int *ierr )
{
  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  throw std::runtime_error("Zoltan: Callback_Num_Dimensions needs to be implemented.");
  //// Replace with possibly TopologicalMetaData::spatial_dimension
  //int  nd = zdata->object_coord_ref()->length_deprecated(0);

  //ThrowAssert( nd > 0 && nd < 4 );

  *ierr = ZOLTAN_OK;
  return 0 /* nd */;

}

void Callback_Centroid_Coord( void *data,
				     int Num_gid_entries,
				     int Num_lid_entries,
				     ZOLTAN_ID_PTR global_id,
				     ZOLTAN_ID_PTR local_id,  //{Iterator or index}
				     double *geom,
				     int *ierr )
{
  //int i = 0;

  std::vector<double> temp(3,0.0);

  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  //ThrowAssert( 2 == Num_gid_entries ); // global_id()'s Region ID index
  //ThrowAssert( 2 == Num_lid_entries ); // local_id()'s  Local element ID index

  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }

  //int lid = local_id[  0 ]; // Local Element ID

  throw std::runtime_error("Zoltan: Callback_Centroid_Coord needs to be implemented.");
  //const Fmwk::MeshObj & target_obj = * zdata->mesh_object( lid );
  //const mesh::Field              & coor = * zdata->object_coord_ref();
  //int                        nd   =   coor.length_deprecated(0);

  /*
   * Obtain the centroid coordinates of the element by averaging all
   * the nodal coordinates of the nodes associated with the element.
   * Use GeomDecomp friend function, obj_to_point( , , )
   *
   */

  //rebalance::GeomDecomp::obj_to_point( target_obj, coor, temp );

  //for ( ; i < nd ; i++ ) geom[ i ] = (double) temp[ i ];

  *ierr = ZOLTAN_OK;
}



//void getNeighbors( const Fmwk::MeshObj & obj,
//                   std::set<const Fmwk::MeshObj*> & nodes ) {
//
//  nodes.clear();
//
//  const Fmwk::MeshObj::const_iterator
//    iElemEnd = obj.end_relation(Fmwk::MeshObj::ELEMENT,
//                                Fmwk::Relation::USED_BY);
//  Fmwk::MeshObj::const_iterator
//    iElem = obj.begin_relation(Fmwk::MeshObj::ELEMENT,
//                               Fmwk::Relation::USED_BY);
//
//  for ( ; iElem != iElemEnd; ++iElem ) {
//    Fmwk::MeshObj * elem = iElem->obj;
//
//    const Fmwk::MeshObj::const_iterator
//      iNodeEnd = elem->end_relation(Fmwk::MeshObj::NODE,
//                                    Fmwk::Relation::USES);
//    Fmwk::MeshObj::const_iterator
//      iNode = elem->begin_relation(Fmwk::MeshObj::NODE,
//                                   Fmwk::Relation::USES);
//
//    for ( ; iNode != iNodeEnd; ++iNode ) {
//      Fmwk::MeshObj * node = iNode->obj;
//      if (&obj != node) {
//        nodes.insert( node );
//      }
//    }
//  }
//}

//int numEdges( const Fmwk::MeshObj & obj ) {
//
//  std::set<const Fmwk::MeshObj*> nodes;
//
//  getNeighbors( obj, nodes );
//  return nodes.size();
//}

int Callback_Num_Edges( void *data,
                        int Num_gid_entries,
                        int Num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                        int *ierr )
{
  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  //ThrowAssert( 2 == Num_gid_entries ); // global_id()'s Region ID index
  //ThrowAssert( 2 == Num_lid_entries ); // local_id()'s  Local node ID index

  *ierr = ZOLTAN_OK;

  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  //int lid = local_id[  0 ]; // Local Element ID

  throw std::runtime_error("Zoltan: Callback_Num_Edges needs to be implemented.");
  //const Fmwk::MeshObj & target_obj = * zdata->mesh_object( lid );

  return 0 /* numEdges( target_obj ) */;

}

void Callback_Edge_List( void *data,
                         int Num_gid_entries,
                         int Num_lid_entries,
                         ZOLTAN_ID_PTR global_id,
                         ZOLTAN_ID_PTR local_id,  //{Iterator or index}

                         // Returned
                         ZOLTAN_ID_PTR nbor_global_id, //owning processor
                         int *nbor_procs,
                         int wgt_dim,
                         float *ewgts,

                         int *ierr )
{
  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  //ThrowAssert( 2 == Num_gid_entries ); // global_id()'s Region ID index
  //ThrowAssert( 2 == Num_lid_entries ); // local_id()'s  Local node ID index

  *ierr = ZOLTAN_OK;

  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }

  //int lid = local_id[  0 ]; // Local Node ID

  throw std::runtime_error("Zoltan: Callback_Edge_List needs to be implemented.");
  //const Fmwk::MeshObj & target_obj = * zdata->mesh_object( lid );

  //std::set<const Fmwk::MeshObj*> nodes;
  //getNeighbors( target_obj, nodes );

  //int counter(0);
  //for ( std::set<const Fmwk::MeshObj*>::iterator i = nodes.begin();
  //      i != nodes.end(); ++i ) {
  //  nbor_global_id[counter*2+0] = 0; // region_id

  //  const Fmwk::MeshObj & n = **i;
  //  nbor_global_id[counter*2+1] = n.global_id();

  //  nbor_procs[counter] = n.owner_processor_rank();

  //  if ( wgt_dim ) {
  //    ewgts[counter] = 1;
  //  }
  //  ++counter;
  //}
}

}

static Teuchos::ParameterList *Name_Conversion =NULL;
static Teuchos::ParameterList *Value_Conversion=NULL;

double rebalance::Zoltan::zoltan_version()   const { return static_zoltan_version();  }

//: ===========
//: Constructor
//: ===========

//Diag::Writer &
//rebalance::Zoltan::verbose_print(
//  Diag::Writer &                dout) const
//{
//  if (dout.shouldPrint()) {
//    dout << "Fmwk::Zoltan" << sierra::Diag::push << dendl;
//    GeomDecomp::verbose_print(dout).dendl();
//    dout.m(LOG_MEMBERS) << "parameter_entry_Name, " << parameter_entry_Name << dendl;
//    dout.m(LOG_MEMBERS) << "zoltan_version, " << static_zoltan_version() << dendl;
//
//    dout << sierra::Diag::pop;
//  }
//
//  return dout;
//}

//namespace {
//void Merge_Parameters(std::vector <std::pair<std::string, std::string> > &Str_Zoltan_Params,
//		      const Teuchos::ParameterList &Zoltan_Params) {
//  Fmwk::Teuchos::ParameterList Merged_Zoltan_Params   ;
//  Fmwk::Teuchos::ParameterList Converted_Zoltan_Params;
//
//  rebalance::Zoltan::merge_default_values (Zoltan_Params,
//				      Merged_Zoltan_Params);
//
//  rebalance::Zoltan::convert_names_and_values(Merged_Zoltan_Params,
//					 Converted_Zoltan_Params);
//
//  convert_param_to_string (Converted_Zoltan_Params,
//			   Str_Zoltan_Params);
//  return;
//}
//}

rebalance::Zoltan rebalance::Zoltan::create_default(ParallelMachine pm, const Teuchos::ParameterList & rebal_region_parameters)
{
  const std::string &param_name = stk::rebalance::GeomDecomp::zoltan_parameters_name();
  const bool exist = rebal_region_parameters.isParameter(param_name);
  const std::string parameter_set = ( exist ? (rebal_region_parameters.get<const std::string>(param_name)) : "");

  return rebalance::Zoltan(pm, parameter_set);
}

rebalance::Zoltan::Zoltan(ParallelMachine pm, const std::string &Parameters_Name) :
  GeomDecomp(pm),
  zoltan_id(NULL)
{
  /* Get the default Zoltan parameter set from the Domain.
     This will be created the first time and returned every
     time after that.
  */
  //Teuchos::ParameterList &domain_parameters =
  //  Fmwk::Domain::singleton()->parameters().set_nested(GeomDecomp::zoltan_parameters_name());

  /* Determine if the default set of parameters already exists. */
  //if ( !domain_parameters.get_nested(GeomDecomp::default_parameters_name()) ) {
  //  init_default_parameters();
  //}

  /* If name is empty, use default values */
  std::string Default_Name =
    (Parameters_Name.empty()) ? default_parameters_name() : Parameters_Name ;

  //const Teuchos::ParameterList *Zoltan_Params = domain_parameters.get_nested(Default_Name);
  //if ( !Zoltan_Params ) {
  //  throw RuntimeError() << "The Zoltan parameter set '" << Default_Name << "' does not exist." << std::endl << StackTrace;
  //}

  /* Save this library name for future reference. */
  //parameter_entry_Name = Zoltan_Params->nested_name();
  //ThrowAssert (!parameter_entry_Name.empty());

  std::vector <std::pair<std::string, std::string> > Str_Zoltan_Params;
  //Merge_Parameters(Str_Zoltan_Params, *Zoltan_Params);

  init(Str_Zoltan_Params);
}


void rebalance::Zoltan::init( const vector< pair<std::string,std::string> >
		   &dynamicLoadRebalancingParameters ) {
  if (0==static_zoltan_version()) {
    const double v = init_zoltan_library();
    static_zoltan_version(v);
  }

  /**
   * Create a zoltanID object
   */

  zoltan_id = Zoltan_Create( comm_ );
  //if ( zoltan_id == NULL ) {
  //  throw RuntimeError() << " P" << Env::parallel_rank()
  //      		 << ": (FATAL ERROR) Zoltan_Create() returned NULL" << std::endl
  //      		 << StackTrace;
  //}

  /**
   * Set up dynamic load rebalancing
   */

  //ThrowAssert (!dynamicLoadRebalancingParameters.empty());

  vector<pair<std::string,std::string> >::const_iterator
    P  =  dynamicLoadRebalancingParameters.begin(),
    PE =  dynamicLoadRebalancingParameters.end();

  for ( ; PE != P ; P++ ) {

    //char * label = const_cast<char*>( P->first.c_str() ) ;
    //char * value = const_cast<char*>( P->second.c_str() ) ;

    //if (ZOLTAN_OK != (Zoltan_Set_Param(zoltan_id,label,value)))
    //{
    //  throw RuntimeError() << " P" << Env::parallel_rank()
    //    		   << ": FATAL ERROR in Zoltan_Set_Param when setting parameter "
    //    		   << label << " to value '" << value << "'" << std::endl << StackTrace;
    //}
    //if (case_strcmp(label, "NUM_GID_ENTRIES") == 0)
    //  ThrowAssert (::num_gid_entries() == (unsigned) atoi(value) );
    //if (case_strcmp(label, "NUM_LID_ENTRIES") == 0)
    //  ThrowAssert (::num_lid_entries() == (unsigned) atoi(value) );
    //if (case_strcmp(label, "OBJ_WEIGHT_DIM") == 0)
    //  ThrowRequire (::wdim() == (unsigned) atoi(value) );
  }

  /**
   * Register the Zoltan/SIERRA "call-back" (querry) functions.
   */
  //if ( ZOLTAN_OK != register_callbacks() )
  //  throw RuntimeError() << "zoltan->Register_Callbacks, error code "<< ZOLTAN_FATAL << std::endl << StackTrace;

#if STK_GEOMDECOMP_DEBUG>=2
  {
    debug_print_decomp_export( zoltan, zoltan_id );
  }
#endif
  return;

}

double rebalance::Zoltan::init_zoltan_library () {
  float version = 0.0;
  //if ( Zoltan_Initialize( 0, NULL, &version) != ZOLTAN_OK )
  //  throw std::runtime_error() << "Zoltan_Initialize() != ZOLTAN_OK " << std::endl << StackTrace;

  static_zoltan_version(version);
  std::ostringstream s;
  s << version;

  //sierra::ProductRegistry::instance().addTPL("Zoltan", s.str());

  return version;
}



//: ==========
//: Destructor
//: ==========

rebalance::Zoltan::~Zoltan()
{
  if ( zoltan_id != NULL ) {
    Zoltan_Destroy( &zoltan_id );
  }
  zoltan_id = NULL ;
  if(Name_Conversion) {
    delete Name_Conversion;
    Name_Conversion = NULL;
  }
  if(Value_Conversion) {
    delete Value_Conversion;
    Value_Conversion = NULL;
  }
}

int rebalance::Zoltan::point_assign( double *position,
			  int  *proc_id ) const
{
  int status = Zoltan_LB_Point_Assign(zoltan_id, position, proc_id );

  if (status != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_LB_Point_Asssign returned status code " + status);
    return 1;
  }
  //: Sanity check for Zoltan's point assignment Processor_ID

  //ThrowAssert( ( *proc_id >= 0       ) &&
  //        ( *proc_id <  Env::parallel_size() ));
  return 0;
}


int rebalance::Zoltan::box_assign(double min[],
		       double max[],
		       std::vector<int> &procs) const
{
  /* Allocate maximum array size needed to hold processor
     numbers.
  */
  int num_procs;
  int *procbuf = new int[parallel_machine_size(comm_)];
  if (ZOLTAN_OK != Zoltan_LB_Box_Assign   (zoltan_id,
					   min[0],  min[1], min[2],
					   max[0],  max[1], max[2],
					   procbuf, &num_procs)) {
    delete [] procbuf;
    return 1;
  }
  procs.resize(num_procs);
  for (int i=0; i<num_procs; ++i) procs[i] = procbuf[i];
  delete [] procbuf;
  return 0;
}


const std::string &rebalance::Zoltan::parameter_entry_name() const
{
  //ThrowAssert (!parameter_entry_Name.empty());
  return parameter_entry_Name;
}


//: Load Balance calls or Load "Re-partitioning" calls

int rebalance::Zoltan::register_callbacks()
{
  /*
   * Register the Zoltan/SIERRA "call-back" (querry) functions.
   * Use ONLY THE BARE ESSENTIALS for decompositions:
   *
   *    Zoltan_Set_Num_Obj_Fn
   *    Zoltan_Set_First_Obj_Fn
   *    Zoltan_Set_Next_Obj_Fn
   *    Zoltan_Set_Num_Geom_Fn
   *    Zoltan_Set_Geom_Fn
   */

  /*
   * flag for what data is to be registered with the zoltan callback
   * functions in combination with the "static_cast<CLASS>(data)"
   * statement used in the zoltan interface routine,
   * Fmwk_Zoltaninterface.C
   *
   */

  if ( Zoltan_Set_Num_Obj_Fn( zoltan_id,
			      Callback_Num_Elements,
			      this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Num_Obj_Fn using Callback_Num_Elements failed to register");
  }
  if ( Zoltan_Set_Obj_List_Fn( zoltan_id, Callback_Element_List,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Next_Obj_Fn using Callback_Element_List");
  }
  if ( Zoltan_Set_First_Obj_Fn( zoltan_id, Callback_First_Element,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_First_Obj_Fn using Callback_fFirst_Element");
  }
  if ( Zoltan_Set_Next_Obj_Fn( zoltan_id, Callback_Next_Element,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Next_Obj_Fn using Callback_Next_Element");
  }
  if ( Zoltan_Set_Num_Geom_Fn( zoltan_id, Callback_Num_Dimensions,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Num_Geom_Fn using Callback_Num_Dimensions");
  }
  if ( Zoltan_Set_Geom_Fn( zoltan_id, Callback_Centroid_Coord,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Geom_Fn using Callback_Centroid_Coord");
  }
  if ( Zoltan_Set_Num_Edges_Fn( zoltan_id, Callback_Num_Edges,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Num_Edges_Fn using Callback_Num_Edges");
  }
  if ( Zoltan_Set_Edge_List_Fn( zoltan_id, Callback_Edge_List,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Edge_List_Fn using Callback_Edge_List");
  }

  return 0;

}


int  rebalance::Zoltan::evaluate( int    print_stats,
		       int*   nobj,
		       double*  obj_wgt,
		       int*   ncuts,
		       double*  cut_wgt,
		       int*   nboundary,
		       int*   nadj      )
{
  int
    ierr        = 0,
    z_nobj      = 0,
    z_ncuts     = 0,
    z_nboundary = 0,
    z_nadj      = 0;
  float
    z_obj_wgt   = 0,
    z_cut_wgt   = 0;

  if ( Zoltan_LB_Eval( zoltan_id,
		       print_stats,
		       &z_nobj,
		       &z_obj_wgt,
		       &z_ncuts,
		       &z_cut_wgt,
		       &z_nboundary,
		       &z_nadj ) != ZOLTAN_OK ) {
    ierr = 1;
  }
  *nobj      = z_nobj;
  *obj_wgt   = z_obj_wgt;
  *ncuts     = z_ncuts;
  *cut_wgt   = z_cut_wgt;
  *nboundary = z_nboundary;
  *nadj      = z_nadj;

  return ierr;

}

int rebalance::Zoltan::determine_new_partition (bool &RebalancingNeeded)
{
  //: Transfer export global ID lists and owning processors
  //: to SIERRA Framework's data structures

  /**
   * Perform the dynamic load rebalancing.  This just determines
   * the redistribution.  It does not actually redistribute
   * the objects.
   */

  int length_gid, length_lid;

  int           new_decomp;
  int           num_imported;
  ZOLTAN_ID_PTR import_gids;
  ZOLTAN_ID_PTR import_lids;
  int          *import_procs=NULL;
  int           num_exported;
  ZOLTAN_ID_PTR export_gids;
  ZOLTAN_ID_PTR export_lids;
  int          *export_procs=NULL;



      /** Zoltan_LB_Balance():
       *
       * lb               Pointer to the load-balancing structure, created by
       *                  Zoltan_Create, to be used in this invocation
       *                  of the load-balancing routine.
       * new_decomp       Set to 1 or .TRUE. if the decomposition was
       *                  changed by the load-balancing method; 0 or
       *                  .FALSE. otherwise.
       * num_gid_entries  Upon return, the number of array entries used to
       *                  describe a single global ID.  This value is the
       *                  maximum value over all processors of the parameter
       *                  NUM_GID_ENTRIES.
       * num_lid_entries  Upon return, the number of array entries used to
       *                  describe a single local ID.  This value is the
       *                  maximum value over all processors of the parameter
       *                  NUM_LID_ENTRIES.
       * num_imported     Upon return, the number of objects that are now
       *                  assigned to this processor that were assigned to
       *                  other processors in the old decomposition (i.e.,
       *                  the number of objects to be imported to this
       *                  processor). If the value returned is -1, no import
       *                  information has been returned and all import arrays
       *                  below are NULL (see the RETURN_LISTS parameter for
       *                  more information).
       * import_gids      Upon return, an array of num_import global IDs of
       *                  objects to be imported to this processor.
       *                  (size = num_import * num_gid_entries)
       * import_lids      Upon return, an array of num_import local IDs of
       *                  objects to be imported to this processor.
       *                  (size = num_import * num_lid_entries)
       * import_procs     Upon return, an array of size num_import listing
       *                  the processor IDs of the processors that owned the
       *                  imported objects in the previous decomposition
       *                  (i.e., the source processors).
       * num_exported     Upon return, the number of objects that were
       *                  assigned to this processor in the previous
       *                  decomposition that are now assigned to other
       *                  processors (i.e., the number of objects that must
       *                  be exported from this processor to other processors)
       *                  If the value returned is -1, no export information
       *                  has been returned and all export arrays below are
       *                  NULL (see the RETURN_LISTS parameter for more
       *                  information).
       * export_gids      Upon return, an array of num_export global IDs of
       *                  objects to be exported from this processor.
       *                  (size = num_export * num_gid_entries)
       * export_lids      Upon return, an array of num_export local IDs of
       *                  objects to be exported from this processor.
       *                  (size = num_export * num_lid_entries)
       * export_procs     Upon return, an array of size num_export listing
       *                  the processor IDs of processors that will own the
       *                  exported objects in the new decomposition (i.e.,
       *                  the destination processors).
       */

  int status = Zoltan_LB_Balance( zoltan_id,        &new_decomp,
				  &length_gid     , &length_lid,
				  &num_imported,    &import_gids,
				  &import_lids,     &import_procs,
				  &num_exported,    &export_gids,
				  &export_lids,     &export_procs );
  if (status != ZOLTAN_OK) {
    throw std::runtime_error("Zoltan_Balance() returned error code " + status);
  }

  //ThrowAssert (::num_gid_entries() == (unsigned) length_gid);
  //ThrowAssert (::num_lid_entries() == (unsigned) length_lid);

  //: Initialize destination processor IDs (dest_proc_ids)
  reset_dest_proc_data();

  int actual_exported = 0;
  if ( new_decomp && ( num_exported != -1 ) ) {
    const unsigned parallel_rank = parallel_machine_rank(comm_);
    /* New Decomposition was generated */
    for (int j=0; j < num_exported; ++j ) {

      //: Get exported region, local, global, and processor ids
      //const unsigned rid = export_gids[ j*::num_gid_entries() ];  // Region ID variable
      const unsigned lid = export_lids[ j*::num_lid_entries() ];  // Local  ID variable
      const unsigned pid = export_procs[ j ];                     // Exported Processor ID (i.e., MPI "rank" )

      //ThrowRequire(rid == 0);

      //: Sanity check for exported (obj) Global_ID
      //ThrowAssert( (int) export_gids[ j*::num_gid_entries() + 1 ] ==
      //        globalID(lid) );

      //: Sanity check for exported Processor_ID
      //ThrowAssert( 0 <= (int) pid && (int) pid < Env::parallel_size() );
      //: Set exported procesor_id corresponding to global_id
      if (parallel_rank != pid) {
	++actual_exported;
	set_destination_proc(lid, pid);
      }
    }
  }

  RebalancingNeeded = 0 ;
  if (new_decomp) {
    int rebalneeded=0;
    stk::all_reduce_sum(comm_, &actual_exported, &rebalneeded, 1);
    if (rebalneeded)  RebalancingNeeded = 1;
  }

  /**
   * Clean up after zoltan
   */
  if ( ZOLTAN_OK !=
       Zoltan_LB_Free_Data( &import_gids, &import_lids, &import_procs,
			    &export_gids, &export_lids, &export_procs )) {
      //throw RuntimeError() << " P" << Env::parallel_rank()
      //  		   << ": FATAL ERROR in Zoltan_LB_Free_Data error code " << ZOLTAN_FATAL << std::endl
      //  		   << StackTrace;
  }
  return EXIT_SUCCESS;
}

void rebalance::Zoltan::convert_names_and_values(const Teuchos::ParameterList &from, Teuchos::ParameterList &to)
{
  /* First time through, fill the conversion tables. */
  if (!Name_Conversion) {
    Name_Conversion = new Teuchos::ParameterList;
    //ThrowAssert (Name_Conversion);
    //fill_name_conversion (*Name_Conversion);
  }
  if (!Value_Conversion) {
    Value_Conversion = new Teuchos::ParameterList;
    //ThrowAssert (Value_Conversion);
    //fill_value_conversion(*Value_Conversion);
  }

  // NOTE: "LOAD BALANCING METHOD" default
  // is also hard coded in fill_default_value();
  std::string algorithm;
  //const ParamValue *param = from.get("LOAD BALANCING METHOD");
  const std::string keyname("LOAD BALANCING METHOD");
  if( from.isParameter(keyname) )
    algorithm = from.get<std::string>(keyname); //(const std::string &) param->get();
  else
    algorithm = "0";

  //const Teuchos::ParameterList & General = Name_Conversion->sublist("General");
  //if (!General) {
  //  throw RuntimeError()
  //    << "Unable to find required parameter set General, which contains the General parameters"
  //    << std::endl << StackTrace;
  //}

  //const Teuchos::ParameterList & Algorithm = Name_Conversion->sublist(algorithm);
  //if (!Algorithm) {
  //  throw RuntimeError() << "Unable to find required parameter set " << algorithm
  //      		 << ", which contains the algorithm parameters"
  //      		 << std::endl << StackTrace;
  //}

  // TODO - convert this loop to use Teuchos::ParameterList
  //Parameters::const_iterator
  //  from_iter  = from.begin(),
  //  from_end   = from.end();

  ///* Iterate over all of the input parameters to find proper
  //   Zoltan names and Zoltan parameter values. */
  //for (; from_iter != from_end; ++from_iter) {
  //  std::string
  //    from_name = (*from_iter).name(),
  //    to_name;

  //  /* Check to see if this is a general parameter name
  //     and if not, check if it is an algorithm specific name. */
  //  const ParamValue *name = General->get(from_name);
  //  if (name) {
  //    to_name = (const std::string &) name->get();
  //  } else {
  //    name = Algorithm->get(from_name);
  //    //if (!name) {
  //    //  throw RuntimeError() << "Unable to find the parameter " << from_name << std::endl << StackTrace;
  //    //}
  //    to_name = (const std::string &) name->get();
  //  }

  //  /* Now convert the parameter value to the correct form.
  //     Only a couple of parameters have parameter conversion.
  //     The ones converted are nested in Value_Conversion.
  //  */
  //  std::string to_value = (*from_iter).get();
  //  Teuchos::ParameterList *Value_Convert = Value_Conversion->get_nested(from_name);
  //  if (Value_Convert) to_value = (const std::string &) Value_Convert->get(to_value)->get();
  //  if (!to_name.empty()) to.set(to_name, to_value);
  //}
}

//void rebalance::Zoltan::merge_default_values(const Teuchos::ParameterList &from,
//					Teuchos::ParameterList &to)
//{
//  Teuchos::ParameterList Default_Values;
//  fill_default_value(Default_Values);
//  to.copy_from(Default_Values.nested("General") );
//  to.copy_from(from);
//}
