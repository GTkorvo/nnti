/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef _DR_CONST_H
#define _DR_CONST_H


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include "zoltan.h"

/*****************************************************************************
 *  Definitions for the LB library driver program.
 *****************************************************************************/

#define DRIVER_NAME "zdrive"
#define VER_STR "1.0"


/* If it doesn't get defined in stdio.h then use this as a default */
#ifndef FILENAME_MAX
#define FILENAME_MAX    1024
#endif

#define MAX_NP_ELEM	27 /* max nodes per element */
#define MAX_DIM		 3 /* max number of coordinate dimensions */
#define MAX_CPU_WGTS	10 /* max number of cpu weights */

enum LISTS {  /* NULL lists to pass to Zoltan_Migrate */
  NONE = 0,
  IMPORT_LISTS,
  EXPORT_LISTS
};

enum DATA_TYPE {
  MESH = 0,
  GRAPH,
  HYPERGRAPH
};

/*
 * Structure used to describe an element. Each processor will
 * allocate an array of these structures.
 */
struct Element_Description
{
  int      border;	/* set to 1 if this element is a border element */
  int      globalID;	/* Global ID of this element, the local ID is the
			   position in the array of elements              */
  int      elem_blk;    /* element block number which this element is in */
  int      my_part;     /* Partition to which the element is assigned; 
                           default is processor number. */
  int      perm_value;  /* Value for this element in permutation vector 
                           generated by Zoltan_Order.   
                           Default is -1 (no ordering done). */
  int      invperm_value;/* Value for this element in inverse permutation vector
                           generated by Zoltan_Order.   
                           Default is -1 (no ordering done). */
  float    cpu_wgt[MAX_CPU_WGTS]; /* the computational weight(s) associated with the elem */
  float    mem_wgt;	/* the memory weight associated with the elem */
  double   avg_coord[3];/* average coordinates (centroid) for the element */
  float  **coord;	/* array for the coordinates of the element.
                             for Nemesis meshes, nodal coordinates are stored;
                             for Chaco graphs with geometry, one set of coords
                                 is stored. */
  int     *connect;	/* list of nodes that make up this element, the node
			   numbers in this list are global and not local    */
  int     *adj;		/* list of adjacent elements .
                           For Nemesis input, the list is ordered by
                           side number, to encode side-number info needed to
                           rebuild communication maps.  Value -1 represents 
                           sides with no neighboring element (e.g., along mesh
                           boundaries).  Chaco doesn't have "sides," so the 
                           ordering is irrelevent for Chaco input. */
  int     *adj_proc;	/* list of processors for adjacent elements */
  float   *edge_wgt;	/* edge weights for adjacent elements */
  int      nadj;	/* number of entries in adj */
  int      adj_len;	/* allocated length of adj/adj_proc/edge_wgt arrays */
};
typedef struct Element_Description ELEM_INFO;
typedef struct Element_Description *ELEM_INFO_PTR;

/*
 * structure for general mesh information
 */
/* Structure used to store information about the mesh */
struct Mesh_Description
{
  enum DATA_TYPE data_type;     /* Type of data stored in this data structure,
                                   based on input file type.
                                   Valid types are MESH, GRAPH, or HYPERGRAPH.*/
  int     vwgt_dim;             /* number of weights per element.            */
  int     ewgt_dim;             /* number of weights per graph edge.         */
  int     num_nodes;		/* number of nodes on this processor         */
  int     num_elems;		/* number of elements on this processor      */
  int     num_dims;		/* number of dimensions for the mesh         */
  int     num_el_blks;		/* number of element blocks in the mesh      */
  int     num_node_sets;	/* number of node sets in the mesh           */
  int     num_side_sets;	/* number of side sets in the mesh           */
  char  **eb_names;		/* element block element names               */
  int    *eb_etypes;            /* element block element types               */
  int    *eb_ids;		/* element block ids                         */
  int    *eb_cnts;		/* number of elements in each element block  */
  int    *eb_nnodes;		/* number of nodes per element in each
				   element block                           
                                      for Nemesis meshes, this value depends
                                          on element type;
                                      for Chaco graphs, only one "node" per
                                          element.                           */
  int    *eb_nattrs;		/* number of attributes per element in each
				   element block                             */
  int     necmap;               /* number of elemental communication maps.   */
  int    *ecmap_id;             /* IDs of each elemental communication map.  */
  int    *ecmap_cnt;            /* number of elements in each elemental
                                   communication map.                        */
  int    *ecmap_elemids;        /* element ids of elements for all elemental
                                   communication maps. (local numbering)     */
  int    *ecmap_sideids;        /* side ids of elements for all elemental 
                                   communication maps.                       */
  int    *ecmap_neighids;       /* elements ids of neighboring elements 
                                   for all elemental communication maps. 
                                   (global numbering)                        */
  int     elem_array_len;	/* length that the ELEM_INFO array is
				   allocated for. Need to know this when array
				   is not completely filled during migration */
  ELEM_INFO_PTR elements;       /* array of elements that are in the mesh.   */
  int     gnhedges;             /* for hypergraphs, the number of global
                                   hyperedges.*/
  int     nhedges;              /* for hypergraphs, the number of local 
                                   hyperedges.*/
  int     hewgt_dim;            /* for hypergraphs, the number of weights per
                                   hyperedge.                                */
  int    *hgid;                /* Global numbering for hyperedges; numbering
                                   is derived implicitly from order hyperedges
                                   are read from file. Numbering is 0-based. */
  int    *hindex;               /* for hypergraphs, an entry for each 
                                   hyperedge, giving the starting index into
                                   hvertex for that hyperedge.               */ 
  int    *hvertex;              /* for hypergraphs, an array listing the 
                                   vertices making up each hyperedge.        */
  float  *hewgts;               /* for hypergraphs, an array of hyperedge
                                   weights; size = hewgt_dim * nhedges;      */

};
typedef struct Mesh_Description  MESH_INFO;
typedef struct Mesh_Description *MESH_INFO_PTR;

/* Structure for the problem description. */
struct Parameter_Description {
  char Name[128];     /* parameter name */
  char Val[128];      /* parameter value */
  int  Index;         /* index for vector params */
};
typedef struct Parameter_Description Parameter_Pair;

struct Problem_Description
{
  char   method[32];                 /* this is the method string that will
                                        be passed unchanged to Zoltan        */
  int num_params;                    /* number of parameters read.           */
  Parameter_Pair *params;            /* parameter array to be passed to
                                        Zoltan.  Parameters are specified as
                                        pairs of strings:
                                          param_str = value_str              */
};
typedef struct Problem_Description  PROB_INFO;
typedef struct Problem_Description *PROB_INFO_PTR;

/* Structure for driver flags for various test options. */
struct Test_Flags {
  int DDirectory;           /* Exercises data directories */
  int Local_Partitions;     /* Sets NUM_LOCAL_PARTITIONS parameter in various
                               ways. */
  int Drops;                /* Exercises point- and box-assign. */
  int Multi_Callbacks;      /* Exercises list-based callback functions. */
  int Gen_Files;            /* Exercise output file generation. */
  int Null_Lists;           /* Exercises null import or export lists to
                               Zoltan_Migrate. */
};

/* Structure for output flags for various types of output. */
struct Output_Flags {
  int Gnuplot;
  int Nemesis;
  int Plot_Partitions;
  int Mesh_Info_File;
};

/* Global variables for driver */
extern int Debug_Driver;
extern int Debug_Chaco_Input;
extern int Number_Iterations;
extern int Driver_Action;
extern int Chaco_In_Assign_Inv;
extern struct Test_Flags Test;
extern struct Output_Flags Output;

extern double Total_Partition_Time;

#define DEBUG_TRACE_START(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER ENTERING %s\n", (proc), yo);
#define DEBUG_TRACE_END(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER LEAVING %s\n", (proc), yo);
#define DEBUG_TRACE_DETAIL(proc,yo,str) \
  if (Debug_Driver > 2) \
    printf("%d DRIVER %s: %s\n", proc,yo, str);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif /* _DR_CONST_H */
