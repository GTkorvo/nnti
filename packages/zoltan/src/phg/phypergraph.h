/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef ZOLTAN_HYPERGRAPH_H
#define ZOLTAN_HYPERGRAPH_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include "zz_const.h"
#include "params_const.h"
#include "phg_hypergraph.h"
#include "phg_util.h"

#define MIN(A,B)  (((A) < (B)) ? (A) : (B))
#define MAX(A,B)  (((A) > (B)) ? (A) : (B))

#define RANDOM_SEED   123456789   /* (time ((time_t*)NULL)) */
#define EPS           1e-6        /* small value, like machine epsilon */



typedef int *Matching;  /* length |V|, matching information of vertices */
/* If a vertex i is not being contracted (matched) with other vertices,  
 * Matching[i] == i.  If vertices i, j, and k are being contracted together to
 * form one new vertex, Matching[i] == j; Matching[j] == k; and Matching[k] == i; 
 * The cycle (of any length) describes the contraction. */
 
typedef int *LevelMap;  /* length |V|, mapping of fine vertices onto coarse vertices */
typedef int *Partition; /* length |V|, partition ID for each vertex */

typedef struct {
   PHGraph   *hg;
   Partition part;
   Matching  match;
   LevelMap  lmap;
} PHGraphLevel;

/* Function types for options to hypergraph partitioning */
struct PHGPartParamsStruct;  /* Forward declaration */
typedef int ZOLTAN_PHG_MATCHING_FN   (ZZ*, PHGraph*, Matching, int*);
typedef int ZOLTAN_PHG_SERIALPARTITION_FN(ZZ*, PHGraph*, int, Partition,
                                     struct PHGPartParamsStruct*);
typedef int ZOLTAN_PHG_REFINEMENT_FN  (ZZ*, PHGraph*, int, Partition,
                                     struct PHGPartParamsStruct*, float);

/* Function types for edge-weight scaling functions. Placeholder for now; */
/* if these functions end up having the same argument list for each type, */
/* do not need separate types here or separate pointers in HGPartParams.  KDD */
typedef int ZOLTAN_PHG_MATCHING_EWS_FN (ZZ*, PGraph*);


/* Parameters to the hypergraph functions */
struct PHGPartParamsStruct {
  float bal_tol;                       /* Balance tolerance in % of average */
  int kway;                             /* 1 -> direct kway, 0->recursive bisection */
  int redl;                             /* Reduction limit (constant). */

  char redm_str[MAX_PARAM_STRING_LEN];  /* Reduction method string. */
  ZOLTAN_PHG_MATCHING_FN *matching;      /* Pointers to Matching, Packing and */

  char redmo_str[MAX_PARAM_STRING_LEN];  /* Matching optimization string*/
  ZOLTAN_PHG_MATCHING_FN *matching_opt;  /* Pointers to Matching, Packing and */
  int ews;                               /* type of hyperedge weight scaling */
                        
  char global_str[MAX_PARAM_STRING_LEN]; /* Global partitioning string and */
  ZOLTAN_PHG_SERIALPARTITION_FN *global_part;/* pointer to Global partitioning fn */

  char local_str[MAX_PARAM_STRING_LEN];   /* Local refinement string and */
  ZOLTAN_PHG_REFINEMENT_FN *local_ref;    /* pointer to Local refinement fn */

  int check_graph;      /* Flag indicating whether the input hypergraph should 
                         * be checked for errors. */
  int output_level;     /* Flag indicating amount of output from HG algorithms.
                         * See levels PHG_DEBUG_* below.  */
};
typedef struct PHGPartParamsStruct PHGPartParams;

/* Hypergraph output levels: */
#define PHG_DEBUG_NONE 0
#define PHG_DEBUG_LIST 1
#define PHG_DEBUG_ALL  2
#define PHG_DEBUG_PRINT 3
#define PHG_DEBUG_PLOT 4


int Zoltan_PHG_Set_Part_Options  (ZZ*, PHGPartParams*);
int Zoltan_PHG_HPart_Lib    (ZZ*, PHGraph*, int, Partition, PHGPartParams*, int);
int Zoltan_PHG_HPart_Info   (ZZ*, PHGraph*, int, Partition, PHGPartParams*);
double Zoltan_PHG_hcut_size_total (PHGraph*, Partition);
double Zoltan_PHG_hcut_size_links (ZZ*, PHGraph*, Partition);
double Zoltan_PHG_HPart_balance   (ZZ*, PHGraph*, int, Partition);

/* Scale Edge Weight */
int Zoltan_PHG_Scale_HGraph_Weight (ZZ*, PHGraph*, float*, int);
int Zoltan_PHG_Scale_Graph_Weight  (ZZ*, PGraph*,  float*, int);

/* Matching functions */
int Zoltan_PHG_Matching (ZZ*, PHGraph*, Matching, PHGPartParams*, int*);
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams*);

/* Coarsening */
int Zoltan_PHG_Coarsening   (ZZ*, PHGraph*, Packing, PHGraph*, int*);

/* Serial Partitioning functions */
int Zoltan_PHG_SerialRefinement (ZZ*, PHGraph*, int, Partition, PHGPartParams*);
ZOLTAN_PHG_GLOBAL_PART_FN *Zoltan_PHG_Set_Global_Part_Fn(char*);

/* Refinement functions */ /* KDD Placeholder for later. */
int Zoltan_PHG_Refinement (ZZ*, PHGraph*, int, Partition, PHGPartParams*);
ZOLTAN_PHG_REFINEMENT_FN *Zoltan_PHG_Set_Local_Ref_Fn(char*);

/* Sorting */
void Zoltan_quicksort_pointer_dec_float_int (int*, float*, int*, int, int);
void Zoltan_quicksort_pointer_dec_float     (int*, float*, int,  int);
void Zoltan_quicksort_pointer_inc_int_int   (int*, int*,   int*, int, int);
void Zoltan_quicksort_list_inc_int          (int*, int,    int);
void Zoltan_quicksort_pointer_inc_int_mult  (int*, int,    int,  int*, int*);

/* Heap datastructure */
typedef struct {
   int    space;
   int    n;
   int   *ele;
   int   *pos;
   float *value;
} HEAP;

#define Zoltan_PHG_heap_empty(H)         (((H)->n)==0)
#define Zoltan_PHG_heap_not_empty(H)     (((H)->n)!=0)
#define Zoltan_PHG_heap_max_value(H)     ((H)->value[(H)->ele[0]])
#define Zoltan_PHG_heap_peek_max(H)      ((H)->ele[0])
#define Zoltan_PHG_heap_count(H)         ((H)->n)

int  Zoltan_PHG_heap_init         (ZZ*, HEAP*, int);
void Zoltan_PHG_heap_clear        (HEAP*);
void Zoltan_PHG_heap_free         (HEAP*);
int  Zoltan_PHG_heap_check        (HEAP*);
int  Zoltan_PHG_heap_input        (HEAP*, int, float);
int  Zoltan_PHG_heap_make         (HEAP*);
int  Zoltan_PHG_heap_change_value (HEAP*, int, float);
int  Zoltan_PHG_heap_extract_max  (HEAP*);
int  Zoltan_PHG_heap_extract      (HEAP*, int);

int  Zoltan_PHG_move_vertex (PHGraph*, int, int, int, int*, int**, double*, HEAP*);
void Zoltan_PHG_Plot(int, int, int, int*, int*, int*, char*);
void Zoltan_PHG_Free_Structure (ZZ *zz);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
