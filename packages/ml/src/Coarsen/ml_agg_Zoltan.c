/* ************************************************************************* */
/* Author        : Marzio Sala (SNL)                                         */
/* Date          : September 2004                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_agg_Zoltan.h"
#include "ml_agg_METIS.h"
#include "ml_viz_stats.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif

int ML_BuildReorderedDecomposition( int starting_decomposition[],
					   int reordered_decomposition[],
					   int Nrows, int Naggregates,
					   int nodes_per_aggre[],
					   USR_COMM comm );
extern ML_Operator * ML_BuildQ( int StartingNumElements,
				int ReorderedNumElements,
				int num_PDE_eqns, int,
				int * reordered_decomposition,
				double * StartingNullSpace,
				double * ReorderedNullSpace,
				int ComputeNewNullSpace,
				double * StartingBdry, double * ReorderedBdry,
				USR_COMM mpi_communicator,
				ML_Comm *ml_communicator );

extern int ML_ApplyQ(int StartingNumElements,
		     int ReorderedNumElements,
		     int NumVectors,
		     double* StartingVectors,
		     double* ReorderedVectors);

extern void ML_DestroyQ(void);
  
int ML_DecomposeGraph_with_Zoltan(ML_Operator *Amatrix,
				  int N_parts,
				  int graph_decomposition[],
				  double bdry_nodes[],
				  double [], double [], double [],
				  int);

extern int OPTIMAL_VALUE;
extern int PARMETIS_DEBUG_LEVEL;
extern int OPTIMAL_LOCAL_COARSE_SIZE;

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#if defined(HAVE_ML_ZOLTAN)
#endif

/* ======================================================================== */
/*!
 \brief Reorder the local graph for the coarser level matrix using a Zoltan.

 \param Amatrix (In) -
 pointer to the amalgamated ML_Operator to decompose

 \param N_parts (In) -
 total number of parts in which the graph has to be decomposed

 \param graph_decomposition (Out) -
 integer vector, of size Amatrix->getrow->Nrows, that in output
 will contain the global partition number of each local row.

 \param bdry_nodes (In) -

 \param N_nonzeros () -

 \param current_level (In) -
 
 
*/
/* ------------------------------------------------------------------------ */

int ML_DecomposeGraph_with_Zoltan(ML_Operator *Amatrix,
				  int N_parts,
				  int graph_decomposition[],
				  double bdry_nodes[], double old_x[], 
				  double old_y[], double old_z[],
				  int current_level)
{

  int i, j,jj,  count;
  int Nrows;
  ML_Comm * comm = Amatrix->comm;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  int *global_numbering = NULL;
  int N_procs = Amatrix->comm->ML_nprocs;
  int mypid = Amatrix->comm->ML_mypid;
  int * offsets = NULL;
  int ncon = 1;
  int ok;
  int * nodes_per_aggre = NULL, * nodes_per_aggre2 = NULL;
  double t0;
  double debug_starting_time;
  int mod;
  int Nrows_global;
  int N_dimensions;

  /* ------------------- execution begins --------------------------------- */

  t0 = GetClock();
  
  /* ********************************************************************** */
  /* some general variables                                                 */
  /* ********************************************************************** */
  
  if (old_z)
    N_dimensions = 3;
  else if (old_y)
    N_dimensions = 2;
  else
    N_dimensions = 1;

  Nrows = Amatrix->getrow->Nrows;
  
  if( N_parts <= 0 ) N_parts = 1;

#ifdef FIXME
  /* allocates memory for global_ordering */

  ML_build_global_numbering(Amatrix, comm, &global_numbering);

#ifdef ML_MPI
  MPI_Allreduce( &Nrows, &Nrows_global, 1, MPI_INT, MPI_SUM, comm );
#else
  Nrows_global = Nrows;
#endif
#endif

  /* ********************************************************************** */
  /* no need to call Zoltan if only one global aggregate is required.       */
  /* ********************************************************************** */

  if (N_parts == 1) {
    for( i=0 ; i<Nrows ; i++ ) {
      graph_decomposition[i] = 0;
    }
    return 1;
  }

  /* RAY:
   * here call to Zoltan or whatever you may want to create
   * graph_decomposition. At this point, I simply decompose it. */
  
  if (N_parts > Nrows)
    N_parts = Nrows;

  for (i = 0 ; i < Nrows ; i++) {
    graph_decomposition[i] = i % N_parts;
  }

  /* ------------------- that's all folks --------------------------------- */

  t0 = GetClock() - t0;

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("Zoltan (level %d) : time required = %e\n",
	   current_level,
	   t0 );
    
  }
  
#ifdef FIXME
  if (global_numbering)
    ML_free(global_numbering);
#endif

  /* returns the *global* number of partitions */
  return(N_parts);
  
} /* ML_DecomposeGraph_with_Zoltan */

/* ======================================================================== */
/*!
 \brief create non-smoothed aggregates using Zoltan. 

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_CoarsenZoltan(ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
			       ML_Operator **Pmatrix, ML_Comm *comm)
{
   unsigned int nbytes, length;
   int     i, j, jj, k, Nrows, exp_Nrows,  N_bdry_nodes;
   int     diff_level, Nrows_global;
   int     aggr_count, index, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, nullspace_dim;
   int     Ncoarse, count;
   int     *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     *aggr_cnt_array = NULL;
   int     level, index3, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  *new_val = NULL, epsilon;
   double  *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct ML_CSR_MSRdata *csr_data;
   double                *starting_amalg_bdry, *reordered_amalg_bdry;
   int                   Nghost;
   int                   allocated = 0, *rowi_col = NULL, rowi_N;
   double                *rowi_val = NULL;
   int Nnonzeros2 = 0;
   int optimal_value;
   ML_Operator * Pmatrix2 = NULL;
   
#if defined(OUTPUT_AGGREGATES) || defined(DDEBUG) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
   FILE *fp;
   char fname[80];
   static int level_count = 0;
   double *d2temp;
   int agg_offset, vertex_offset;
#endif
   ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
   ML_Aggregate_Options * aggr_options;
   int Nprocs;
   int * starting_offset = NULL, * reordered_offset = NULL;
   int desired_aggre_per_proc;
   int * nodes_per_aggre = NULL;
   int * starting_decomposition = NULL;
   int * reordered_decomposition = NULL;
   ML_Operator * QQ = NULL;
   ML_Operator *Pstart = NULL;
   int starting_aggr_count;
   char str[80], * str2;
   double * new_nullspace_vect = NULL;
   int * graph_decomposition = NULL;
   double debug_starting_time;
   int N_dimensions = 0;
   double* old_x = NULL;
   double* old_y = NULL;
   double* old_z = NULL;
   double* new_x = NULL;
   double* new_y = NULL;
   double* new_z = NULL;
   double* next_level_x = NULL;
   double* next_level_y = NULL;
   double* next_level_z = NULL;
   double* old_nodal_coord = NULL;
   double* new_nodal_coord = NULL;
   double* next_level_nodal_coord = NULL;
   int iaggre;
   
   /* ------------------- execution begins --------------------------------- */

   sprintf( str, "Zoltan (level %d) :", ml_ag->cur_level );
   
   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */
   
   mypid                   = comm->ML_mypid;
   Nprocs                  = comm->ML_nprocs;
   epsilon                 = ml_ag->threshold;
   num_PDE_eqns            = ml_ag->num_PDE_eqns;
   nullspace_dim           = ml_ag->nullspace_dim;
   nullspace_vect          = ml_ag->nullspace_vect;
   Nrows                   = Amatrix->outvec_leng;

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s num_PDE_eqns = %d\n",
	    str,
	    num_PDE_eqns);
   }

#ifdef ML_MPI
   MPI_Allreduce( &Nrows, &Nrows_global, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
#else
   Nrows_global = Nrows;
#endif
   
   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */
     
   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("*ML*ERR* : Nrows must be multiples of num_PDE_eqns.\n");
      exit(EXIT_FAILURE);
   }
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   Nghost = Amatrix->getrow->pre_comm->total_rcv_length;
   
   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && 8 < ML_Get_PrintLevel() )
   {
      printf("%s current eps = %e\n",
	     str,
	     epsilon);
      if( epsilon != 0.0 ) {
	fprintf( stderr,
		 "WARNING: Zoltan may not work with dropping!\n"
		 "WARNING: Now proceeding -- with fingers crossed\n" );
      }
   }
   
   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   Nrows /= num_PDE_eqns;
   Nrows_global /= num_PDE_eqns;
   Nghost /= num_PDE_eqns;
   exp_Nrows = Nrows+Nghost;

   /* record the Dirichlet boundary. unamalg_bdry[i] == 1 ==> the nodes is */
   /* a boundary node, unamalg_bdry[i] == 0 ==> the node is not            */

   nbytes = sizeof(double)*(exp_Nrows + 1);
   starting_amalg_bdry = (double *) ML_allocate(nbytes);
   if( starting_amalg_bdry == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough space to allocate %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   for (i = Nrows ; i < exp_Nrows; i++) starting_amalg_bdry[i] = 0.0;

   Nnonzeros2 = 0, N_bdry_nodes = 0;
   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);

      if (rowi_N > 1) {
        starting_amalg_bdry[i] = 0.0;
        Nnonzeros2 += rowi_N;
      } else {
	starting_amalg_bdry[i] = 1.0;
	N_bdry_nodes++;
      }
   }

   if( rowi_col != NULL ) ML_free(rowi_col );
   if( rowi_val != NULL ) ML_free(rowi_val );
   
   i = ML_Comm_GsumInt(comm, N_bdry_nodes);
   
   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s # bdry (block) nodes = %d, # (block) nodes = %d\n",
	    str,
	    i, Nrows_global);
   }
   
   /* communicate the boundary information */
   
   ML_exchange_bdry(starting_amalg_bdry,Amatrix->getrow->pre_comm,
		    Nrows,comm, ML_OVERWRITE,NULL);

   /* ********************************************************************** */
   /* allocate memory for starting_decomposition and call Zoltan to          */
   /* decompose the local                                                    */
   /* graph into the number of parts specified by the user with a call       */
   /* ML_Aggregate_Set_LocalNumber( ml, ag, level, Nparts)                   */
   /* ********************************************************************** */

   /* FIXME: only Nrows??? */
   nbytes = (Nrows+Nghost) * sizeof(int);

   if ( nbytes > 0 ) starting_decomposition = (int *)ML_allocate(nbytes);
   else              starting_decomposition = NULL;
   
   if( starting_decomposition == NULL && nbytes > 0 ) {
     
     fprintf( stderr,
	      "*ML*ERR* not enough memory to allocated %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   /* ********************************************************************** */
   /* Retrive the user's defined data. If not set, pick up default settings  */
   /* ********************************************************************** */
     
   aggr_options = (ML_Aggregate_Options *) ml_ag->aggr_options;
   
   if( aggr_options == NULL ) {

     if( mypid == 0 && 8 < ML_Get_PrintLevel() ) {
       printf("%s Using default values\n",
	      str );
     }
    
     optimal_value = ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate();
     starting_aggr_count = (int)1.0*Nrows_global/optimal_value;
     if( starting_aggr_count < 1 ) starting_aggr_count = 1;
     /* 'sto xxxcio e` un po` piu` difficile, non ne ho idea, piglio
	a caso... */
     desired_aggre_per_proc = ML_max( 128, Nrows );
       
   } else {

     /* ******************************************************************** */
     /* Retrive the user's defined choice to define the number of aggregates */
     /* For local number, it is ok.                                          */
     /* If global number of aggregates or nodes per aggregate have been      */
     /* specified, compute the local one (evenly dividing this global number)*/
     /* For those two latter cases, I suppose that the input value is the    */
     /* same on all processes (but I don't check ... )                       */
     /* ******************************************************************** */

     switch( aggr_options[ml_ag->cur_level].choice ) {

     case ML_NUM_LOCAL_AGGREGATES:
       i = aggr_options[ml_ag->cur_level].Naggregates_local;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d local aggregates (on proc 0)\n",
		 str,
		 i );
       }
#ifdef ML_MPI
       MPI_Allreduce(&i,&starting_aggr_count,1,MPI_INT,MPI_SUM,comm->USR_comm);
#else
       starting_aggr_count = i;
#endif
       break;
       
     case ML_NUM_GLOBAL_AGGREGATES:
       
       starting_aggr_count = aggr_options[ml_ag->cur_level].Naggregates_global;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d global aggregates\n",
		 str,
		 starting_aggr_count );
       }

       break;
       
     case ML_NUM_NODES_PER_AGGREGATE:

       starting_aggr_count = aggr_options[ml_ag->cur_level].Nnodes_per_aggregate;

       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d nodes per aggregate\n",
		 str,
		 starting_aggr_count );
       }
       
       if( starting_aggr_count >= Nrows_global) {

	 i = starting_aggr_count;

	 starting_aggr_count = Nrows/OPTIMAL_VALUE;
	 if( starting_aggr_count == 0) starting_aggr_count = 1;
	 
	 if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "%s WARNING : # nodes per aggregate (%d) > # nodes (%d)\n"
		    "%s WARNING : now proceeding with # aggregates = %d\n",
		    str,
		    i,
		    Nrows,
		    str,
		    starting_aggr_count);
	 }

       } else {

	 starting_aggr_count = Nrows_global/starting_aggr_count;
	 if( starting_aggr_count == 0 ) starting_aggr_count = 1;

       }
       
       break;
       
     } /* switch */

     /* reorder_flag = aggr_options[ml_ag->cur_level].reordering_flag; */

     desired_aggre_per_proc = aggr_options[ml_ag->cur_level].desired_aggre_per_proc;

     if( desired_aggre_per_proc <= 0 )
       desired_aggre_per_proc = OPTIMAL_LOCAL_COARSE_SIZE;
     
   } /* if( aggr_options == NULL )*/

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s Objective : %d aggregates over %d nodes\n",
	    str,
	    starting_aggr_count,
	    Nrows_global);
     printf("%s Objective : %d aggregates on each process\n",
	    str,
	    desired_aggre_per_proc );
   } 
   
#ifdef FIXME
   N_dimensions = ml_ag->N_dimensions;
   old_nodal_coord = ml_ag->nodal_coord[abs(diff_level)];

   /* RAY: this is the place for Zoltan... 
    * Vectors old_x, old_y and old_z contains the coordinates
    * for each (block) node. They still refer to the current
    * distribution of rows. Later, these vectors have to be redifined
    * following the shipping of information.
    */
   old_x = old_nodal_coord;
   if (N_dimensions > 1 && old_x)
     old_y = old_nodal_coord + Nrows;
   else
     old_y = 0;
   if (N_dimensions > 2 && old_x)
     old_z = old_nodal_coord + 2 * Nrows;
   else
     old_z = 0;
#endif

   /* Amatrix is the *amalgamated* matrix.
    *
    * starting_aggr_count is the number of *global* partition
    * I would like to create. The function returns the number
    * that it could actually create.
    *
    * starting_decomposition[i] will contain the *global*
    * partition ID for *local* row i. `local' means that it still
    * refers to the non-reordered (starting) row layout.
    *
    * starting_amalg_bdry is a vector containing the bc.
    *
    * old_x, old_y, and old_z are passed to the function,
    * so that they can be used with (or in substitution to) the graph.
    *
    * See also notes in the function itself. */
   starting_aggr_count =
     ML_DecomposeGraph_with_Zoltan(Amatrix, starting_aggr_count,
				   starting_decomposition,
				   starting_amalg_bdry,
				   old_x, old_y, old_z,
				   ml_ag->cur_level);
   
   /* From now on the code should not change because of Zoltan
    * until the next Zoltan comment (marked with `RAY')... */

   if( starting_aggr_count <= 0 ) {
     fprintf( stderr,
	      "*ML*ERR* Something went *very* wrong in Zoltan...\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) 
     printf("%s Using %d aggregates (globally)\n",
	    str,
	    starting_aggr_count );
   
   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
     printf("%s # aggre/ # (block) rows = %7.3f %%  (= %d/%d)\n",
	    str,
	    100.0*starting_aggr_count/Nrows_global,
	    starting_aggr_count,
	    Nrows_global);
   }
   
   /* ********************************************************************** */
   /* compute operator complexity                                            */
   /* ********************************************************************** */
   
   Nnonzeros2 = ML_Comm_GsumInt( comm, Nnonzeros2);

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
     printf("%s Total (block) nnz = %d ( = %5.2f/(block)row)\n",
	    str,
	    Nnonzeros2,1.0*Nnonzeros2/Nrows_global);
   
   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = Nnonzeros2;
      ml_ag->operator_complexity = Nnonzeros2;
   }
   else ml_ag->operator_complexity += Nnonzeros2;

   /* FIXME: erase meeeeeeeee
      fix aggr_index for num_PDE_eqns > 1 
   
   for (i = Nrows - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
         aggr_index[i*num_PDE_eqns+j] = aggr_index[i];
      }
   }
   */
   
   /* ********************************************************************** */
   /* I allocate room to copy aggr_index and pass this value to the user,    */
   /* who will be able to analyze and visualize this after the construction  */
   /* of the levels. This way, the only price we have to pay for stats and   */
   /* viz is essentially a little bit of memory.                             */
   /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
   /* I set the pointers using the ML_Aggregate_Info structure. This is      */
   /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
   /* ********************************************************************** */
   
   if( ml_ag->aggr_viz_and_stats != NULL ) {

     graph_decomposition = (int *) ML_allocate(sizeof(int)*Nrows );

     if( graph_decomposition == NULL ) {
       fprintf( stderr,
		"*ML*ERR* Not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		(int)sizeof(int)*Nrows,
		__FILE__,
	      __LINE__ );
       exit( EXIT_FAILURE );
     }

     for( i=0 ; i<Nrows ; i++ )
       graph_decomposition[i] = starting_decomposition[i];

     aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *) (ml_ag->aggr_viz_and_stats);
     aggr_viz_and_stats[ml_ag->cur_level].graph_decomposition = graph_decomposition;
     aggr_viz_and_stats[ml_ag->cur_level].Nlocal = Nrows;
     aggr_viz_and_stats[ml_ag->cur_level].Naggregates = starting_aggr_count;
     aggr_viz_and_stats[ml_ag->cur_level].local_or_global = ML_GLOBAL_INDICES;
     aggr_viz_and_stats[ml_ag->cur_level].is_filled = ML_YES;
     
   }

   /* ********************************************************************** */
   /* Compute the new distribution, so that `desired_aggre_per_proc' aggre   */
   /* are stored on each processor (up to the maximum number of aggregates). */
   /* - starting_offset : decomposition of the unknowns for the finer grid   */
   /*                     before redistribution                              */
   /* - reordered_offset : decomposition of the unknowns for the finer grid  */
   /*                      as will be after redistribution, as done by       */
   /*                      operator QQ                                       */
   /* - Nrows, new_Nrows : number of local rows for the finer grid before    */
   /*                      and after redistribution                          */
   /* ********************************************************************** */

   starting_offset  = (int *)ML_allocate( sizeof(int) * (Nprocs+1));
   reordered_offset = (int *)ML_allocate( sizeof(int) * (Nprocs+1));
   nodes_per_aggre = (int *) ML_allocate( sizeof(int) * starting_aggr_count );

   if( starting_offset == NULL || reordered_offset == NULL
       || nodes_per_aggre == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough memory\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
   ML_DecomposeGraph_BuildOffsets( Nrows, starting_offset, Nprocs,
				   Amatrix->comm->USR_comm);
   
   /* ********************************************************************** */
   /* Compute how many nodes are contained in each aggregate. This will be   */
   /* done for all the aggregates (so some communications will occur).       */
   /* ********************************************************************** */
   
   ML_CountNodesPerAggre( Nrows, starting_decomposition,
			  starting_aggr_count, nodes_per_aggre,
			  Amatrix->comm->USR_comm );

   /* ********************************************************************** */
   /* Compute how many aggregates will be stored on this process. This is    */
   /* based on the `desired_aggre_per_proc', so that the first processes will*/
   /* have about this number (and then maybe some processes will have none). */
   /* This is used to determine a reorderd offset, so that each processor    */
   /* will hold the rows of the matrix required to form the given aggregates */
   /* This new row decomposition is hold in `reordered_decomposition'        */
   /* ********************************************************************** */

   aggr_count = ML_BuildReorderedOffset(starting_offset,
					desired_aggre_per_proc,
					Nprocs, nodes_per_aggre,
					starting_aggr_count,
					reordered_offset, mypid );

   new_Nrows = reordered_offset[mypid+1] - reordered_offset[mypid];
   
   i = 0;
   if( new_Nrows > 0 ) i++;

#ifdef ML_MPI
   MPI_Reduce( &i, &j, 1, MPI_INT, MPI_SUM, 0, comm->USR_comm);
#else
   j = i;
#endif

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf( "%s Processes with at least 1 row at next level = %d\n",
	     str,
	     j );
   } 
   
   reordered_decomposition = (int *) ML_allocate( sizeof(int) * (Nrows+1) );
   if( reordered_decomposition == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough memory to allocate %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      (int)sizeof(int) * (Nrows+1),
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   ML_BuildReorderedDecomposition(starting_decomposition,
				  reordered_decomposition, Nrows,
				  starting_aggr_count, nodes_per_aggre,
				  Amatrix->comm->USR_comm );
   
   /* ********************************************************************** */
   /* Finally, built the operator QQ, moving from the reorderd decomposition */
   /* to the starting one. So, QQ will be applied to vectors (or matrices)   */
   /* whose row decomposition is the reordered one. I need QQ because we have*/
   /* P_tent = QQ * \hat{P}_tent, where P_tent is the tentative prolongator  */
   /* as used by ML (from the next level to this level, in the starting      */
   /* decomposition), and \hat{P}_tent the one based on the reordered dec.   */
   /* ********************************************************************** */

   if( nullspace_vect == NULL /*&& diff_level == 0*/ ) {
     new_nullspace_vect = nullspace_vect;
     i = 0;
   } else {
     nbytes = sizeof(double) * (new_Nrows * num_PDE_eqns * nullspace_dim );

     if( nbytes == 0 ) new_nullspace_vect = NULL;
     else {
       new_nullspace_vect = (double *) ML_allocate( nbytes );
       if( new_nullspace_vect == NULL ) {
	 fprintf( stderr,
		  "*ML*ERR* Not enough memory to allocate %d bytes\n"
		  "*ML*ERR* (file %s, line %d)\n",
		  nbytes,
		  __FILE__,
		  __LINE__ );
	 exit( EXIT_FAILURE );
       }
     }
       
     i = 1;
   }

   reordered_amalg_bdry = (double *) ML_allocate(sizeof(double)*(new_Nrows+1));
   
#ifdef ML_WITH_EPETRA
   if( nullspace_dim != num_PDE_eqns ) {
     printf("---------> Never tested with nullspace_dim != num_PDE_eqns\n"
	    "---------> Memory allocation within  ML_BuildQ to be checked...\n" );
   }
   QQ = ML_BuildQ(Nrows, new_Nrows, num_PDE_eqns, nullspace_dim,
		  reordered_decomposition,
		  nullspace_vect, new_nullspace_vect, i,
		  starting_amalg_bdry, reordered_amalg_bdry,
		  Amatrix->comm->USR_comm,
		  comm );
#else
   if( mypid == 0 ) 
     fprintf( stderr,
	      "*ML*ERR* Sorry, you cannot redistribute matrices within the Zoltan\n"
	      "*ML*ERR* aggregation without epetra. Please recompile using epetra...\n" );
   exit( EXIT_FAILURE );
#endif

   if (starting_decomposition != NULL) {
     ML_free( starting_decomposition );
     starting_decomposition = NULL;
   }
   if (reordered_decomposition != NULL) {
     ML_free( reordered_decomposition );
     reordered_decomposition = NULL;
   }
   if (starting_amalg_bdry != NULL) {
     ML_free( starting_amalg_bdry );
     starting_amalg_bdry = NULL;
   }
   if (starting_offset != NULL) {
     ML_free( starting_offset );
     starting_offset = NULL;
   }
   if (reordered_offset != NULL) {
     ML_free( reordered_offset );
     reordered_offset = NULL;
   }

   /* ********************************************************************** */
   /* Now reallocating aggr_index so that we can build the prolongator       */
   /* as if all the aggregates were local. Need some reallocations here.     */
   /* ********************************************************************** */

   nbytes = sizeof(int) * new_Nrows * num_PDE_eqns;
   
   if ( nbytes > 0 ) ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
   else              aggr_index = NULL;

   /* k is the ID of the first aggregate on this proc */
#ifdef ML_MPI
   MPI_Scan( &aggr_count, &k, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
   k -= aggr_count;
#else
   k = 0;
#endif

   if( new_Nrows != 0 ) {

     jj = 0;
     for( i=0 ; i<aggr_count ; i++ ) {
       for( j=0 ; j<nodes_per_aggre[i+k] ; j++ ) {
	 aggr_index[jj] = i;
	 jj++;
       }
     }

     if( new_Nrows != jj ) {
       fprintf( stderr,
		"*ML*ERR* something went wrong in coarsening with Zoltan:\n"
		"*ML*ERR* new_Nrows = %d, jj = %d\n"
		"*ML*ERR* (file %s, line %d)\n",
		new_Nrows, jj,
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }
   }
   
   /* ================================================================= */
   /* RAY: If previous level was not the finest, build the              *
    * coordinates for this level's nodes (that is, the aggregates of    *
    * the previous level).                                              *
    *                                                                   *
    * NOTE: I suppose at this point that all the nodes belonging to any *
    * aggregate are LOCAL. This is always verified by the redistributor.*
    *                                                                   *
    * NOTE2: variables with `old' refer to coordinates on the previos   *
    * level with the previous (old) distribution.                       *
    * Variables with `new' refer to coordinates on the previous (old)   *
    * level, with a new distribution.                                   *
    * Variables with `next_level' are self-explicatory. ;)              */
   /* ================================================================= */
   
#ifdef FIXME
   i = N_dimensions * sizeof(double) * (new_Nrows);
   if (i)
     new_nodal_coord = (double*) ML_allocate(i);
   else
     new_nodal_coord = 0;

   ml_ag->nodal_coord[abs(diff_level)] = new_nodal_coord;
   
#ifdef HAVE_ML_EPETRA
   ML_ApplyQ(Nrows, new_Nrows, N_dimensions,
	     old_nodal_coord, new_nodal_coord);
#else
   puts("I must die now. Sorry.");
   exit(0);
#endif

   i = sizeof(double) * N_dimensions * (aggr_count);
   if (i)
     next_level_nodal_coord = (double*)ML_allocate(i);
   else
     next_level_nodal_coord = 0;
   /* set up the nodal coordinates for the next level. 
    * NOTE: this would create a small over(ab)use of memory
    * as the last level will still create the nodal coordinates
    * for the next level (that will never exist). However,
    * the coarse grid should be small enough at that point. */
   ml_ag->nodal_coord[abs(diff_level) + 1] = next_level_nodal_coord;

   new_x = new_nodal_coord;
   if (N_dimensions > 1 && new_nodal_coord)
     new_y = new_nodal_coord + new_Nrows;
   else
     new_y = 0;
   if (N_dimensions > 2 && new_nodal_coord)
     new_z = new_nodal_coord + 2 * new_Nrows;
   else
     new_z = 0;

   next_level_x = next_level_nodal_coord;
   if (N_dimensions > 1 && next_level_nodal_coord)
     next_level_y = next_level_nodal_coord + aggr_count;
   if (N_dimensions > 2 && next_level_nodal_coord)
     next_level_z = next_level_nodal_coord + 2 * aggr_count;

   /* zero-out new coordinates */
   for (i = 0 ; i < aggr_count ; ++i) {
     next_level_x[i] = 0.0;
     if (next_level_y) 
       next_level_y[i] = 0.0;
     if (next_level_z) 
       next_level_z[i] = 0.0;
   }

   /* sum up all the contributions for each node in the aggregate */
   for (i = 0 ; i < new_Nrows ; ++i) {
     iaggre = aggr_index[i];
     if (iaggre < 0 || iaggre > aggr_count)
       continue;
     next_level_x[iaggre] += new_x[i];
     if (next_level_y)  
       next_level_y[iaggre] += new_y[i];
     if (next_level_z)  
       next_level_z[iaggre] += new_z[i];
   }

   for (i = 0 ; i < aggr_count ; ++i) {
     next_level_x[i] = next_level_x[i] / nodes_per_aggre[i];
     if (next_level_y)  
       next_level_y[i] = next_level_y[i] / nodes_per_aggre[i];
     if (next_level_z) 
       next_level_z[i] = next_level_z[i] / nodes_per_aggre[i];
   }

   /* =================================== */
   /* RAY: END OF COORDINATES COMPUTATION */
   /* =================================== */

   /* free memory */
   if (old_nodal_coord) 
     ML_free(old_nodal_coord);
#endif

   if( nodes_per_aggre != NULL ) {
     ML_free( nodes_per_aggre );
     nodes_per_aggre = NULL;
   }

   /* *************************** */
   /* TO ADD: OPERATOR COMPLEXITY */
   /* *************************** */

   for (i = new_Nrows - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
         aggr_index[i*num_PDE_eqns+j] = aggr_index[i];
      }
   }

   if ( mypid == 0 && 8 < ML_Get_PrintLevel())
   {
      printf("Calling ML_Operator_UnAmalgamateAndDropWeak\n");
      fflush(stdout);
   }

   ML_Operator_UnAmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   
   new_Nrows  *= num_PDE_eqns;
   Nrows_global  *= num_PDE_eqns;

   /* count the size of each aggregate. Now all aggregates are local */
   
   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*aggr_count);
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < new_Nrows; i++) 
      if (aggr_index[i] >= 0) 
         aggr_cnt_array[aggr_index[i]]++;

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */
   
   Ncoarse = aggr_count;
   
   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = new_Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < new_Nrows; i+=num_PDE_eqns ) 
   {
      if ( aggr_index[i] >= 0 )
      {
         for ( j = 0; j < num_PDE_eqns; j++ ) 
            ml_ag->aggr_info[level][i+j] = aggr_index[i];
         if (aggr_index[i] >= count) count = aggr_index[i] + 1;
      }
      /*else
       *{
       *   printf("%d : CoarsenMIS error : aggr_index[%d] < 0\n",
       *          mypid,i);
       *   exit(1);
       *}*/
   }
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */

   /* ============================================================= */
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= Ncoarse ) 
      {
         printf("*ML*WRN* index out of bound %d = %d (%d)\n"
		"*ML*WRN* (file %s, line %d)\n",
		i, aggr_index[i], 
                Ncoarse,
		__FILE__,
		__LINE__ );
      }
   }
   nbytes = ( new_Nrows + 1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
   nbytes = new_Nrows * nullspace_dim * sizeof(int); 
   ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
   nbytes = new_Nrows * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */
   
   nbytes = Ncoarse * nullspace_dim * nullspace_dim * sizeof(double);
   if( nbytes != 0 ) {
     ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
     for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
       new_null[i] = 0.0;
   } else {
     new_null = NULL;
   }
   
   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= new_Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* trying this when a Dirichlet row is taken out */
   j = 0;
   new_ia[0] = 0;
   for (i = 0; i < new_Nrows; i++) {
      if (aggr_index[i] != -1) j += nullspace_dim;
      new_ia[i+1] = j;
   }

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLs");
   for (i = 0; i < aggr_count; i++) 
   {
      rows_in_aggs[i] = (int *) ML_allocate(aggr_cnt_array[i]*sizeof(int));
      aggr_cnt_array[i] = 0;
      if (rows_in_aggs[i] == NULL) 
      {
         printf("*ML*ERR* couldn't allocate memory in CoarsenZoltan\n");
         exit(1);
      }
   }
   for (i = 0; i < new_Nrows; i+=num_PDE_eqns) 
   {
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
      {
         for (j = 0; j < num_PDE_eqns; j++)
         {
            index = aggr_cnt_array[aggr_index[i]]++; 
            rows_in_aggs[aggr_index[i]][index] = i + j;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   else             qr_tmp = NULL;
   nbytes = nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");
   else             tmp_vect = NULL;
   
   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&work, nbytes, "AGw");
   else             work = NULL;
  
   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */
     
   for (i = 0; i < aggr_count; i++) 
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */

      length = aggr_cnt_array[i];

      if (new_nullspace_vect == NULL) 
      {
         for (j = 0; j < (int) length; j++)
         {
            index = rows_in_aggs[i][j];
            for (k = 0; k < nullspace_dim; k++)
            {
              if ( reordered_amalg_bdry[index/num_PDE_eqns] == 1.0) qr_tmp[k*length+j] = 0.;
               else
               {
                  if (index % num_PDE_eqns == k) qr_tmp[k*length+j] = 1.0;
                  else                           qr_tmp[k*length+j] = 0.0;
               }
            }
         }
      }
      else 
      {
	
	for (k = 0; k < nullspace_dim; k++) 
         {
	   
            for (j = 0; j < (int) length; j++)
            {
               index = rows_in_aggs[i][j];
	       
               if ( reordered_amalg_bdry[index/num_PDE_eqns] == 1.0) qr_tmp[k*length+j] = 0.;
               else {
                  if (index < new_Nrows) {
		    qr_tmp[k*length+j] = new_nullspace_vect[k*new_Nrows+index];
                  }
                  else {
		    fprintf( stderr,
			     "*ML*ERR* error in QR factorization within Zoltan aggregation\n"
			     "*ML*ERR* (file %s, line %d)\n",
			     __FILE__,
			     __LINE__ );
		    exit( EXIT_FAILURE );
                  }
               }
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */
      if (aggr_cnt_array[i] >= nullspace_dim) {

	DGEQRF_F77(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
			  &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0)
	  pr_error("ErrOr in CoarsenZoltan : "
		   "dgeqrf returned a non-zero %d %d\n",
		   aggr_cnt_array[i],i);

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGx");
	  }
	else lwork=(int) work[0];
		 
	/* ---------------------------------------------------------- */
	/* the upper triangle of qr_tmp is now R, so copy that into   */
	/* the new nullspace                                          */
	/* ---------------------------------------------------------- */

	for (j = 0; j < nullspace_dim; j++)
	  for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
		 
	/* ---------------------------------------------------------- */
	/* to get this block of P, need to run qr_tmp through another */
	/* LAPACK function:                                           */
	/* ---------------------------------------------------------- */

	if ( aggr_cnt_array[i] < nullspace_dim ){
	  printf("*ML*ERR* in dorgqr on %d row (dims are %d, %d)\n",
		 i,aggr_cnt_array[i],
                 nullspace_dim);
	  printf("*ML*ERR* performing QR on a MxN matrix where M<N.\n");
	}
	DORGQR_F77(&(aggr_cnt_array[i]), &nullspace_dim,
			  &nullspace_dim, qr_tmp, &(aggr_cnt_array[i]),
			  tmp_vect, work, &lwork, &info);
	if (info != 0) {
	  printf("*ML*ERR* in dorgqr on %d row (dims are %d, %d)\n",
		 i,aggr_cnt_array[i],
                 nullspace_dim);
	  pr_error("Error in CoarsenZoltan: dorgqr returned a non-zero\n");
	}

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGy");
	  }
	else lwork=(int) work[0];

	/* ---------------------------------------------------------- */
	/* now copy Q over into the appropriate part of P:            */
	/* The rows of P get calculated out of order, so I assume the */
	/* Q is totally dense and use what I know of how big each Q   */
	/* will be to determine where in ia, ja, etc each nonzero in  */
	/* Q belongs.  If I did not assume this, I would have to keep */
	/* all of P in memory in order to determine where each entry  */
	/* should go                                                  */
	/* ---------------------------------------------------------- */

	for (j = 0; j < aggr_cnt_array[i]; j++)
	  {
	    index = rows_in_aggs[i][j];
	    if ( index < new_Nrows )
	      {
		index3 = new_ia[index];
		for (k = 0; k < nullspace_dim; k++) 
		  {
		    new_ja [index3+k] = i * nullspace_dim + k;
		    new_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
		  }
	      }
	    else 
	      {
		fprintf( stderr,
			 "*ML*ERR* error in QR factorization within Zoltan\n" );
		exit( EXIT_FAILURE );
	      }
	  }
      }
      else {
	/* We have a small aggregate such that the QR factorization can not */
	/* be performed. Instead let us copy the null space from the fine   */
        /* into the coarse grid nullspace and put the identity for the      */
	/* prolongator????                                                  */
	for (j = 0; j < nullspace_dim; j++)
	  for (k = 0; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
	for (j = 0; j < aggr_cnt_array[i]; j++) {
	  index = rows_in_aggs[i][j];
	  index3 = new_ia[index];
	  for (k = 0; k < nullspace_dim; k++) {
	    new_ja [index3+k] = i * nullspace_dim + k;
	    if (k == j) new_val[index3+k] = 1.;
	    else new_val[index3+k] = 0.;
	  }
	}
      }

   } /* for( i over aggregates ) */
   
   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   if( new_null != NULL ) {
     ML_memory_free( (void **) &new_null);
     new_null = NULL;
   }
   
   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   Pstart = ML_Operator_Create( Amatrix->comm );
			       
   ML_Operator_Set_ApplyFuncData( Pstart, nullspace_dim*Ncoarse, new_Nrows, 
                                  csr_data, new_Nrows, NULL, 0);
   Pstart->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;

   Pstart->getrow->pre_comm = ML_CommInfoOP_Create();
   
   ML_Operator_Set_Getrow((Pstart), new_Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc(Pstart, CSR_matvec);
   Pstart->max_nz_per_row = 1;

   Pmatrix2 = ML_Operator_Create( Amatrix->comm );
   
   ML_2matmult(QQ, Pstart, Pmatrix2, ML_CSR_MATRIX );
   
   ML_Operator_Set_1Levels(Pmatrix2, (*Pmatrix)->from, (*Pmatrix)->to);
   ML_Operator_Set_BdryPts(Pmatrix2, (*Pmatrix)->bc);
   str2 = (char *)ML_allocate(80*sizeof(char));
   sprintf(str2,"%s",(*Pmatrix)->label);
   ML_Operator_Set_Label( Pmatrix2,str2);
   
   ML_free(str2);

   ML_Operator_Clean( *Pmatrix );

   memcpy((void *) *Pmatrix, (void *)Pmatrix2, sizeof(ML_Operator));
   /* FIXME : am I ok  ????? */
   ML_free(Pmatrix2);
      
   /* ********************************************************************** */
   /* I have to destroy the tentative local matrix, and the redistribution   */
   /* matrix QQ. This is actually an ML_Operator on the top of an Epetra     */
   /* object. So, I call ML_DestroyQ, which is a CPP function, to delete the */
   /* memory interally used by Epetra.                                       */
   /* ********************************************************************** */

   ML_Operator_Destroy( &Pstart ); 
   ML_Operator_Destroy( &QQ );

#ifdef ML_WITH_EPETRA
   ML_DestroyQ( );
#endif
   
   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &aggr_index);
   if( aggr_cnt_array != NULL ) ML_free(aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   if( qr_tmp != NULL ) {
     ML_memory_free((void**)&qr_tmp);
     qr_tmp = NULL;
   }
   if( tmp_vect != NULL ) {
     ML_memory_free((void**)&tmp_vect);
     tmp_vect = NULL;
   }
   if( work != NULL ) {
     ML_memory_free((void**)&work);
     work = NULL;
   }

   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }

   if( reordered_amalg_bdry != NULL ) {
     ML_free( reordered_amalg_bdry );
     reordered_amalg_bdry = NULL;
   }
   if( new_nullspace_vect != NULL ) {
     ML_free( new_nullspace_vect );
     new_nullspace_vect = NULL;
   }

   /* ------------------- that's all folks --------------------------------- */

   return Ncoarse*nullspace_dim;

} /* ML_Aggregate_CoarsenZoltan */
