/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* -- SuperLU routine (version 1.1) --                                  */
/* Univ. of California Berkeley, Xerox Palo Alto Research Center,       */
/* and Lawrence Berkeley National Lab.                                  */
/* November 15, 1997                                                    */
/* ******************************************************************** */
/* ML interface to SuperLU :                                            */
/*   Author : Charles Tong (LLNL) and Ray Tuminaro (SNL)                */
/*   Date   : April, 1998                                               */
/* ******************************************************************** */

#ifdef SUPERLU
#include "dsp_defs.h"
#include "util.h"
#elif DSUPERLU
#include "mpi.h"
#include <malloc.h>
#include "superlu_ddefs.h"
#endif

#include "ml_comm.h"
#include "ml_struct.h"
#include "ml_solver.h"
#define REPLACE 1
#define NEQU    0

extern void ML_SuperLU_Set_Tile( int nprocs, int* tsz, int* stile, 
                                 int* mtile, int* ltile);
/*
extern int heap_info(int*, int*, int*, int*);
extern int get_heap_info(int*, int*, int*, int*);
*/

/* ******************************************************************** */
/* This subroutine calls the SuperLU subroutine to perform LU           */
/* factorization of a given matrix                                      */
/* ******************************************************************** */

int ML_SuperLU_Solve(void *vsolver,int ilen,double *x,int olen,double *rhs)
{
#ifdef SUPERLU
   int            i, n, info, flag, *perm_r, *perm_c, permc_spec;
   int            N_local, offset, *etree, panel_size, lwork;
   double         *local_x, *local_rhs, *R, *C, *ferr, *berr;
   double         rpg, rcond;
   char           fact[1], equed[1], trans[1], refact[1];
   void           *work=NULL;
   ML_Comm        *comm;
   ML_Solver      *solver;
   SuperMatrix    *A, *L, *U, B, X;
   factor_param_t iparam;
   mem_usage_t    mem_usage;

#ifdef ML_TIMING_DETAILED
   double         t0;
#endif

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen ) 
   {
      printf("ML_SuperLU_Solve error : lengths not matched.\n");
      exit(1);
   }
   solver   = (ML_Solver *) vsolver;

   comm     = (ML_Comm *) solver->void_params1;
   N_local  = (int) solver->dble_params1[0];
   offset   = (int) solver->dble_params1[1];
   n        = (int) solver->dble_params1[2];
   flag     = solver->reuse_flag;

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* ------------------------------------------------------------- */

   if ( flag == -999 ) 
   {
      A = (SuperMatrix *) solver->Mat1;
      if ( A != NULL )
      {
         ML_memory_free((void*)&(((NRformat *) A->Store)->colind));
         ML_memory_free((void*)&(((NRformat *) A->Store)->rowptr));
         ML_memory_free((void*)&(((NRformat *) A->Store)->nzval));
         SUPERLU_FREE( ((NRformat *) A->Store)->nzval);
         SUPERLU_FREE( A->Store );
         ML_memory_free((void**) &A);
      }
      solver->Mat1 = NULL;
      L = (SuperMatrix *) solver->Mat2;
      if ( L != NULL ) 
      {
         Destroy_SuperNode_Matrix(L);
         ML_memory_free((void**) &L);
      }
      solver->Mat2 = NULL;
      U = (SuperMatrix *) solver->Mat3;
      if ( U != NULL )
      {
         SUPERLU_FREE( ((NRformat *) U->Store)->colind);
         SUPERLU_FREE( ((NRformat *) U->Store)->rowptr);
         SUPERLU_FREE( ((NRformat *) U->Store)->nzval);
         SUPERLU_FREE( U->Store );
         ML_memory_free((void**) &U);
      }
      solver->Mat3 = NULL;
      perm_r = (int *) solver->int_params1;
      if (perm_r != NULL) ML_memory_free((void**) &(solver->int_params1));
      solver->int_params1 = NULL;
      perm_c = (int *) solver->int_params2;
      if (perm_c != NULL) ML_memory_free((void**) &(solver->int_params2));
      solver->int_params2 = NULL;
#ifdef ML_TIMING_DETAILED
      if ( comm->ML_mypid == 0 )
         printf("Total SuperLU solve time = %e\n", solver->dble_data);
#endif
      return 0;
   } 
   else if ( flag == 0 ) 
   {
      A = (SuperMatrix *) solver->Mat1;
      ML_memory_alloc((void**) &perm_c, n * sizeof(int), "LU2" );
      ML_memory_alloc((void**) &perm_r, 2 * n * sizeof(int), "LU3" );
      solver->int_params1 = perm_r;
      solver->int_params2 = perm_c;
      permc_spec = 2;
      get_perm_c(permc_spec, A, perm_c);
      ML_memory_alloc((void**) &L, sizeof(SuperMatrix), "LU4" );
      ML_memory_alloc((void**) &U, sizeof(SuperMatrix), "LU5" );
      solver->Mat2 = (void *) L;
      solver->Mat3 = (void *) U;
   } 
   else 
   {
      A = (SuperMatrix *) solver->Mat1;
      L = (SuperMatrix *) solver->Mat2;
      U = (SuperMatrix *) solver->Mat3;
      perm_r = (int *) solver->int_params1;
      perm_c = (int *) solver->int_params2;
   }

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   ML_memory_alloc((void**) &local_x,   n*sizeof(double),"LU2" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) comm, local_rhs, &i, n);

   /* ------------------------------------------------------------- */
   /* create right hand side vector conforming to SuperLU format    */
   /* ------------------------------------------------------------- */

   dCreate_Dense_Matrix(&B, n, 1, local_rhs, n, DN, _D, GE);
   dCreate_Dense_Matrix(&X, n, 1, local_x,   n, DN, _D, GE);

   /* ------------------------------------------------------------- */
   /* perform LU decomposition and then solve                       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &etree,   n*sizeof(int),"LU3" );
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   if ( flag == 0 )   *fact = 'N';
   else               *fact = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R    = (double *) SUPERLU_MALLOC(n * sizeof(double));
   C    = (double *) SUPERLU_MALLOC(n * sizeof(double));
   ferr = (double *) SUPERLU_MALLOC(sizeof(double));
   berr = (double *) SUPERLU_MALLOC(sizeof(double));

#ifdef ML_TIMING_DETAILED
   t0 = GetClock();
#endif

   dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
          equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
          ferr, berr, &mem_usage, &info);

#ifdef ML_TIMING_DETAILED
   t0 = GetClock() - t0;
   solver->dble_data += t0;
#endif

   if ( info != 0 && info != n+1 )
   {
      printf("ML_SuperLU_Solve : error coming from dgssvx %d\n", info);
      exit(1);
   } 
   else if ( solver->reuse_flag == 0 )
   {
#ifdef ML_DEBUG_SUPERLU
      if ( rcond != 0.0 && offset == 0 )
         printf("ML_SuperLU_Solve : condition number = %e\n", 1.0/rcond);
      else if ( offset == 0 )
         printf("ML_SuperLU_Solve : Recip. condition number = %e\n", rcond);
#endif
   }

   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_x[i+offset];

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   ML_memory_free( (void **) &local_x );
   ML_memory_free( (void **) &local_rhs );
   ML_memory_free( (void **) &etree );
   Destroy_SuperMatrix_Store(&B);
   Destroy_SuperMatrix_Store(&X);
   solver->reuse_flag = 1;
#elif DSUPERLU
   int                flag, N_local, offset;
   double            *local_rhs;
   ML_Comm           *mlcomm;
   ML_Solver         *solver;
   superlu_options_t  options;
   superlu_options_t *optionsptr;
   SuperMatrix       *A;
   ScalePermstruct_t *ScalePermstruct;
   SuperLUStat_t      stat;
   LUstruct_t        *LUstruct;
   gridinfo_t        *mygrid;
   ML_Lugrid         *lugrid_tiles;
   double             berr;
   int_t              i, n, nprow, npcol, nprocs; /* short integers */
   int                iam, info, nrhs, color, key;
   int                q, g, l, k, mygroup;
   int                stile, mtile, ltile, ntile, tsz, tsz2;
   int_t             *usermap;
   /* In-lining superlu_gridmap */
   MPI_Group mpi_base_group, superlu_grp;
   int                mycol, myrow, j;
   /* heap info arguments */
   int fragments, total_free, largest_free, total_used;

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen ) 
   {
      printf("ML_SuperLU_Solve error : lengths not matched.\n");
      exit(1);
   }
   solver          = (ML_Solver *) vsolver;
   A               = (SuperMatrix *) solver->Mat1;
   mlcomm          = (ML_Comm *) solver->void_params1;
   iam             = mlcomm->ML_mypid;
   nprocs          = mlcomm->ML_nprocs;
   N_local         = (int) solver->dble_params1[0];
   offset          = (int) solver->dble_params1[1];
   i               = (int) solver->dble_params1[2];
   n               = i;   /* n is a short integer */
   flag            = solver->reuse_flag;
   optionsptr      = &options;
   LUstruct        = (LUstruct_t *) solver->LUspl;
   ScalePermstruct = (ScalePermstruct_t *) solver->PERMspl;

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* Fetching the factorization (flag=1) is unneccessary           */
   /* ------------------------------------------------------------- */

   if ( flag == -999 ) 
   {
     if( iam == 0 )printf("ml_superlu: clean up\n"); /* dmd */

     /* deallocate storage and clean up */
     info = flag;
     if ( A != NULL ) 
     {
        ML_memory_free((void*)&(((NCformat *) A->Store)->rowind));
        ML_memory_free((void*)&(((NCformat *) A->Store)->colptr));
        ML_memory_free((void*)&(((NCformat *) A->Store)->nzval));
        SUPERLU_FREE( A->Store );
        ML_memory_free((void**) &A);
     }
     solver->Mat1 = NULL;
     mygroup = (int ) solver->ML_subgroup;
     lugrid_tiles    = solver->gridtiles;
     mygrid = &((lugrid_tiles[mygroup]).grid);
     Destroy_LU(n, mygrid, LUstruct);
     ScalePermstructFree(ScalePermstruct);
     LUstructFree(LUstruct);
     free(LUstruct);
     free(ScalePermstruct);
     solver->PERMspl = NULL;
     solver->LUspl = NULL;
     superlu_gridexit(mygrid);
     free (lugrid_tiles);
     solver->gridtiles = NULL;
     return 0;
   } 
   else if ( flag == 0 ) 
   {
      ML_SuperLU_Set_Tile(nprocs, &tsz ,&stile, &mtile, &ltile);
      ntile = stile + mtile + ltile;
      tsz2 = tsz * tsz;
      usermap = (int_t *) malloc( tsz2*sizeof(int_t) );
      lugrid_tiles = (ML_Lugrid *) malloc( ntile*sizeof(ML_Lugrid) );
      k = 0;
      for( g=0 ; g < ntile; g++)
      {
         if( g < stile )
         {
            tsz2 = (tsz-1)*(tsz-1);
            nprow = tsz-1;
            npcol = tsz-1;
         }
         else if( g < stile+mtile)
         {
            tsz2 = tsz*(tsz-1);
            nprow = tsz;
            npcol = tsz-1;
         }
         else
         {
            tsz2 = tsz*tsz;
            nprow = tsz;
            npcol = tsz;
         }
         for( l=0; l<tsz2; l++)
         {
            usermap[l] = l+k;
            if( iam == l+k ) mygroup = g;
         }
         k = k + tsz2;
         /* in-lining
          * superlu_gridmap( MPI_COMM_WORLD, 
          * nprow, npcol, usermap, nprow, &((lugrid_tiles[g]).grid));
          */
         (lugrid_tiles[g]).grid.nprow = nprow;
         (lugrid_tiles[g]).grid.npcol = npcol;
         MPI_Comm_group( MPI_COMM_WORLD, &mpi_base_group );
         MPI_Group_incl( mpi_base_group, tsz2, usermap, &superlu_grp );
         MPI_Comm_create(MPI_COMM_WORLD,superlu_grp,&(lugrid_tiles[g].grid.comm));
         if ( lugrid_tiles[g].grid.comm == MPI_COMM_NULL ) 
         {
            lugrid_tiles[g].grid.comm = MPI_COMM_WORLD;
            lugrid_tiles[g].grid.iam  = iam;
         } 
         else 
         {   /* mygroup=g and iam%mtile=npcol*myrow + mycol */
            MPI_Comm_rank(lugrid_tiles[g].grid.comm,&(lugrid_tiles[g].grid.iam));
            myrow = lugrid_tiles[g].grid.iam / npcol;
            mycol = lugrid_tiles[g].grid.iam % npcol;
            MPI_Comm_split(lugrid_tiles[g].grid.comm, 
                           myrow, mycol, &(lugrid_tiles[g].grid.rscp.comm));
            MPI_Comm_split(lugrid_tiles[g].grid.comm, 
                           mycol, myrow, &(lugrid_tiles[g].grid.cscp.comm));
            lugrid_tiles[g].grid.rscp.Np  = npcol;
            lugrid_tiles[g].grid.rscp.Iam = mycol;
            lugrid_tiles[g].grid.cscp.Np  = nprow;
            lugrid_tiles[g].grid.cscp.Iam = myrow;
         }
      } /* end for group g */
      if( nprocs != k )
      {
         printf("Error nprocs %d  k %d \n", nprocs, k);
         exit(-1);
      }
      free (usermap);
      solver->ML_subgroup = mygroup;
      /*
       * Fact = DOFACT Trans = NOTRANS Equil = EQUI RowPerm = LargeDiag
       * ColPerm = COLAMD ReplaceTinyPivot = REPLACE IterRefine = DOUBLE
       */
      set_default_options(optionsptr);
      optionsptr->Equil = NEQU;
      optionsptr->IterRefine = NOREFINE;
      /*
       * Future possiblities to experiment with include
       * optionsptr->RowPerm = NOROWPERM;
       * optionsptr->ColPerm = COLAMD;   (default)
       * optionsptr->ColPerm = MMD_AT_PLUS_A;
       * optionsptr->ColPerm = NATURAL;
       * ... and Equil
       */
      ScalePermstruct = ( ScalePermstruct_t *) malloc( sizeof( ScalePermstruct_t));
      ScalePermstructInit(n, n, ScalePermstruct);
      LUstruct = ( LUstruct_t *) malloc( sizeof( LUstruct_t) );
      LUstructInit(n, n, LUstruct);
      solver->PERMspl = (void *) ScalePermstruct;
      solver->LUspl = (void *) LUstruct;
      solver->gridtiles = lugrid_tiles;
      /* rst: mygrid is a pointer to a structure, not a structure.
       *  mygrid = ( gridinfo_t *) malloc( sizeof( gridinfo_t) );
       */
   } 
   else 
   {
      /* Indicate that the factored form of A is supplied. */
      /* Reset options */
      optionsptr->Fact = FACTORED;
      optionsptr->Trans = NOTRANS;
      optionsptr->Equil = NEQU;
      optionsptr->RowPerm = MY_PERMR;
      optionsptr->ColPerm = MY_PERMC;
      optionsptr->ReplaceTinyPivot = REPLACE;
      optionsptr->IterRefine = NOREFINE;
      lugrid_tiles    = solver->gridtiles;
      mygroup = (int ) solver->ML_subgroup;
   }
   mygrid = &((lugrid_tiles[mygroup]).grid);

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   nrhs = 1;
   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) mlcomm, local_rhs, &i, n);

   /* ------------------------------------------------------------- */
   /* perform LU decomposition and then solve                       */
   /* ------------------------------------------------------------- */

   info = flag;
   PStatInit(&stat);
   pdgssvx_ABglobal(optionsptr, A, ScalePermstruct, local_rhs, n, 
                    nrhs, mygrid, LUstruct, &berr, &stat, &info);

/*
if( iam == 0 ){
heap_info(&fragments, &total_free, &largest_free, &total_used);
printf("memory usage: fragments %d free: total %d, largest %d, total_used %d\n",
          fragments, total_free, largest_free, total_used);
}
*/

   if ( flag == 0 ) 
   {
     if ( A != NULL ) 
     {
/*
        ML_memory_free((void*)&(((NCformat *) A->Store)->rowind));
        ML_memory_free((void*)&(((NCformat *) A->Store)->colptr));
        ML_memory_free((void*)&(((NCformat *) A->Store)->nzval));
        SUPERLU_FREE( A->Store );
*/
        /* to satisfy pdgssvx_ABglobal argument check, postpone
         * ML_memory_free((void**) &A);
         */
      }
   }

   solver->reuse_flag = 1;
   PStatFree(&stat);
   if( info != 0 )
   {
      if( iam == 0 )printf("Error: ml_superlu    info = %d\n",info);
      return(-1);
   }
   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_rhs[i+offset];
   ML_memory_free( (void **) &local_rhs );
#else
   printf("ML_SuperLU_Solve : SuperLU not used.\n");
#endif
   return 0;
}

/* ******************************************************************** */
/* This subroutine calls the SuperLU subroutine to solve a given        */
/* subproblem where the matrix and right hand side are both residing    */
/* in the local processor (domain decomposition)                        */
/* ******************************************************************** */

int ML_SuperLU_SolveLocal(void *vsolver, double *x, double *rhs)
{
#ifdef ML_SUPERLU2

   SuperMatrix *A, *L, *U, B;
   int         i, n, info, flag, *perm_r, *perm_c, permc_spec;
   ML_Solver   *solver;

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   solver  = (ML_Solver *) vsolver;
   A       = (SuperMatrix *) solver->Mat1;
   n       = (int) solver->dble_params1[0];
   flag    = solver->reuse_flag;
   if ( flag != 0 ) 
   {
      for ( i = 0; i < n; i++ ) x[i] = rhs[i];
   }

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* ------------------------------------------------------------- */

   if ( flag == 0 ) 
   {
      ML_memory_alloc((void**) &perm_c, n * sizeof(int), "LU6" );
      ML_memory_alloc((void**) &perm_r, n * sizeof(int), "LU7" );
      solver->int_params1 = perm_r;
      solver->int_params2 = perm_c;
      permc_spec = 1;
      get_perm_c(permc_spec, A, perm_c);
      ML_memory_alloc((void**) &L, sizeof(SuperMatrix), "LU8" );
      ML_memory_alloc((void**) &U, sizeof(SuperMatrix), "LU9" );
      solver->Mat2 = (void *) L;
      solver->Mat3 = (void *) U;
   } 
   else 
   {
      perm_r = (int *) solver->int_params1;
      perm_c = (int *) solver->int_params2;
      L = (SuperMatrix *) solver->Mat2;
      U = (SuperMatrix *) solver->Mat3;
   }

   /* ------------------------------------------------------------- */
   /* create right hand side vector conforming to SuperLU format    */
   /* ------------------------------------------------------------- */

   dCreate_Dense_Matrix(&B, n, 1, x, n, DN, _D, GE);
   info = flag;
   dgssv(A, perm_c, perm_r, L, U, &B, &info);

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   Destroy_SuperMatrix_Store(&B);
   solver->reuse_flag = 1;
#else
   printf("ML_SuperLU_SolveLocal : SuperLU not used.\n");
#endif
   return 0;
}

/* ************************************************************************* */
/* clean up                                                                  */
/* ------------------------------------------------------------------------- */

int ML_CSolve_Clean_SuperLU( void *vsolver, ML_CSolveFunc *func)
{
   ML_Solver   *solver;

#ifdef SUPERLU
   SuperMatrix *Amat;

   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);

   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL )
   {
      SUPERLU_FREE( ((NRformat *) Amat->Store)->colind);
      SUPERLU_FREE( ((NRformat *) Amat->Store)->rowptr);
      SUPERLU_FREE( ((NRformat *) Amat->Store)->nzval);
      SUPERLU_FREE( Amat->Store );
      ML_memory_free(  (void**) &(solver->Mat1) );
      solver->Mat1 = NULL;
   }
#elif DSUPERLU
   SuperMatrix *Amat;

   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);
   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL)
   {
      Destroy_CompCol_Matrix(Amat);
      ML_memory_free((void**) &Amat);
   }
   solver->Mat1 = NULL;

#else
   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);
#endif
   ML_Solver_Destroy( &solver );
   return 0;
}

/* ******************************************************************** */
/* An array of processors decomposes into tiles of three sizes:         */
/* (n-1)^2 , n (n-1) or n^2.  * For any natural numbers, p and n,       */
/* there exist integers i,j,k such that 0<=i,j<n and                    */
/* p = i*(n-1)^2 + j*n*(n-1) + k*n^2 = (i+j+k)*(n^2) -n*(2*i+j)+i       */
/* If p >= (2*n-1)*(n-1)^2 , then k >= 0                                */
/* ******************************************************************** */

void ML_SuperLU_Set_Tile( int p, int* n, int* i, int* j, int* k)
{
   int l,q,r,s;
   double quotient, cuberoot; 
   if( p < 12) 
   {
      *n = 2;
      *i = p;
      *j = 0;
      *k = 0;
   }
   else if( p < 54) 
   {
      /* p = l + 4s + 12 q */
      l = p % 4;
      r = (p-l)/4;
      s = r % 3;
      q = (r-s)/3;
      if( l == 0 )
      {
         *n = 2;
         *i = 0;
         *j = 0;
         *k = s + 3*q;
      }
      else if( l == 1 )
      {
         *n = 3;
         *i = s + 3*q -2;
         *j = 0;
         *k = 1; 
      }
      else if( l == 2 )
      {
         *n = 3;
         *i = s + 3*q -1;
         *j = 1;
         *k = 0; 
      }
      else
      { 
         *n = 3;
         *i = s + 3*q -3;
         *j = 1;
         *k = 1; 
      }
   }
   else
   {
      quotient = (double) p / 2;
      cuberoot = (double) 1 / 3;
      *n = (int) floor(pow(quotient,cuberoot));
      *i = p % (*n);
      q = (p - *i )/ *n;
      r = (q + 2 * *i ) % (*n);
      if( r > 0 )
         *j = *n - r;
      else
         *j = 0;
      *k = ( (q + 2 * *i + *j ) / (*n) ) - *i - *j;
   }
} /* end of ML_SuperLU_Set_Tile */

