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
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "phg.h"

#ifdef ZOLTAN_PARKWAY

void Zoltan_ParaPartKway(int numVertices, int numHedges, const int *vWeights, const int *hEdgeWts, const int *pinList, const int *offsets, int numParts, double constraint, int *k_1cut, const int *options, int *pVector, const char *outFile, MPI_Comm comm);    
    
static int scale_round_weights(float *, int *, int, int, int);

#define ZOLTAN_PARKWAY_ERROR(str, err) \
  {ZOLTAN_PRINT_ERROR(zz->Proc, yo, str); ierr = err; goto End;}

#endif  /* ZOLTAN_PARKWAY */


/*****************************************************************************/
int Zoltan_PHG_ParKway(
  ZZ        *zz,
  HGraph    *hg,
  int       nparts,           /* # of desired partitions */
  Partition partvec,          /* Output:  partition assignment vector */
  PHGPartParams *hgp          /* Input: hypergraph parameters */  
)
{
    int ierr = ZOLTAN_OK;
    char *yo = "Zoltan_HG_ParKway";

#ifndef ZOLTAN_PARKWAY
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "ParKway method selected but Zoltan is not"
                       "built and linked with ParKway.");
    return ZOLTAN_FATAL;
#else

    int options[26];                      /* ParKway options */
    int *ivwgts = NULL, *iewgts = NULL;   /* ParKway expects integer weights. */
    int *pvector=NULL;                    /* partvec for "local" vertices */
    int cut;                              /* Diagnostics from ParKway */
    double constraint;                    /* imbalance ratio */
    int i, anVtx, nVtx;                   /* counter and local vertex cnt (for first k-1 parts)*/
    PHGComm *hgc=hg->comm;
    int *disp=NULL, *recv_size=NULL;      /* for allgatherv */
    
    /* ParKway expects integer weights; convert if weights are provided. */
    ivwgts = (int *) ZOLTAN_MALLOC(hg->nVtx  * sizeof(int));
    iewgts = (int *) ZOLTAN_MALLOC(hg->nEdge * sizeof(int));    
    if (!ivwgts || !iewgts)
        ZOLTAN_PARKWAY_ERROR("Memory error.", ZOLTAN_MEMERR);
    
    
    if (hg->VtxWeightDim > 1) { 
        ZOLTAN_PARKWAY_ERROR("ParKway supports Vtx_Weight_Dim == 0 or 1 only.",
                             ZOLTAN_FATAL);
    } else if (hg->VtxWeightDim == 1) 
        scale_round_weights(hg->vwgt, ivwgts, hg->nVtx, hg->VtxWeightDim, 0);
    else 
        for (i=0; i<hg->nVtx; ++i)
            ivwgts[i] = 1;

        
    if (hg->EdgeWeightDim > 1) {
        ZOLTAN_PARKWAY_ERROR("ParKway supports Edge_Weight_Dim == 0 or 1 only.",
                             ZOLTAN_FATAL);
    } else if (hg->EdgeWeightDim == 1) 
        scale_round_weights(hg->ewgt, iewgts, hg->nEdge, hg->EdgeWeightDim, 0);
    else
        for (i=0; i<hg->nEdge; ++i)
            iewgts[i] = 1;
    
    
    
    anVtx = hg->nVtx / hgc->nProc;
    nVtx = (hgc->myProc==hgc->nProc-1) ? hg->nVtx-(anVtx*(hgc->nProc-1)) : anVtx;
    uprintf(hgc, "anVtx=%d nVtx=%d nProc=%d\n", anVtx, nVtx, hgc->nProc);

    pvector = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int));
    disp = (int *) ZOLTAN_MALLOC(hgc->nProc * sizeof(int));
    recv_size = (int *) ZOLTAN_MALLOC(hgc->nProc * sizeof(int));
    if ((nVtx && !pvector) || !disp || !recv_size)
        ZOLTAN_PARKWAY_ERROR("Memory error.", ZOLTAN_MEMERR);


    /* ----- Set ParKway's options --------------- */
    options[0] = 1; /*0 -> all options use default, else user define*/
    options[1] = 0; /*0 -> use default seed (or if using SPRNG, it chooses seed), else use options[1] as seed*/
    options[2] = 2; /*0 -> no disp info, 1 -> some, 2 -> lots*/
    options[3] = 0; /*0 -> do not write partition to disk, 1 -> do write*/
    options[4] = 1; /*number of parallel runs*/
    options[5] = 0; /*perform random vertex shuffle: 0 -> no, 1 -> yes*/
    options[6] = 200; /*numParts*options[5] -> min number of vertices in coarse hypergraph*/
    options[7] = 7; /*[7] and [8] specify reduction ratio in parallel coasrsening*/
    options[8] = 4; /*r = [7]/[8]*/
    options[9] = 3; /*vertex visit order: 3 -> random, 1/2 inc/dec by vertex id, 4/5 inc/dec by vertex wt*/
    options[10] = 1; /*use cluster weight to determine connectivity: 0 -> no, 1 -> yes  */
    options[11] = 3; /*matching request resolution order: 3 -> random, 2 -> as they arrive */
    options[12] = 1; /*number serial partitioning runs*/

    options[13] = 1; /* serial partitioning routine, 1-3 RB, 4 khmetis, 5 patoh, see manual*/
    if (!strcasecmp(hgp->parkway_serpart, "patoh"))
        options[13] = 5;
    else if (!strcasecmp(hgp->parkway_serpart, "hmetis"))
        options[13] = 4;

    options[14] = 2; /*serial coarsening algorithm (only if [13] = RB, see manual)*/
    options[15] = 2; /*num bisection runs in RB (only if [13] = RB, see manual)*/
    options[16] = 10; /*num initial partitioning runs in RB (only if [13] = RB, see manual)*/
    options[17] = 2; /*hmetis_PartKway coarsening option, vals 1-5, see manual (only if [13] = 4)  */ 
    options[18] = 2; /*hmetis_PartKway refinement option, vals 0-3, see manual (only if [13] = 4)*/
    options[19] = 1; /*patoh_partition parameter settings, vals 1-3, see manual (only if [13] = 5)*/
    options[20] = 1; /*parallel uncoarsening algorithm, 1 simple, 2 only final V-Cycle, 3 all V-Cycle*/
    options[21] = 1000; /*limit on number of V-Cycle iterations (only if [21] = 2/3)*/
    options[22] = 0; /*min allowed gain for V-Cycle (percentage, see manual, only if [21] = 2/3) */
    options[23] = 70; /*percentage threshold used to reject partitions from a number of runs (see manual) */
    options[24] = 70; /*reduction in [23] as partitions propagate by factor [24]/100 (see manual) */
    options[25] = 100; /*early exit criterion in parallel refinement, will exit if see ([25]*num vert)/100 consecutive -ve moves */
    
    constraint = hgp->bal_tol-1.0;

    Zoltan_Print_Sync_Start (zz->Communicator, 1);
    uprintf(hgc, "nVtx=%d nEdge=%d nparts=%d constrain=%.3lf\n", nVtx, hg->nEdge, nparts, constraint);
    for (i=0; i<hg->nEdge; ++i) {
        int j;

        uprintf(hgc, "%d (%d): ", i, iewgts[i]);
        for (j=hg->hindex[i]; j<hg->hindex[i+1]; ++j)
            printf("%d (%d), ", hg->hvertex[j], ivwgts[hg->hvertex[j]]);
        printf("\n");
    }
    
    Zoltan_Print_Sync_End(zz->Communicator, 1);
      
    Zoltan_ParaPartKway(nVtx, hg->nEdge, &ivwgts[hgc->myProc*anVtx], iewgts,
                 hg->hindex, hg->hvertex, nparts,
                 constraint, &cut, options, pvector, NULL, zz->Communicator);

    Zoltan_Print_Sync_Start (zz->Communicator, 1);    
    uprintf(hgc, "ParaPartKway cut=%d\n pvector:\n", cut);


    for (i=0; i<nVtx; ++i) {
        uprintf(hgc, "%d -> %d\n", i, pvector[i]);
        if (pvector[i]<0 || pvector[i]>=nparts)
            errexit("parkway: invalid pvector[%d]=%d", i, pvector[i]); 
    }
    Zoltan_Print_Sync_End(zz->Communicator, 1);
    
    /* after partitioning Zoltan needs partvec exist on all procs for nProc_x=1 */       
    disp[0] = 0; 
    for (i = 1; i < hgc->nProc; ++i)
        disp[i] = disp[i-1] + anVtx;

    memset(partvec, 0xff, sizeof(int)*hg->nVtx);

    
    MPI_Allgatherv(pvector, 1, MPI_INT, 
                  partvec, recv_size, disp, MPI_INT, zz->Communicator);

    /*
    for (i = 0; i < hgc->nProc; ++i)
        if (recv_size[i]!=anVtx) {
            errexit("recv_size[%d]=%d but anVtx=%d\n", i, recv_size[i], anVtx);
        }
    if (recv_size[i]!=nVtx) {
        errexit("recv_size[%d]=%d but nVtx=%d\n", i, recv_size[i], nVtx);
    }
    */

    Zoltan_Print_Sync_Start (zz->Communicator, 1);    
    uprintf(hgc, "partvec:\n");


    for (i=0; i<hg->nVtx; ++i) {
        uprintf(hgc, "%d -> %d\n", i, partvec[i]);
        if (partvec[i]<0 || partvec[i]>=nparts)
            errexit("parkway: invalid partvec[%d]=%d", i, partvec[i]);
    }

    Zoltan_Print_Sync_End(zz->Communicator, 1);
    
    for (i=0; i<hg->nVtx; ++i)
        partvec[i] = 0;
            
    
  /* HERE:  Check whether imbalance criteria were met. */

End:

    Zoltan_Multifree(__FILE__,__LINE__, 5, &ivwgts, &iewgts, &pvector, &disp, &recv_size);
    uprintf(hgc, "Terminating parKway...\n");
    
#endif
  return ierr;
}

    
/*****************************************************************************/

#ifdef ZOLTAN_PARKWAY

    
#define INT_EPSILON (1e-5)

static int scale_round_weights(
  float *fwgts, 
  int *iwgts, 
  int n, 
  int dim,
  int mode
)
{
/* Convert floating point weights to integer weights.
 * This routine is stolen from scale_round_weights in parmetis_jostle.c.
 * Because it needs to run only serially, and because it uses only 
 * integers (not idxtype), it has been largely duplicated here.
 */

  int i, j, tmp, ierr; 
  int max_wgt_sum = INT_MAX/8;
  int *nonint;
  float *scale, *sum_wgt, *max_wgt;
  char msg[256];
  static char *yo = "scale_round_weights";

  ierr = ZOLTAN_OK;

  if (mode == 0) {
    /* No scaling; just convert to int */
    for (i=0; i<n*dim; i++){
      iwgts[i] = (int) ceil((double) fwgts[i]);
    }
  }
  else{
    /* Allocate local arrays */
    nonint = (int *)ZOLTAN_MALLOC(dim*sizeof(int));
    scale = (float *)ZOLTAN_MALLOC(3*dim*sizeof(float));
    sum_wgt = scale + dim;
    max_wgt = sum_wgt + dim;
    if (!(nonint && scale)){
      ZOLTAN_PRINT_ERROR(0, yo, "Out of memory.");
      ZOLTAN_FREE(&nonint);
      ZOLTAN_FREE(&scale);
      return ZOLTAN_MEMERR;
    }

    /* Initialize */
    for (j=0; j<dim; j++){
      nonint[j] = 0;
      sum_wgt[j] = 0;
      max_wgt[j] = 0;
    }

    /* Compute local sums of the weights */
    /* Check if all weights are integers */
    for (i=0; i<n; i++){
      for (j=0; j<dim; j++){
        if (!nonint[j]){ 
          /* tmp = (int) roundf(fwgts[i]);  EB: Valid C99, but not C89 */
          tmp = (int) floor((double) fwgts[i] + .5); /* Nearest int */
          if (fabs((double)tmp-fwgts[i*dim+j]) > INT_EPSILON){
            nonint[j] = 1;
          }
        }
        sum_wgt[j] += fwgts[i*dim+j];
        if (fwgts[i*dim+j] > max_wgt[j])
          max_wgt[j] = fwgts[i*dim+j]; 
      }
    }

    /* Calculate scale factor */
    for (j=0; j<dim; j++){
      scale[j] = 1.;
      /* Scale unless all weights are integers (not all zero) */
      if (nonint[j] || (max_wgt[j] <= INT_EPSILON) 
                    || (sum_wgt[j] > max_wgt_sum)){
        if (sum_wgt[j] == 0){
          ierr = ZOLTAN_WARN;
          sprintf(msg, "All weights are zero in component %1d", j);
          ZOLTAN_PRINT_WARN(0, yo, msg);
        }
        else /* sum_wgt[j] != 0) */
          scale[j] = max_wgt_sum/sum_wgt[j];
      }
    }

    /* If mode==2, let the scale factor be the same for all weights */
    if (mode==2){
      for (j=1; j<dim; j++){
        if (scale[j]<scale[0])
          scale[0] = scale[j];
      }
      for (j=1; j<dim; j++){
        scale[j] = scale[0];
      }
    }

    /* Convert weights to positive integers using the computed scale factor */
    for (i=0; i<n; i++){
      for (j=0; j<dim; j++){
        iwgts[i*dim+j] = (int) ceil((double) fwgts[i*dim+j]*scale[j]);
      }
    }

    ZOLTAN_FREE(&nonint);
    ZOLTAN_FREE(&scale);
  }
  return ierr;
}


#endif /* ZOLTAN_PARKWAY */


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
