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

#include <stdlib.h>
#include "phg.h"

#undef  USE_SUBROUNDS

static ZOLTAN_PHG_MATCHING_FN pmatching_local; /* function for local matching */
static ZOLTAN_PHG_MATCHING_FN pmatching_ipm;   /* inner product matching */
static ZOLTAN_PHG_MATCHING_FN pmatching_alt_ipm;   /* alternating ipm */

static void phasethreereduce (void*, void*, int*, MPI_Datatype*);
static int Zoltan_PHG_match_isolated(ZZ* zz, HGraph* hg, Matching match, 
                                     int small_degree);

typedef struct triplet {
    int candidate;      /* gno of candidate vertex */
    int partner;        /* gno of best match found so far */
    float ip;           /* total inner product between candidate and partner */
} Triplet;

static Triplet *Tmp_Best = NULL;  /* Temp buf used in MPI_Allreduce fn */
static HGraph *HG_Ptr;



/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist=1;
    
    if (!strcasecmp(hgp->redm_str, "no"))
        hgp->matching = NULL;
    else if (!strncasecmp(hgp->redm_str, "l-", 2))  {
        HGPartParams hp;

        strcpy(hp.redm_str, hgp->redm_str+2);
        strcpy(hp.redmo_str, hgp->redmo_str);
        if (!Zoltan_HG_Set_Matching_Fn(&hp)) {
            exist = 0;
            hgp->matching = NULL;
        } else {   
            hgp->matching = pmatching_local; 
            hgp->locmatching = hp.matching;
            hgp->matching_opt = hp.matching_opt;
        }
    } else if (!strcasecmp(hgp->redm_str, "c-ipm"))
        hgp->matching = pmatching_ipm;   
    else if (!strcasecmp(hgp->redm_str, "ipm"))
        hgp->matching = pmatching_ipm;
    else if (!strcasecmp(hgp->redm_str, "alt-ipm"))
        hgp->matching = pmatching_alt_ipm;
    else {
        exist = 0;
        hgp->matching = NULL;
    }
    
    return exist;
}


/*****************************************************************************/
int Zoltan_PHG_Matching (
  ZZ *zz,
  HGraph *hg,
  Matching match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   ierr = ZOLTAN_OK;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hgp->edge_scaling) {
     if (hg->nEdge && !(new_ewgt = (float*) 
                      ZOLTAN_MALLOC(hg->nEdge * sizeof(float))))
         MEMORY_ERROR;
 
     Zoltan_PHG_Scale_Edges (zz, hg, new_ewgt, hgp);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
  }

  /* Create/update scale vector for vertices for inner product */
  if (hgp->vtx_scaling) 
     Zoltan_PHG_Scale_Vtx (zz, hg, hgp);
  
  
  /* Do the matching */
  if (hgp->matching) {
    /* first match isolated vertices */
    Zoltan_PHG_match_isolated(zz, hg, match, 0);
    /* now do the real matching */
    ierr = hgp->matching (zz, hg, match, hgp);
    /* clean up by matching "near-isolated" vertices of degree 1 */
    /* only useful in special cases (e.g. near-diagonal matrices).
       in many cases there is a slight increase in cuts, so 
       turn it off for now.                                      */
    /* Zoltan_PHG_match_isolated(zz, hg, match, 1); */
  }

End: 

  /* Restore the old edge weights if scaling was used. */
  if (hgp->edge_scaling)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return ierr;
}


static int Zoltan_PHG_match_isolated(
  ZZ *zz,
  HGraph *hg,
  Matching match,
  int small_degree /* 0 or 1; 0 corresponds to truely isolated vertices */
)
{
    int v=-1, i, unmatched=0, *ldeg, *deg;
#ifdef _DEBUG
    int cnt=0;
#endif
    static char *yo = "Zoltan_PHG_match_isolated";
    int ierr = ZOLTAN_OK;

    if (hg->nVtx) {
      if (!(ldeg = (int*)  ZOLTAN_MALLOC(2*hg->nVtx*sizeof(int))))
        MEMORY_ERROR;
      deg = ldeg + hg->nVtx;
      /* match isolated vertices.
         UVCUVC: right now we match in the natural order,
         I don't think we need random matching but if needed
         we can match in random order. */
      for (i=0; i<hg->nVtx; ++i)
          ldeg[i] = hg->vindex[i+1] - hg->vindex[i];
      MPI_Allreduce(ldeg, deg, hg->nVtx, MPI_INT, MPI_SUM, hg->comm->col_comm);
      
      if (small_degree>0){
          /* Only match on procs with many unmatched vertices */
          unmatched= 0;
          for (i=0; i<hg->nVtx; ++i)
              if (match[i]==i) unmatched++;
      }
      if ((small_degree==0) || (unmatched > 0.8*hg->nVtx))
          for (i=0; i<hg->nVtx; ++i){
              if ((match[i]==i) && (deg[i] <= small_degree)) { 
#ifdef _DEBUG
                  ++cnt;
#endif
                  /* match with previous unmatched vertex */
                  /* EBEB For degree-1 vertices, we could be more clever
                     and match vertices that share a common neighbor */
                  if (v==-1)
                      v = i;
                  else {
                      match[v] = i;
                      match[i] = v;
                      v = -1;
                  }
              }
          }
#ifdef _DEBUG
      if (cnt)
          uprintf(hg->comm, "Local H(%d, %d, %d) and there were %d isolated vertices\n", hg->nVtx, hg->nEdge, hg->nPins, cnt);           
#endif
End:
      ZOLTAN_FREE(&ldeg);
    }
    return ierr;
}

static int pmatching_local(
  ZZ *zz,
  HGraph *hg,
  Matching match,
  PHGPartParams *hgp
)
{
    int limit=hg->nVtx, err=ZOLTAN_OK;
    PHGComm *hgc=hg->comm;
    int root_matchcnt, root_rank;
    
    err = hgp->locmatching (zz, hg, match, &limit);
    
    /* Optimization */
    if (hgp->matching_opt) 
        err = hgp->matching_opt (zz, hg, match, &limit);
    
    /* find the index of the proc in column group with the best match
       (max #matches); it will be our root proc */
    Zoltan_PHG_Find_Root(hg->nVtx-limit, hgc->myProc_y, hgc->col_comm,
                         &root_matchcnt, &root_rank);
    
    MPI_Bcast(match, hg->nVtx, MPI_INT, root_rank, hgc->col_comm);
    
    return err;
}


/**************************************************************************
  Alternating ipm method. Alternate between full ipm and fast method. 
 *************************************************************************/
static int pmatching_alt_ipm(
  ZZ *zz,
  HGraph* hg,
  Matching match,
  PHGPartParams *hgp
)
{
  int ierr = ZOLTAN_OK;
  char redm_orig[MAX_PARAM_STRING_LEN];
  static int level=0;
  static int old_nvtx=0;

  strcpy(redm_orig, hgp->redm_str); /* save original parameter string */

  if (hg->nVtx > old_nvtx){
    /* larger hgraph; must have started new bisection v-cycle */
    level= 0;
  }

  /* first level is 0 */
  if ((level&1) == 0)  /* alternate even-odd levels */
    strcpy(hgp->redm_str, hgp->redm_fast); /* fast method is c-ipm for now */
  else
    strcpy(hgp->redm_str, "ipm");  

  ierr = pmatching_ipm(zz, hg, match, hgp);  /* only works for ipm and c-ipm! */

  ++level;  /* we don't have access to level data, so keep track this way */
  old_nvtx = hg->nVtx;

  /* set redm parameter back to original */
  strcpy(hgp->redm_str, redm_orig);
  
  return ierr;
}

/****************************************************************************
 * inner product matching (with user selectable column variant, c-ipm)
 * Parallelized version of the serial algorithm (see hg_match.c)
 * Based on conversations with Rob Bisseling by Aaron Becker, UIUC, summer 2004
 * completed by R. Heaphy
 */
               
#define ROUNDS_CONSTANT 8     /* controls the number of candidate vertices */ 
#define IPM_TAG        28731  /* MPI message tag, arbitrary value */
#define HEADER_COUNT    4          /* Phase 2 send buffer header size in ints */

/* these thresholds need to become parameters in phg - maybe ??? */
#define PSUM_THRESHOLD 0.0    /* ignore inner products (i.p.) < threshold */
#define TSUM_THRESHOLD 0.0    /* ignore inner products (i.p.) < threshold */
                                 
/* Forward declaration for a routine that encapsulates the common calls to use
** the Zoltan unstructured communications library for the matching code */
static int communication_by_plan (ZZ* zz, int sendcnt, int* dest, int* size, 
 int scale, int* send, int* reccnt, int* recsize, int* nRec, int** rec,
 MPI_Comm comm, int tag);

/* Actual inner product calculations between candidates (in rec buffer) */
/* and local vertices.  Not inlined because inline is not universal, yet.  */
#define INNER_PRODUCT1(ARG)\
  for (i = 0; i < count; r++, i++)\
    for (j = hg->hindex[*r]; j < hg->hindex[*r + 1]; j++) {\
      if (cmatch[hg->hvertex[j]] == hg->hvertex[j])  {\
        if (sums[hg->hvertex[j]] == 0.0)\
          index[m++] = hg->hvertex[j];\
        sums[hg->hvertex[j]] += (ARG);\
      }\
    }
            
/* Mostly identical inner product calculation to above for c-ipm variant. Here */
/* candidates are a subset of local vertices and are not in a separate buffer  */     
#define INNER_PRODUCT2(ARG)\
   for (i = hg->vindex[candidate_gno]; i < hg->vindex[candidate_gno+1]; i++)  {\
     edge = hg->vedge[i];\
     for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++)  {\
       if (cmatch[hg->hvertex[j]] == hg->hvertex[j])    {\
         if (sums[hg->hvertex[j]] == 0.0)\
           index[m++] = hg->hvertex[j];\
         sums[hg->hvertex[j]] += (ARG);\
       }\
     }\
   }  

/* simple macro to start timer */
#define MACRO_TIMER_START(arg, message, sync) \
  if (hgp->use_timers > 3)  {\
    if (timer[arg] < (arg))\
      timer[arg] = Zoltan_Timer_Init(zz->ZTime, sync, message);\
    ZOLTAN_TIMER_START(zz->ZTime, timer[arg], hg->comm->Communicator);\
  }   

/* simple corresponding macro to stop timer */
#define MACRO_TIMER_STOP(arg) \
  if (hgp->use_timers > 3) \
    ZOLTAN_TIMER_STOP(zz->ZTime, timer[arg], hg->comm->Communicator);

/* convenience macro to encapsulate resizing a buffer when necessary */
#define MACRO_REALLOC(new_size, old_size, buffer)  {\
  old_size = new_size;\
  if (!(buffer = (int*) ZOLTAN_REALLOC (buffer, old_size * sizeof(int)))) {\
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.");\
    err = ZOLTAN_MEMERR;\
    goto fini;\
    }\
  } 

  
/****************************************************************************/
/* Because this calculation is done in two locations it has been converted to a
** subroutine to assure consistancy. Inline is not yet always available!
** ERIK: need to calculate nCandidates based on # of unmatched vertices     */
static int calc_nCandidates (int num_vtx, int procs)
{
  /* Constant 2 below because each match pairs 2 vertices */
  return num_vtx ? 1 + num_vtx/(2 * procs * ROUNDS_CONSTANT) : 0;
}
 

/****************************************************************************/
static int pmatching_ipm (ZZ *zz,
  HGraph* hg,
  Matching match,
  PHGPartParams *hgp)
{
  int i, j, k, n, m, round, vindex, kstart, *r, *s;        /* loop counters  */
  int lno, count, old_kstart;                              /* temp variables */
  int candidate_gno;                             /* gno of current candidate */
  int sendcnt, sendsize, reccnt, recsize, msgsize;         /* temp variables */
  int nRounds;                /* # of matching rounds to be performed;       */
                              /* identical on all procs in hgc->Communicator.*/
  int nCandidates;            /* # of candidates on this proc; identical     */
                              /* on all procs in hgc->col_comm.              */
  int nTotal;                 /* on a given proc, total # of candidates for
                               which to compute inner products. When using 
                               SPARSE_CANDIDATES, nTotal may differ among procs
                               in a column; otherwise, nTotal is the same on
                               all procs in the communicator.  */
  int max_nTotal;            /* max within proc column of nTotal. */
  int total_nCandidates;     /* Sum of nCandidates across row. */
  int *send = NULL,          /* working buffers, may be reused. */
      *dest = NULL,
      *size = NULL,
      *rec = NULL,
      *index = NULL,
      *aux = NULL;
  int *visit = NULL,       /* fixed usage arrays, candidate visit order */
      *cmatch = NULL,      /* working copy of match array */
      *select = NULL,      /* current selected candidates */
      *permute = NULL,    /* reorder of candidates after global communication */
      *edgebuf = NULL;    /* holds received candidates for processing */
  int nSend, nDest, nSize, nRec, nIndex, nEdgebuf; /* currently allocated size
                                                      of the corresponding
                                                      working buffers */
  float bestsum;      /* holds current best inner product */
  float *sums = NULL, /* holds candidate's inner products with each local vtx */
        *f = NULL;    /* used to stuff floating value into integer message */
  PHGComm *hgc = hg->comm;
  int err = ZOLTAN_OK, old_row, row;
  int max_nPins, max_nVtx;       /* Global max # pins/proc and vtx/proc */
  int **rows = NULL;             /* used only in merging process */
  int bestlno, partner_gno, nselect, edge;
  Triplet *master_data = NULL;
  int *master_procs = NULL;
  int cFLAG;                    /* if set, do only a column matching, c-ipm */
  Triplet *global_best = NULL;
  MPI_Op phasethreeop;
  MPI_Datatype phasethreetype;
  int candidate_index;
  int first_candidate_index;
  static int timer[7] = {-1, -1, -1, -1, -1, -1, -1};
  char *yo = "pmatching_ipm";
  
  
  ZOLTAN_TRACE_ENTER (zz, yo);
  MACRO_TIMER_START (0, "matching setup", 0);
     
  /* this restriction may be removed later, but for now NOTE this test */
  if (sizeof(int) < sizeof(float))  {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Code must be modified before using");
    err = ZOLTAN_FATAL;
    goto fini;
  }

  /* set a flag if user wants a column matching or a full matching */
  cFLAG = strcasecmp(hgp->redm_str, "c-ipm") ? 0 : 1;
  if (!cFLAG) {
    MPI_Type_contiguous (sizeof(Triplet), MPI_CHAR, &phasethreetype);
    MPI_Type_commit (&phasethreetype);
    MPI_Op_create (&phasethreereduce, 1, &phasethreeop);
  }

  /* determine basic working parameters */
  nRounds     = cFLAG ? ROUNDS_CONSTANT : hgc->nProc_x * ROUNDS_CONSTANT;
  nCandidates = calc_nCandidates (hg->nVtx, cFLAG ? 1 : hgc->nProc_x); 
    
  /* determine maximum number of Vtx and Pins for storage allocation */
  /* determine initial sum of all candidates = total_nCandidates==>allocation */
  if (cFLAG)  {
    total_nCandidates = nCandidates;
    max_nVtx  = hg->nVtx;
    max_nPins = hg->nPins;
  }
  else  {
    MPI_Allreduce(&hg->nPins, &max_nPins, 1, MPI_INT,MPI_MAX,hgc->Communicator);
    max_nVtx = total_nCandidates = 0;
    for (i = 0; i < hgc->nProc_x; i++)  {
      count = hg->dist_x[i+1]-hg->dist_x[i];  /* number of vertices on proc i */
      if (count > max_nVtx)
        max_nVtx = count;
      if (i == hgc->myProc_x)
        first_candidate_index = total_nCandidates;
      total_nCandidates += calc_nCandidates (count, hgc->nProc_x);
    }
  }
                 
  /* allocate "complicated" fixed sized array storage */
  nIndex = 1 + MAX(MAX(total_nCandidates, max_nVtx), hgc->nProc_y);
  nDest  = 1 + MAX(MAX(hgc->nProc_x,hgc->nProc_y),
                   MAX(total_nCandidates,max_nVtx));
  nSize  = 1 + MAX(MAX(hgc->nProc_x,hgc->nProc_y),
                   MAX(total_nCandidates,max_nVtx));

  /* These 3 buffers are REALLOC'd iff necessary; this should be very rare  */
  nSend    = max_nPins;   /* nSend/nEdgebuf are used for candidate exchange */
  nRec     = max_nPins;  
  nEdgebuf = max_nPins;   /* list of <candidate_gno, #pins, pin_list> */

  if (hg->nVtx)  
    if (!(cmatch = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
     || !(visit  = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
     || !(aux    = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))     
     || !(sums   = (float*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(float)))) {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.");
       err = ZOLTAN_MEMERR;
       goto fini;
    }

  if (!cFLAG && total_nCandidates && (hgc->myProc_y == 0))   /* Master row */
    if (!(master_data=(Triplet*)ZOLTAN_MALLOC(total_nCandidates*sizeof(Triplet))) ||
     !(global_best=(Triplet*)ZOLTAN_MALLOC(total_nCandidates*sizeof(Triplet)))){
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.");
       err = ZOLTAN_MEMERR;
       goto fini;
    }  

  if (!cFLAG && hgc->myProc_y == 0)
    for (i = 0; i < total_nCandidates; i++) {
      master_data[i].candidate = -1;
      master_data[i].partner   = -1;
      master_data[i].ip        = -1.0;
    }

  if (!cFLAG)
    if (!(edgebuf = (int*) ZOLTAN_MALLOC(nEdgebuf   * sizeof(int)))) {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.");
       err = ZOLTAN_MEMERR;
       goto fini;
    }  
  
  if ((total_nCandidates 
   && !(permute = (int*) ZOLTAN_MALLOC(total_nCandidates * sizeof(int))))
   || (nCandidates && !(select = (int*) ZOLTAN_MALLOC(nCandidates*sizeof(int))))
   || (nSend && !(send = (int*) ZOLTAN_MALLOC(nSend * sizeof(int))))
   || !(dest = (int*) ZOLTAN_MALLOC(nDest * sizeof(int)))
   || !(size = (int*) ZOLTAN_MALLOC(nSize * sizeof(int)))
   || (nRec && !(rec = (int*) ZOLTAN_MALLOC(nRec * sizeof(int))))
   || !(index = (int*)  ZOLTAN_MALLOC(nIndex * sizeof(int)))
   || !(rows = (int**) ZOLTAN_MALLOC((hgc->nProc_y + 1) * sizeof(int*)))) {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.");
     err = ZOLTAN_MEMERR;
     goto fini;
  }
    
  /* match[] is local slice of global matching array.  It uses local numbering 
   * (zero-based). Initially, match[i] = i. After matching, match[i]=i indicates
   * an unmatched vertex. A matching between vertices i & j is indicated by 
   * match[i] = j & match[j] = i.  NOTE: a match to an off processor vertex is
   * indicated by a negative number, -(gno+1), and must use global numbers
   * (gno's) and not local numbers, lno's, which are zero based per processor */
  /* Compute candidates' vertex visit order (selection). Random is default. */
  Zoltan_PHG_Vertex_Visit_Order (zz, hg, hgp, visit);
  
  /* Loop processing ncandidates vertices per column each round.
   * Each loop has 4 phases:
   * Phase 1: send ncandidates vertices for global matching - horizontal comm
   * Phase 2: sum  inner products, find best in column - vertical communication
   * Phase 3: return best sums to owner's column    - horizontal communication
   * Phase 4: return actual match selections        - horizontal communication
   *
   * No conflict resolution required because temp locking prevents conflicts. */

  MACRO_TIMER_STOP (0);
  vindex = 0;                                /* marks position in visit array */
  for (round = 0; round < nRounds; round++) {
    MACRO_TIMER_START (1, "matching phase 1", 0);
    
    /************************ PHASE 1: ***************************************/
    
    for (i = 0; i < total_nCandidates; i++)
       permute[i] = -1;
    
    memcpy (cmatch, match, hg->nVtx * sizeof(int));  /* for temporary locking */
    if (cFLAG)
      /* select upto nCandidates unmatched vertices to locally match */
      for (nTotal=0; nTotal < nCandidates && vindex < hg->nVtx; vindex++)  {
        if (cmatch[visit[vindex]] == visit[vindex])  {         /* unmatched */
          permute[nTotal++] = visit[vindex];    /* select it as a candidate */
          cmatch[visit[vindex]] = -1;           /* mark it as a pending match */
        }
      }
    else  {       
      /* Select upto nCandidates unmatched vertices to globally match. */
      for (sendcnt = 0; sendcnt < nCandidates && vindex < hg->nVtx; vindex++)  {
        if (cmatch[visit[vindex]] == visit[vindex])  {         /* unmatched */
          select[sendcnt++] = visit[vindex];    /* select it as a candidate */
          cmatch[visit[vindex]] = -1;           /* mark it as a pending match */
        }
      }
      nselect = sendcnt;                          /* save for later use */
                        
      /* assure send buffer is large enough by first computing required size */
      sendsize = (HEADER_COUNT-1) * sendcnt;      /* takes care of header */
      for (i = 0; i < sendcnt; i++)  {
        lno = select[i];
        sendsize += (hg->vindex[lno+1] - hg->vindex[lno]);
      }
      if (sendsize > nSend)
        MACRO_REALLOC (sendsize, nSend, send);     /* make send buffer bigger */
    
      /* fill send buff: list of 
       * <candidate_gno, candidate_gno's edge count, list of edge lno's> 
       */
      s = send;
      n = 0;
      candidate_index = first_candidate_index;
      for (i = 0; i < sendcnt; i++)   {
        lno = select[i];
        candidate_gno = VTX_LNO_TO_GNO(hg, lno);
        /* Optimization: Send only vertices that are non-empty locally */
        if ((hg->vindex[lno+1] > hg->vindex[lno]) 
         || (candidate_gno % hgc->nProc_y == hgc->myProc_y)) {
            n++;  /* non-empty vertex */
            *s++ = candidate_gno;
            *s++ = candidate_index;
            *s++ = hg->vindex[lno+1] - hg->vindex[lno];              /* count */
            for (j = hg->vindex[lno]; j < hg->vindex[lno+1]; j++)  
              *s++ = hg->vedge[j];                             /* lno of edge */
        }
        candidate_index++;
      }
      sendsize = s - send;
    
      /* determine actual global number of candidates this round */
      /* n is actual number of local non-empty vertices */
      /* nTotal is the global number of candidate vertices sent in this row */
      MPI_Allreduce (&n, &nTotal, 1, MPI_INT, MPI_SUM, hgc->row_comm);
      MPI_Allreduce (&nTotal, &max_nTotal, 1, MPI_INT, MPI_MAX, hgc->col_comm);
      if (max_nTotal == 0) {
        if (hgp->use_timers > 3)
           ZOLTAN_TIMER_STOP (zz->ZTime, timer[1], hg->comm->Communicator);
        break;                          /* globally all work is done, so quit */
      }
      
      /* communication to determine global size & displacements of rec buffer */
      MPI_Allgather (&sendsize, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm);
     
      /* determine size of the rec buffer & reallocate bigger iff necessary */
      recsize = 0;
      for (i = 0; i < hgc->nProc_x; i++)
        recsize += size[i];          /* compute total size of edgebuf in ints */
      if (recsize > nEdgebuf)
        MACRO_REALLOC (recsize, nEdgebuf, edgebuf);        /* enlarge edgebuf */
    
      /* setup displacement array necessary for MPI_Allgatherv */
      dest[0] = 0;
      for (i = 1; i < hgc->nProc_x; i++)
        dest[i] = dest[i-1] + size[i-1];

      /* communicate vertices & their edges to all row neighbors */
      MPI_Allgatherv(send, sendsize, MPI_INT, edgebuf, size, dest, MPI_INT,
       hgc->row_comm);
         
      /* create random permutation of index into the edge buffer */
      i = 0;
      for (j = 0 ; j < nTotal  &&  i < recsize; j++)   {
        permute[j] = i++;         /* save index of candidate_gno in permute[] */
        i++;                      /* skip candidate_index */
        count = edgebuf[i++];     /* count of edges */
        i += count;               /* skip over count edges */
      }

      /* Communication has grouped candidates by processor, rescramble!     */
      /* Otherwise all candidates from proc column 0 will be matched first, */
      /* Future: Instead of Zoltan_Rand_Perm_Int, we could use              */
      /* Zoltan_PHG_Vertex_Visit_Order() to reorder the candidates          */
      /* but that routine uses a local hg so won't work on the candidates.  */
      if (hgc->nProc_x > 1) {
        Zoltan_Srand_Sync(Zoltan_Rand(NULL),&(hgc->RNGState_col),hgc->col_comm);
#ifdef RTHRTH
        Zoltan_Rand_Perm_Int (permute, total_nCandidates, &(hgc->RNGState_col));
#endif 
      }
    }                           /* DONE:  if (cFLAG) else ...  */
    MACRO_TIMER_STOP (1);
    
    /************************ PHASE 2: ***************************************/
      
    /* for each candidate vertex, compute all local partial inner products */
    kstart = old_kstart = 0;         /* next candidate (of nTotal) to process */
    while (kstart < total_nCandidates)  { 

if (kstart > 0)        
printf ("RTHRTH    kstart > 0   RTHRTH\n");    

      MACRO_TIMER_START (2, "Matching kstart A", 0);
      sendsize = 0;                      /* position in send buffer */
      sendcnt = 0;                       /* count of messages in send buffer */
      s = send;                          /* start at send buffer origin */
      for (k = kstart; k < total_nCandidates; k++)  {  
        if (!cFLAG)  {
          if (permute[k] == -1)
             continue; 
          r     = &edgebuf[permute[k]];
          candidate_gno   = *r++;          /* gno of candidate vertex */
          candidate_index = *r++;          /* candidate_index of vertex */
          count = *r++;                    /* count of following hyperedges */
        }
        else
          candidate_gno = permute[k];  /* need to use next local vertex */
                                  /* here candidate_gno is really a local id */
                  
        /* now compute the row's nVtx inner products for kth candidate */
        m = 0;
        if (!cFLAG) {
          if      ((hg->ewgt == NULL) && (hgp->vtx_scal == NULL))
            INNER_PRODUCT1(1.0)
          else if ((hg->ewgt == NULL) && (hgp->vtx_scal != NULL))
            INNER_PRODUCT1(hgp->vtx_scal[hg->hvertex[j]])
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal == NULL))
            INNER_PRODUCT1(hg->ewgt[*r])
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal != NULL))
            INNER_PRODUCT1(hgp->vtx_scal[hg->hvertex[j]] * hg->ewgt[*r])
                  
        } else /* cFLAG */ {
          if      ((hg->ewgt == NULL) && (hgp->vtx_scal == NULL))
            INNER_PRODUCT2(1.0)
          else if ((hg->ewgt == NULL) && (hgp->vtx_scal != NULL))
            INNER_PRODUCT2(hgp->vtx_scal[hg->hvertex[j]])
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal == NULL))
            INNER_PRODUCT2(hg->ewgt[edge])
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal != NULL))
            INNER_PRODUCT2(hgp->vtx_scal[hg->hvertex[j]] * hg->ewgt[edge])   
        }
          
        /* if local vtx, remove self inner product (useless maximum) */
        if (cFLAG)
          sums[candidate_gno] = 0.0;   /* for cFLAG, candidate_gno is really 
                                          a local id */
        else if (VTX_TO_PROC_X(hg, candidate_gno) == hgc->myProc_x)
          sums[VTX_GNO_TO_LNO(hg, candidate_gno)] = 0.0;
         
        /* count partial sums exceeding PSUM_THRESHOLD */   
        count = 0;
        for (i = 0; i < m; i++)  {
          lno = index[i];
          if (sums[lno] > PSUM_THRESHOLD)
            aux[count++] = lno;      /* save lno for significant partial sum */
          else
            sums[lno] = 0.0;         /* clear unwanted entries */  
        }     
        if (count == 0)
          continue;         /* no partial sums to append to message */

        /* HEADER_COUNT (row, candidate_gno, candidate_index, count of
         * <lno, psum> pairs describing non-zero partial inner products)  */
        msgsize = HEADER_COUNT + 2 * count;
        
        /* iff necessary, resize send buffer to fit at least first message */
        if (sendcnt == 0 && (msgsize > nSend))  {
          MACRO_REALLOC (5*msgsize, nSend, send);  /* make send buffer bigger */
          s = send;    
        }

        if (sendsize + msgsize <= nSend)  {
          /* current partial sums fit, so put them into the send buffer */
          dest[sendcnt]   = candidate_gno % hgc->nProc_y; /* proc to compute 
                                                             total sum */
          size[sendcnt++] = msgsize;          /* size of message */
          sendsize       += msgsize;          /* cummulative size of message */
          
          *s++ = hgc->myProc_y;               /* save my row (for merging) */
          *s++ = candidate_gno;
          *s++ = candidate_index;
          *s++ = count;
          for (i = 0; i < count; i++)  {          
            *s++ = aux[i];                          /* lno of partial sum */
             f = (float*) s++;
            *f = sums[aux[i]];                      /* partial sum */           
            sums[aux[i]] = 0.0;
          }          
        }
        else  {           /* psum message doesn't fit into buffer */
          for (i = 0; i < count; i++)              
            sums[aux[i]] = 0.0;        
          break;
        }  
      }                  /* DONE: loop over k */                    

      MACRO_TIMER_STOP (2);
      MACRO_TIMER_START (3, "Matching kstart B", 0);
     
      /* synchronize all rows in this column to next kstart value */
      old_kstart = kstart;      

      MPI_Allreduce (&k, &kstart, 1, MPI_INT, MPI_MIN, hgc->col_comm);

            
      /* Send inner product data in send buffer to appropriate rows */
      err = communication_by_plan (zz, sendcnt, dest, size, 1, send, &reccnt, 
       &recsize, &nRec, &rec, hgc->col_comm, IPM_TAG);
      if (err != ZOLTAN_OK)
        goto fini;
      
      /* build index into receive buffer pointer for each new row of data */
      old_row = -1;
      k = 0;
      for (r = rec; r < rec + recsize  &&  k < hgc->nProc_y; )  {     
        row = *r++;        
        if (row != old_row)  {
          index[k++] = r - rec;   /* points at candidate_gno, not row */
          old_row = row;
        }
        candidate_gno   = *r++;
        candidate_index = *r++;
        count           = *r++;
        r += (2 * count);
      }
     
      /* save current positions into source rows within rec buffer */
      for (i = 0; i < k; i++)
        rows[i] = &rec[index[i]];
      for (i = k; i <= hgc->nProc_y; i++)
        rows[i] = &rec[recsize];       /* in case no data came from a row */
      
      /* merge partial i.p. sum data to compute total inner products */
      s = send; 
      for (n = old_kstart; n < kstart; n++)  {
        m = 0;        
        /* here candidate_gno is really a local id when cFLAG */
        candidate_gno = (cFLAG) ? permute[n] : edgebuf[permute[n]];
               
        /* Not sure if this test makes any speedup ???, works without! */
        if (candidate_gno % hgc->nProc_y != hgc->myProc_y)
          continue;                      /* this candidate_gno's partial IPs 
                                            not sent to this proc */
        
        /* merge step: look for target candidate_gno from each row's data */
        for (i = 0; i < hgc->nProc_y; i++)  {
          if (rows[i] < &rec[recsize] && *rows[i] == candidate_gno)  {       
            candidate_index = *(++rows[i]);   /* skip candidate index */
            count = *(++rows[i]);
            for (j = 0; j < count; j++)  {
              lno = *(++rows[i]);         
              if (sums[lno] == 0.0)       /* is this first time for this lno? */
                aux[m++] = lno;           /* then save the lno */          
              sums[lno] += *(float*) (++rows[i]);    /* sum the psums */
            }
            rows[i] += 2;                 /* skip past current psum,row */
          }
        }
          
        /* determine how many total inner products exceed threshold */
        count = 0;
        for (i = 0; i < m; i++)
          if (sums[aux[i]] > TSUM_THRESHOLD)
            count++;   

        /* create <candidate_gno, candidate_index (if CANDIDATE_INDICES),
         * count of <lno,tsum> pairs, <lno, tsum>> in send array.         */
        if (count > 0)  {
          if ( (s - send) + ((HEADER_COUNT-1) + 2 * count) > nSend ) 
          {
            sendsize = s - send;
            MACRO_REALLOC (nSend + (HEADER_COUNT-1)+2*count, nSend, send); 
                                   /* enlarge buffer */
            s = send + sendsize;   /* since realloc buffer could move */ 
          }      
          *s++ = candidate_gno;
          *s++ = candidate_index;
          *s++ = count;
        }  
        for (i = 0; i < m; i++)   {
          lno = aux[i];             
          if (sums[lno] > TSUM_THRESHOLD)  {
            *s++ = lno;
             f = (float*) s++;
            *f = sums[lno];
          }  
          sums[lno] = 0.0;  
        }     
      }
      sendsize = s - send;   /* size (in ints) of send buffer */
      
      /* Communicate total inner product results to MASTER ROW */
      MPI_Gather(&sendsize, 1, MPI_INT, size, 1, MPI_INT, 0, hgc->col_comm);

      if (hgc->myProc_y == 0) {
        recsize = 0;
        for (i = 0; i < hgc->nProc_y; i++)
          recsize += size[i];        
          
        dest[0] = 0;
        for (i = 1; i < hgc->nProc_y; i++)
          dest[i] = dest[i-1] + size[i-1];
        
        if (recsize > nRec)
          MACRO_REALLOC (recsize, nRec, rec);      /* make rec buffer bigger */
      }

      MPI_Gatherv(send, sendsize, MPI_INT, rec, size, dest, MPI_INT, 0,
                  hgc->col_comm);
       
      /* Determine best vertex and best sum for each candidate */
      if (hgc->myProc_y == 0) {   /* do following only if I am the MASTER ROW */
        for (r = rec; r < rec + recsize;)  {
          candidate_gno   = *r++;                    /* candidate's GNO */
          candidate_index = *r++;                    /* candidate's index */

          count = *r++;                          /* count of nonzero pairs */
          bestsum = -1.0;                        /* any negative value will do */
          bestlno = -1;                          /* any negative value will do */
          for (i = 0; i < count; i++)  {
            lno =          *r++;
            f   =  (float*) r++;     
            if (*f > bestsum  &&  cmatch[lno] == lno)  {
              bestsum = *f;
              bestlno = lno;
            }      
          }
         
          if (cFLAG && bestsum > TSUM_THRESHOLD)  {
            match[bestlno]       = candidate_gno;
            match[candidate_gno] = bestlno;
            cmatch[bestlno] = -1;         
          }
                        
          if (!cFLAG && bestsum > TSUM_THRESHOLD)  {
            cmatch[bestlno] = -1;  /* mark pending match to avoid conflicts */
            master_data[candidate_index].candidate = candidate_gno;
            master_data[candidate_index].partner = VTX_LNO_TO_GNO (hg, bestlno);
            master_data[candidate_index].ip = bestsum;
          }
        }
      }
        
      if (cFLAG) {  /* Broadcast what we matched so far */
        MPI_Bcast (match, hg->nVtx, MPI_INT, 0, hgc->col_comm); 
      }
      MACRO_TIMER_STOP (3);    

    }            /* DONE: kstart < max_nTotal loop */ 

    if (cFLAG)
      continue;      /* skip phases 3 and 4, continue rounds */ 
    

    /************************ NEW PHASE 3: ********************************/

    MACRO_TIMER_START (4, "Matching Phase 3", 1);

    /* Only MASTER ROW computes best global match for candidates */
    /* EBEB or perhaps we can do this fully distributed? */
    if (hgc->myProc_y == 0) {
      HG_Ptr = hg;
      MPI_Allreduce(master_data, global_best, total_nCandidates, phasethreetype,
                    phasethreeop, hgc->row_comm);

      /* Look through array of "winners" and update match array. */
      /* Local numbers are used for local matches, otherwise
         -(gno+1) is used in the match array.                    */
      for (i = 0; i < total_nCandidates; i++) {
        int cproc, vproc;
        candidate_gno = global_best[i].candidate;

        /* Reinitialize master_data for next round */
        master_data[i].candidate = -1;
        master_data[i].partner   = -1;
        master_data[i].ip        = -1.0;
        if (candidate_gno == -1)
          continue;

        partner_gno = global_best[i].partner;
        cproc = VTX_TO_PROC_X(hg, candidate_gno);
        vproc = VTX_TO_PROC_X(hg, partner_gno);
        if (cproc == hgc->myProc_x) {
          if (vproc == hgc->myProc_x)   {
            int v1 = VTX_GNO_TO_LNO(hg, partner_gno);
            int v2 = VTX_GNO_TO_LNO(hg, candidate_gno);
            match[v1] = v2;
            match[v2] = v1;
          }
          else 
            match[VTX_GNO_TO_LNO(hg, candidate_gno)] = -partner_gno - 1;
        }                         
        else if (vproc == hgc->myProc_x)
          match[VTX_GNO_TO_LNO(hg, partner_gno)] = -candidate_gno - 1;
      }
    } /* End (hgc->myProc_y == 0) */

    /* broadcast match array to the entire column */
    MPI_Bcast (match, hg->nVtx, MPI_INT, 0, hgc->col_comm);
    MACRO_TIMER_STOP (4);                       /* end of phase 3 */
  }                                             /* DONE: loop over rounds */
  
  MACRO_TIMER_START (6, "Matching Cleanup", 0);

  /* optional sanity tests */
  if (zz->Debug_Level > 4 && hgc->myProc_x == 0 && hgc->myProc_y == 0)  {
    int local = 0, global = 0, unmatched = 0;
    for (i = 0; i < hg->nVtx; i++)  {
      if      (match[i] == i)  unmatched++;
      else if (match[i] < 0)   global++;
      else                     local++;
    }
    uprintf (hgc, "%d RTHRTH %d unmatched, %d external, %d local of %d\n",
     hg->info, unmatched, global, local, hg->nVtx);
  }

  if (zz->Debug_Level > 4 && hgc->myProc_x==0 && hgc->myProc_y==0)
    fprintf (stdout, "RTHRTH rounds %d\n", nRounds);

  if (zz->Debug_Level > 4)  {
    /* The following tests that the global match array is a valid permutation */
    /* NOTE:  THESE TESTS ARE NOT MANDATORY; THEY CAN BE EXCLUDED AFTER WE    */
    /* COMPLETE TESTING OF matching_ipm.                                      */

    for (i = 0; i < hg->nVtx; i++)
      if (match[i] < 0)  cmatch[i] = -match[i] - 1;
      else               cmatch[i] = VTX_LNO_TO_GNO(hg, match[i]);

    MPI_Allgather (&hg->nVtx, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm); 

    recsize = 0;
    for (i = 0; i < hgc->nProc_x; i++)
      recsize += size[i];
  
    dest[0] = 0;
    for (i = 1; i < hgc->nProc_x; i++)
      dest[i] = dest[i-1] + size[i-1];
  
    if (nRec < recsize)
      MACRO_REALLOC (recsize, nRec, rec); /* make rec buffer bigger */
    MPI_Allgatherv (cmatch, hg->nVtx, MPI_INT, rec, size, dest, MPI_INT,
     hgc->row_comm);

    if (nSend < recsize)
      MACRO_REALLOC (recsize, nSend, send);  /* make send buffer bigger */
  
    for (i = 0; i < recsize; i++)
      send[i] = 0;
    for (i = 0; i < recsize; i++)
      ++send[rec[i]];

    count = 0;
    for (i = 0; i < recsize; i++)
      if (send[i] != 1)
        count++;
    if (count)    
      uprintf (hgc, "RTHRTH %d FINAL MATCH ERRORS of %d\n", count, recsize); 
  }
  MACRO_TIMER_STOP (6);

fini:
  if (!cFLAG) {
    MPI_Op_free(&phasethreeop);
    MPI_Type_free(&phasethreetype);
    ZOLTAN_FREE(&global_best);
    ZOLTAN_FREE(&Tmp_Best);
  }

  Zoltan_Multifree (__FILE__, __LINE__, 15, &cmatch, &visit, &sums, &send,
   &dest, &size, &rec, &index, &aux, &permute, &edgebuf, &select, &rows,
   &master_data, &master_procs);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}


  
/****************************************************************************/
static int communication_by_plan (ZZ* zz, int sendcnt, int* dest, int* size, 
 int scale, int* send, int* reccnt, int *recsize, int* nRec, int** rec,
 MPI_Comm comm, int tag)
{
   ZOLTAN_COMM_OBJ *plan = NULL;
   int err;
   char *yo = "communication_by_plan";
   
   /* communicate send buffer messages to other row/columns in my comm */  
   err = Zoltan_Comm_Create (&plan, sendcnt, dest, comm, tag, reccnt);
   if (err != ZOLTAN_OK) {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed to create plan");
     return err;
   }
        
   /* resize plan if necessary */
   if (size != NULL) {
     err = Zoltan_Comm_Resize (plan, size, tag+1, recsize);
     if (err != ZOLTAN_OK) {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed to resize plan");
       return err;
     }
     scale = 1;       
   }
   else
     *recsize = *reccnt * scale;
   
   /* realloc rec buffer if necessary */  
   if (*recsize > *nRec)  {   
     *nRec = *recsize;
     if (!(*rec = (int*) ZOLTAN_REALLOC (*rec, *nRec * sizeof(int))))  {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error");
       return ZOLTAN_MEMERR;
     }
   }
   
   /* send messages from send buffer to destinations */      
   err = Zoltan_Comm_Do (plan, tag+2, (char*) send, scale * sizeof(int),
    (char*) *rec);
   if (err != ZOLTAN_OK)  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed in Comm_Do");
     return err;
   }
   
   /* free memory associated with the plan */
   Zoltan_Comm_Destroy (&plan); 
   return ZOLTAN_OK;
}

/***************************************************************************/

/* MPI_Op for computing the maximum inner-product for each candidate */
static void phasethreereduce (
  void *tin, 
  void *tinout, 
  int  *tnum, 
  MPI_Datatype *mytype)
{
  int num = *tnum;
  Triplet *in    = (Triplet*) tin;
  Triplet *inout = (Triplet*) tinout;
  int i;

  for (i = 0; i < num; i++) {
    if (in[i].candidate == -1 && inout[i].candidate == -1) 
      continue;                         /* No values set for this candidate */

    if (in[i].ip > inout[i].ip)
      inout[i] = in[i];                 /* in has larger inner product */
    else if (in[i].ip == inout[i].ip) {
      int in_proc    = VTX_TO_PROC_X (HG_Ptr, in[i].partner);
      int inout_proc = VTX_TO_PROC_X (HG_Ptr, inout[i].partner);
      int cand_proc  = VTX_TO_PROC_X (HG_Ptr, in[i].candidate);

      /* Give preference to partners on candidate's processor */
      if (((in_proc == cand_proc) && (inout_proc == cand_proc))
       || ((in_proc != cand_proc) && (inout_proc != cand_proc))) {
        /* Both partners are on candidate's proc OR
           neither partner is on candidate's proc.
           Break ties by larger partner gno. */
        if (in[i].partner > inout[i].partner)
          inout[i] = in[i];
      }
      else if (in_proc == cand_proc) {
        inout[i] = in[i];   /* Give preference to local partner */
      }
    } 
  }
}


#undef MACRO_REALLOC
#undef MACRO_TIMER_START
#undef MACRO_TIMER_STOP
#undef INNER_PRODUCT
#undef INNER_PRODUCT2
#undef ROUNDS_CONSTANT
#undef IPM_TAG
#undef HEADER_COUNT
#undef PSUM_THRESHOLD
#undef TSUM_THRESHOLD

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif


