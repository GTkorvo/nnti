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

#ifndef __ZOLTAN_MATRIX_PARTITION
#define __ZOLTAN_MATRIX_PARTITION

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*
 * Assumption is we are creating a hypergraph and using Zoltan_PHG
 * to partition objects, but eventually this code should come out
 * of the phg directory and be more general, using other partitioning
 * methods on the sparse matrix.
 */
#include "phg.h"
/*******************/
#include "zz_util_const.h"

/* Query functions set up by Zoltan_Matrix_Partition for the
 * Zoltan library
 */

/* Parameters for Zoltan method SPARSE_MATRIX
 */
static PARAM_VARS MP_params[] = {
  {"LB_APPROACH",         NULL,  "STRING",     0},
  {"NUM_GID_ENTRIES",     NULL,  "INT",     0},
  {NULL,                  NULL,  NULL,     0}
};
    
/***********************************************************************
** Allow user of parallel hypergraph methods to input a sparse matrix, 
** and have Zoltan create the hypergraph that best solves their problem.
**
** See Erik Boman's message to Isorropia-developers at
** http://software.sandia.gov/pipermail/isorropia-developers/2006-August/000065.html
**
** This avoids the necessity of having the user figure out how our
** hypergraph partitioning codes work.
**
** This source file uses the terms "pins" and "non-zeros" interchangeably.
**
** An IJTYPE is the type of a row or column ID.
**
** LB_METHOD = SPARSE_MATRIX
** LB_APPROACH:
**    PHG_ROWS    use Zoltan's PHG to partition the rows (objects are
**                rows, hyperedges are columns)
**                
**    PHG_COLS    use Zoltan's PHG to partition the columns (objects
**                are columns, hyperedges are rows)
**
** If you want to add a new approach to Zoltan_Matrix_Partition, you
** need to write functions analygous to the 5 functions written
** for the PHG_ROWS and PHG_COLS approaches: phg_rc_setup, phg_rc_results,
** phg_rc_get_pins, phg_rc_get_rows, and phg_rc_get_columns.
**
*************************************************************************/

/*
 * General matrix partition functions
 */
static int process_matrix_input(ZZ *zz, ZOLTAN_MP_DATA *mpd);
static int make_mirror(ZOLTAN_MP_DATA *mpd);
static int allocate_copy(IJTYPE **to, IJTYPE *from, IJTYPE len);

static ZOLTAN_MP_DATA *MP_Initialize_Structure();
static int MP_Initialize_Params(ZZ *zz, ZOLTAN_MP_DATA *mpd);

/*
 * Functions specifically to process approaches PHG_ROWS and PHG_COLUMNS
 */
static int phg_rc_setup(ZOLTAN_MP_DATA *mpd);
static int phg_rc_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  int *export_procs, int *export_to_part);
static int phg_rc_get_pins(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nPins, IJTYPE *I, IJTYPE *J, int *ijProcs, int *ijParts);
static int phg_rc_get_rows(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nRows, IJTYPE *rowIDs, int *rowProcs, int *rowParts);
static int phg_rc_get_columns(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nCols, IJTYPE *colIDs, int *colProcs, int *colParts);

/*
 * Functions to create a search structure to locate rows, columns
 * or non-zeros in the ZOLTAN_MP_DATA structure.  To support
 * sparse matrix queries that can be called when Zoltan_Matrix_Partition
 * has completed.  The lookup function Zoltan_Lookup_Obj is global.
 */
static void free_obj_lookup_table(obj_lookup **lu);
static obj_lookup *create_obj_lookup_table(IJTYPE numObjs, IJTYPE *objGIDs);
static obj_lookup *create_obj_lookup_table2(IJTYPE sizeI, IJTYPE sizeJ,
                   IJTYPE *listI, IJTYPE *indexJ, IJTYPE *listJ);

/****************************************************************/
/****************************************************************/
/* API:
 *
 *  Zoltan_Matrix_Partition()     must be called by all processes
 *
 *  Queries (not all queries are defined for all LB_APPROACH):
 *
 *  Zoltan_MP_Get_Row_Assignment()
 *  Zoltan_MP_Get_Column_Assignment()
 *  Zoltan_MP_Get_NonZero_Assignment()
 *
 *  To evaluate partitioning:
 *
 *  Zoltan_Matrix_Partition_Eval()   must be called by all procs
 */

int Zoltan_MP_Get_Row_Assignment(ZZ *zz, int nRows, IJTYPE *rowIDs,
        int *rowProcs, int *rowParts)
{
  char *yo = "Zoltan_MP_Get_Row_Assignment";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  int ierr = ZOLTAN_FATAL;

  /*
   * Row assignment depends on the LB_APPROACH.  We can only
   * return row assignments for rows that this process gave us
   * non-zeroes for.
   */

  if (mpd){
    switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        ierr = phg_rc_get_rows(mpd, nRows, rowIDs, rowProcs, rowParts);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
    }
  }

  return ierr;
}
int Zoltan_MP_Get_Column_Assignment(ZZ *zz, int nCols, IJTYPE *colIDs,
        int *colProcs, int *colParts)
{
  char *yo = "Zoltan_MP_Get_Column_Assignment";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  int ierr = ZOLTAN_FATAL;

  /*
   * Column assignment depends on the LB_APPROACH.  We can only
   * return column assignments for columns that this process gave us
   * non-zeroes for.
   */

  if (mpd){
    switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        ierr = phg_rc_get_columns(mpd, nCols, colIDs, colProcs, colParts);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
    }
  }

  return ierr;
}
int Zoltan_MP_Get_NonZero_Assignment(ZZ *zz, int nNZ, 
        IJTYPE *rowIDs, IJTYPE *colIDs, int *nzProcs, int *nzParts)
{
  char *yo = "Zoltan_MP_Get_NonZero_Assignment";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  int ierr = ZOLTAN_FATAL;

  /*
   * Nonzero assignment depends on the LB_APPROACH.  We can only
   * return nonzero assignments for nonzeroes that were returned
   * by this process in the query functions.
   */

  if (mpd){
    switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        ierr = phg_rc_get_pins(mpd, nNZ, rowIDs, colIDs, nzProcs, nzParts);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
    }
  }

  return ierr;
}
int Zoltan_Matrix_Partition(ZZ *zz)
{
  char *yo = "Zoltan_Matrix_Partition";
  int ierr = ZOLTAN_OK;
  ZZ *zzLib = NULL;
  IJTYPE npins;

  /* The persistent structure we will save at zz->LB.Data_Structure */

  ZOLTAN_MP_DATA *mpd = MP_Initialize_Structure();

  if (mpd == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (zz->LB.Data_Structure){
    
    zz->LB.Free_Structure(zz);
  }
  zz->LB.Data_Structure = mpd;

  /* Process parameters, write to mpd.  */

  ierr = MP_Initialize_Params(zz, mpd);

  if (ierr != ZOLTAN_OK){
    goto End;
  }

  /*
   * Create the Zoltan problem that we will define from the input sparse matrix.
   */

  zzLib = Zoltan_Create(zz->Communicator);
  mpd->zzLib = zzLib;

  /* Call application defined query functions to get the non-zeros.
   * User has a choice of giving us compressed rows or compressed columns.
   * We need global row and column IDs.  Any process can give us any of the
   * non-zeros, but each non-zero must be supplied by only one process.
   */

  ierr = ZOLTAN_FATAL;

  if ((zz->Get_CSC_Size != NULL) && (zz->Get_CSC != NULL)){
    ierr = ZOLTAN_OK;
    mpd->input_type = COL_TYPE;
  }
  else if ((zz->Get_CSR_Size != NULL) && (zz->Get_CSR != NULL)){
    mpd->input_type = ROW_TYPE;
    ierr = ZOLTAN_OK;
  }

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Required query functions are not defined.\n");
    goto End;
  }

  if (mpd->input_type == COL_TYPE){
    zz->Get_CSC_Size(zz->Get_CSC_Size_Data, &mpd->numRC, &npins, &ierr);
  }
  else {
    zz->Get_CSR_Size(zz->Get_CSR_Size_Data, &mpd->numRC, &npins, &ierr);
  }

  mpd->rcGID = (IJTYPE *)ZOLTAN_MALLOC(mpd->numRC * sizeof(IJTYPE));
  mpd->pinIndex = (IJTYPE *)ZOLTAN_MALLOC((mpd->numRC+1) * sizeof(IJTYPE));
  mpd->pinGID = (IJTYPE *)ZOLTAN_MALLOC(npins * sizeof(IJTYPE));
  if ((mpd->numRC && !mpd->rcGID) || !mpd->pinIndex || (npins && !mpd->pinGID)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (mpd->input_type == COL_TYPE){
    zz->Get_CSC(zz->Get_CSC_Data, mpd->numRC, npins, 
                mpd->rcGID, mpd->pinIndex, mpd->pinGID, &ierr);
  }
  else {
    zz->Get_CSR(zz->Get_CSR_Data, mpd->numRC, npins, 
                mpd->rcGID, mpd->pinIndex, mpd->pinGID, &ierr);
  }
 
  mpd->pinIndex[mpd->numRC] = npins;

  /*
   * Determine globally how many rows and columns were returned, and if
   * their IDs are 0 based or 1 based.  (Assume they are contiguous.)
   * Some row/column numbers may not appear in input if they have no pins.
   */

  ierr = process_matrix_input(zz, mpd);

  /*
   * Set the query functions and parameters required for this problem.
   * Precompute hypergraph if necessary.  This is the code that needs
   * to be written if you are implementing a new LB_APPROACH.
   */
  switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        phg_rc_setup(mpd);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
        ierr = ZOLTAN_FATAL;
        goto End;
    }

  /*
   * Call Zoltan_LB_Partition to partition the objects
   */

  int changes, num_gid_entries, num_lid_entries, num_import, num_export;
  unsigned int *import_global_ids, *import_local_ids;
  unsigned int *export_global_ids, *export_local_ids;
  int *import_to_part, *export_to_part;
  int *import_procs, *export_procs;

  ierr = Zoltan_LB_Partition(zzLib,
      &changes, &num_gid_entries, &num_lid_entries, 
      &num_import,
      &import_global_ids, &import_local_ids, &import_procs, &import_to_part,
      &num_export,
      &export_global_ids, &export_local_ids, &export_procs, &export_to_part);

  /*
   * Save the data required to respond to the queries that are
   * supported by each LB_APPROACH.
   */

  switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        /*
         * Save row or column assignment for all rows or columns
         * in my non-zeroes.  Each of my non-zeroes is assigned to
         * the row or column that it's row or column is assigned to.
         *
         * PHG_ROWS - we only give row and nonzero assignments
         * PHG_COLUMNS - we only give column and nonzero assignments
         */
        phg_rc_result(mpd, num_export,
            export_global_ids, export_local_ids, export_procs, export_to_part);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
        ierr = ZOLTAN_FATAL;
        goto End;
    }

End:
  Zoltan_Multifree(__FILE__, __LINE__, 8,
    &import_global_ids, &import_local_ids, &import_procs, &import_to_part,
    &export_global_ids, &export_local_ids, &export_procs, &export_to_part);

  return ierr;
}
/****************************************************************/
/****************************************************************/
/*
 * Called by Zoltan_Destroy()
 */
void Zoltan_MP_Free_Structure(ZZ *zz)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;

  if (mpd != NULL){

    ZOLTAN_FREE(&mpd->rcGID);
    ZOLTAN_FREE(&mpd->pinIndex);
    ZOLTAN_FREE(&mpd->pinGID);

    ZOLTAN_FREE(&mpd->crGID);
    ZOLTAN_FREE(&mpd->mirrorPinIndex);
    ZOLTAN_FREE(&mpd->mirrorPinGID);

    ZOLTAN_FREE(&mpd->vtxGID);
    ZOLTAN_FREE(&mpd->vtxWgt);
    ZOLTAN_FREE(&mpd->hindex);
    ZOLTAN_FREE(&mpd->hvertex);

    ZOLTAN_FREE(&mpd->rowproc);
    ZOLTAN_FREE(&mpd->rowpart);
    ZOLTAN_FREE(&mpd->colproc);
    ZOLTAN_FREE(&mpd->colpart);
    ZOLTAN_FREE(&mpd->pinproc);
    ZOLTAN_FREE(&mpd->pinpart);

    free_obj_lookup_table(&mpd->row_lookup);
    free_obj_lookup_table(&mpd->col_lookup);
    free_obj_lookup_table(&mpd->pin_lookup);

    Zoltan_Destroy(&(mpd->zzLib));  /* we created this PHG problem */

    ZOLTAN_FREE(&zz->LB.Data_Structure);
  }
}

/****************************************************************/
/****************************************************************/
/*
 * Functions used by Zoltan_Matrix_Partition()
 */

static int MP_Initialize_Params(ZZ *zz, ZOLTAN_MP_DATA *mpd)
{
  char *yo="MP_Initialize_Params";
  int ierr = ZOLTAN_OK;
  char approach[MAX_PARAM_STRING_LEN];

  Zoltan_Bind_Param(MP_params, "LB_APPROACH", approach);  
  Zoltan_Bind_Param(MP_params, "NUM_GID_ENTRIES", &mpd->gidLen);  

  strncpy(approach, "PHG_ROWS", MAX_PARAM_STRING_LEN);
  mpd->gidLen = 1;

  ierr = Zoltan_Assign_Param_Vals(zz->Params, MP_params, zz->Debug_Level,
          zz->Proc, zz->Debug_Proc);

  if (ierr == ZOLTAN_OK){
    /* Update this if you add a new LB_APPROACH */
    if (!strcasecmp(approach, "phg_rows"))
      mpd->approach = PHG_ROWS;
    else if (!strcasecmp(approach, "phg_columns"))
      mpd->approach = PHG_COLUMNS;

    if (mpd->gidLen > sizeof(IJTYPE)){
      /*
       * Maybe they want to use long ints, we only handle ints.  Change
       * IJTYPE to make code use long ints.
       */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input row and column IDs are too large");
      ierr = ZOLTAN_FATAL;
    }
  }
  return ierr;
}


static int process_matrix_input(ZZ *zz, ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int i;
  long int lpins, gpins;
  IJTYPE minID, maxID, minPinID, maxPinID;
  IJTYPE npins=0;
  long int vals[2], gvals[2];

  /* get global number of rows, columns and pins, range of row and column IDs */

  if (mpd->numRC > 0)
    npins = mpd->pinIndex[mpd->numRC];

  lpins = (long int)npins;

  MPI_Allreduce(&lpins, &gpins, 1, MPI_LONG, MPI_SUM, zz->Communicator);

  mpd->nNonZeros = (IJTYPE)gpins;

  maxID = maxPinID = -1;

  if (mpd->numRC > 0){
    minID = maxID = mpd->rcGID[0];

    for (i=1; i<mpd->numRC; i++){
      if (mpd->rcGID[i] < minID) minID = mpd->rcGID[i];
      else if (mpd->rcGID[i] > maxID) maxID = mpd->rcGID[i];
    }
    if (npins > 0){
      minPinID = maxPinID = mpd->pinGID[0];
      for (i=1; i<npins; i++){
        if (mpd->pinGID[i] < minPinID) minPinID = mpd->pinGID[i];
        else if (mpd->pinGID[i] > maxPinID) maxPinID = mpd->pinGID[i];
      }
    }
  }
  vals[0] = (long int)maxID;
  vals[1] = (long int)maxPinID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MAX, zz->Communicator);

  maxID = (IJTYPE)gvals[0];
  maxPinID = (IJTYPE)gvals[1];

  if (npins == 0){
    minPinID = maxPinID;
    if (mpd->numRC == 0){
      minID = maxID;
    }
  }
  vals[0] = (long int)minID;
  vals[1] = (long int)minPinID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MIN, zz->Communicator);

  minID = (IJTYPE)gvals[0];
  minPinID = (IJTYPE)gvals[1];

  if (mpd->input_type == ROW_TYPE){
    mpd->rowBaseID = minID;
    mpd->colBaseID = minPinID;
    mpd->nRows = maxID - minID + 1;
    mpd->nCols = maxPinID - minPinID + 1;
  }
  else{
    mpd->rowBaseID = minPinID;
    mpd->colBaseID = minID;
    mpd->nRows = maxPinID - minPinID + 1;
    mpd->nCols = maxID - minID + 1;
  }

  return ierr;
}
/* TODO - need to copy structure too.
 */
static ZOLTAN_MP_DATA *MP_Initialize_Structure()
{
  ZOLTAN_MP_DATA *mpd = 
    (ZOLTAN_MP_DATA *)ZOLTAN_CALLOC(sizeof(ZOLTAN_MP_DATA), 1);

  if (mpd == NULL){
    return NULL;
  }

  return mpd;
}
static int make_mirror(ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int npins, nIDs;

  /*
   * If sparse matrix was supplied in CSC, create another in CSR format,
   * and vice versa.
   */
  ZOLTAN_FREE(&mpd->crGID);
  ZOLTAN_FREE(&mpd->mirrorPinIndex);
  ZOLTAN_FREE(&mpd->mirrorPinGID);

  nIDs = mpd->numRC;
  npins = mpd->pinIndex[nIDs];

  mpd->numCR = mpd->numRC;
  allocate_copy(&mpd->crGID, mpd->pinGID, npins);
  allocate_copy(&mpd->mirrorPinIndex, mpd->pinIndex, nIDs+1);
  allocate_copy(&mpd->mirrorPinGID, mpd->rcGID, nIDs);

  /* TODO  Convert_to_CSR thinks some of these fields are ints
   *         when IJTYPEs might not be ints.  FIX, but OK for now.
   *
   */

  ierr = Zoltan_Convert_To_CSR(mpd->zzLib, npins, (int *)mpd->pinIndex,
             /* The following get overwritten with mirror */
             (int *)&(mpd->numCR),
             (ZOLTAN_ID_PTR *)&(mpd->mirrorPinGID),
             (int **)&(mpd->mirrorPinIndex),
             (ZOLTAN_ID_PTR *)&(mpd->crGID));

  return ierr;
}
static int allocate_copy(IJTYPE **to, IJTYPE *from, IJTYPE len)
{
  *to = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * len);
  memcpy(*to, from, sizeof(IJTYPE) * len);
  return ZOLTAN_OK;
}
/****************************************************************/
/****************************************************************/
/* Functions to support a hash table mapping global IDs to local
 * IDs, which are indices into an array.  The global ID can be
 * a single IJTYPE or an IJTYPE pair.
 */
static void free_obj_lookup_table(obj_lookup **lu)
{
  obj_lookup *l = *lu;

  if (l == NULL) return;

  ZOLTAN_FREE(&l->htTop);
  ZOLTAN_FREE(&l->ht);
  ZOLTAN_FREE(lu);
}
/*
 * Look up array index for global ID or for (i,j) pair.
 * This is a global function, in case other parts of Zoltan
 * library want to use the lookup tables.
 */
int Zoltan_Lookup_Obj(obj_lookup *lu, IJTYPE I, IJTYPE J)
{
  struct obj_node *hn;
  unsigned int hashVal;
  ZOLTAN_ID_TYPE ids[2];

  ids[0] = I;
  if (lu->key_size > 1) ids[1] = J;

  hashVal = Zoltan_Hash(ids, lu->key_size, lu->table_size);

  for (hn=lu->ht[hashVal]; hn != NULL; hn = hn->next){
    if (hn->i == I){
      if ((lu->key_size == 1) || (hn->j == J)){
        return hn->objLID;
      }
    }
  }
  return -1;
}
/*
** Create a hash table mapping the list of GIDs in objGIDs to their
** array index in the objGIDs list.  GIDs in objGIDs are unique.
*/
static obj_lookup *create_obj_lookup_table(IJTYPE numObjs, IJTYPE *objGIDs)
{
  IJTYPE i;
  unsigned int hashVal;
  struct obj_node *hn;
  obj_lookup *lu = NULL;
  int size_id;

  lu = (obj_lookup *)ZOLTAN_MALLOC(sizeof(obj_lookup));
  if (!lu){
    return NULL;
  }

  lu->key_size = 1;
  lu->table_size = numObjs / 4;
  if (lu->table_size < 4) lu->table_size = 4;

  lu->ht = 
  (struct obj_node **)ZOLTAN_CALLOC(sizeof(struct obj_node*) , lu->table_size);

  hn = lu->htTop = 
  (struct obj_node *)ZOLTAN_MALLOC(sizeof(struct obj_node) * numObjs);

  if (lu->table_size && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  size_id = sizeof(IJTYPE) / sizeof(ZOLTAN_ID_TYPE);
  if (size_id < 1) size_id = 1;

  for (i=0; i<numObjs; i++){

    hn->i = objGIDs[i];
    hn->j = 0;
    hn->objLID = i;

    hashVal = Zoltan_Hash((ZOLTAN_ID_PTR)&objGIDs[i], size_id, lu->table_size);

    hn->next = lu->ht[hashVal];
    lu->ht[hashVal] = hn;

    hn++;
  }

  return lu;
}
/*
** Create a hash table mapping the (I,J) pairs to their location in
** listJ.  IDs in listI are unique.  Number of pairs is sizeJ.
*/
static obj_lookup *create_obj_lookup_table2(IJTYPE sizeI, IJTYPE sizeJ,
                   IJTYPE *listI, IJTYPE *indexJ, IJTYPE *listJ)
{
  IJTYPE i, j;
  unsigned int hashVal;
  struct obj_node *hn;
  obj_lookup *lu = NULL;
  ZOLTAN_ID_TYPE ids[2];

  if (sizeof(IJTYPE) > sizeof(ZOLTAN_ID_TYPE)){
    ZOLTAN_PRINT_ERROR(0, "create_obj_lookup_table2", 
                       "code rewrite required.\n");
    exit(1);
  }

  lu = (obj_lookup *)ZOLTAN_MALLOC(sizeof(obj_lookup));
  if (!lu){
    return NULL;
  }

  lu->key_size = 2;
  lu->table_size = sizeJ / 4;
  if (lu->table_size < 4) lu->table_size = 4;

  lu->ht = 
  (struct obj_node **)ZOLTAN_CALLOC(sizeof(struct obj_node*) , lu->table_size);

  hn = lu->htTop = 
  (struct obj_node *)ZOLTAN_MALLOC(sizeof(struct obj_node) * sizeJ);

  if (lu->table_size && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  for (i=0; i<sizeI; i++){
    for (j = indexJ[i]; j < indexJ[i+1]; j++){
      hn->i = listI[i];
      hn->j = listJ[j];
      hn->objLID = j;

      ids[0] = listI[i];
      ids[1] = listJ[j];

      hashVal = Zoltan_Hash(ids, 2, lu->table_size);

      hn->next = lu->ht[hashVal];
      lu->ht[hashVal] = hn;

      hn++;
    }
  }
  return lu;
}

/****************************************************************/
/****************************************************************/
/* Functions that support LB_APPROACH = PHG_ROWS or PHG_COLS
 *
 * The objects to be partitioned are the matrix rows or columns,
 * and we are using Zoltan's PHG to do it.
 *
 * PHG_ROWS: the objects are the sparse matrix rows, and the
 *  hyperedges are the columns of the sparse matrix.
 *
 * PHG_COLS: the objects are the sparse matrix columns, and
 *  the hyperedges are the rows.
 */

ZOLTAN_NUM_OBJ_FN phg_rc_get_num_obj;
ZOLTAN_OBJ_LIST_FN phg_rc_get_obj_list;
ZOLTAN_HG_SIZE_CS_FN phg_rc_get_size_matrix;
ZOLTAN_HG_CS_FN phg_rc_get_matrix;
static int get_proc_part(ZOLTAN_MP_DATA *mpd, int num_export,
  unsigned int *gids, int *procs, int *parts,
  IJTYPE nobj, IJTYPE *objList, int *objProcs, int *objParts);
static int phg_rc_obj_to_proc(ZOLTAN_MP_DATA *mpd, IJTYPE objID);
static int phg_rc_my_objects(ZOLTAN_MP_DATA *mpd, IJTYPE *nobj, IJTYPE **objIDs);
/*
 *
 */
static int phg_rc_get_pins(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nPins, IJTYPE *I, IJTYPE *J, int *ijProcs, int *ijParts)
{
  if (mpd->approach == PHG_ROWS){
    /* The pin belongs to the partition its row belongs to */
    return phg_rc_get_rows(mpd, nPins, I, ijProcs, ijParts);
  }
  else{
    /* The pin belongs to the partition its column belongs to */
    return phg_rc_get_columns(mpd, nPins, J, ijProcs, ijParts);
  }
}
static int phg_rc_get_columns(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nCols, IJTYPE *colIDs, int *colProcs, int *colParts)
{
  int idx, i;
  char *yo = "phg_rc_get_columns";
  int nTotalCols = ((mpd->input_type == ROW_TYPE) ? mpd->numCR : mpd->numRC);

  if (mpd->approach != PHG_COLUMNS){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, 
     "chosen LB_APPROACH does not support obtaining column partitions\n");
    return ZOLTAN_FATAL; 
  }

  for (i=0; i<nCols; i++){ 
    idx = Zoltan_Lookup_Obj(mpd->col_lookup, colIDs[i], 0);
    if ((idx < 0) || (idx >= nTotalCols)){
      ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Unable to determine column partition\n");
     return ZOLTAN_FATAL; 
    }
    colProcs[i] = mpd->colproc[idx];
    colParts[i] = mpd->colpart[idx];
  }
  return ZOLTAN_OK;
}
static int phg_rc_get_rows(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nRows, IJTYPE *rowIDs, int *rowProcs, int *rowParts)
{
  int idx, i;
  char *yo = "phg_rc_get_rows";
  int nTotalRows = ((mpd->input_type == ROW_TYPE) ? mpd->numRC : mpd->numCR);

  if (mpd->approach != PHG_ROWS){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, 
     "chosen LB_APPROACH does not support obtaining row partitions\n");
    return ZOLTAN_FATAL; 
  }

  for (i=0; i<nRows; i++){ 
    idx = Zoltan_Lookup_Obj(mpd->row_lookup, rowIDs[i], 0);
    if ((idx < 0) || (idx >= nTotalRows)){
      ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Unable to determine row partition\n");
     return ZOLTAN_FATAL; 
    }
    rowProcs[i] = mpd->rowproc[idx];
    rowParts[i] = mpd->rowpart[idx];
  }
  return ZOLTAN_OK;
}
static int phg_rc_obj_to_proc(ZOLTAN_MP_DATA *mpd, IJTYPE objID)
{
  int nProcs = mpd->zzLib->Num_Proc;

  return objID % nProcs;
}
static int phg_rc_my_objects(ZOLTAN_MP_DATA *mpd, IJTYPE *nobj, IJTYPE **objIDs)
{
  IJTYPE i, nmyids=0, *obj=NULL;
  IJTYPE baseid, lastid;
  int ierr = ZOLTAN_OK;
  int me = mpd->zzLib->Proc;

  if (mpd->approach == PHG_ROWS){
    baseid = mpd->rowBaseID;
    lastid = mpd->rowBaseID +  mpd->nRows - 1;
  }
  else{
    baseid = mpd->colBaseID;
    lastid = mpd->colBaseID +  mpd->nCols - 1;
  }

  for (i=baseid; i<=lastid; i++){
    if (phg_rc_obj_to_proc(mpd, i) == me){
      nmyids++;
    }
  }

  obj = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * nmyids);
  if (nmyids && !obj){
    ierr = ZOLTAN_MEMERR;
  }
  else{
    *objIDs = obj;
    *nobj = nmyids;

    for (i=baseid; i<=lastid; i++){
      if (phg_rc_obj_to_proc(mpd, i) == me){
        *obj++ = i;
      }
    }
  }

  return ierr;
}

static int phg_rc_setup(ZOLTAN_MP_DATA *mpd)
{
  char *yo = "phg_rc_setup";
  int ierr = ZOLTAN_OK;

  Zoltan_Set_Param(mpd->zzLib, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(mpd->zzLib, "NUM_LID_ENTRIES", "1");

  Zoltan_Set_Param(mpd->zzLib, "LB_METHOD", "HYPERGRAPH");
  Zoltan_Set_Param(mpd->zzLib, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(mpd->zzLib, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(mpd->zzLib, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(mpd->zzLib, "ADD_OBJ_WEIGHT", "PINS");
  Zoltan_Set_Param(mpd->zzLib, "EDGE_WEIGHT_DIM", "0");

  /* Request export list that has proc/part for every object */
  Zoltan_Set_Param(mpd->zzLib, "RETURN_LISTS", "PARTITION");

  Zoltan_Set_Num_Obj_Fn(mpd->zzLib, phg_rc_get_num_obj, mpd);
  Zoltan_Set_Obj_List_Fn(mpd->zzLib, phg_rc_get_obj_list, mpd);
  Zoltan_Set_HG_Size_CS_Fn(mpd->zzLib, phg_rc_get_size_matrix, mpd);
  Zoltan_Set_HG_CS_Fn(mpd->zzLib, phg_rc_get_matrix, mpd);

  ierr = make_mirror(mpd);

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "make_mirror error\n");
  }
  else{
    ierr = phg_rc_my_objects(mpd, &mpd->nMyVtx, &mpd->vtxGID);
  }

  return ierr;
}

int phg_rc_get_num_obj(void *data, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  return mpd->nMyVtx;
}

void phg_rc_get_obj_list(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts, int *ierr)
{
  IJTYPE i;
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;
 
  for (i=0; i<mpd->nMyVtx; i++){
    gids[i] = (ZOLTAN_ID_TYPE)mpd->vtxGID[i];
    lids[i] = (ZOLTAN_ID_TYPE)i;
  }
}

void phg_rc_get_size_matrix(void *data, int *num_lists, int *num_pins,
  int *format, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  if (((mpd->approach == PHG_ROWS) && (mpd->input_type == ROW_TYPE)) ||
      ((mpd->approach == PHG_COLUMNS) && (mpd->input_type == COL_TYPE)) ){
    *num_lists = mpd->numCR;
    *num_pins  = mpd->mirrorPinIndex[mpd->numCR];
  }
  else{
    *num_lists = mpd->numRC;
    *num_pins  = mpd->pinIndex[mpd->numRC];
  }
  *format = ZOLTAN_COMPRESSED_EDGE;
}

void phg_rc_get_matrix(void *data, int num_gid_entries, int num_vtx_edge,
  int num_pins, int format, 
  ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, ZOLTAN_ID_PTR pin_GID, int *ierr)
{
  int i;
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;
  IJTYPE *edgeGID = (IJTYPE *)vtxedge_GID;
  IJTYPE *pinGID = (IJTYPE *)pin_GID;

  if (((mpd->approach == PHG_ROWS) && (mpd->input_type == ROW_TYPE)) ||
      ((mpd->approach == PHG_COLUMNS) && (mpd->input_type == COL_TYPE)) ){
    for (i=0; i<mpd->numCR; i++){
      edgeGID[i] = mpd->crGID[i];
      vtxedge_ptr[i] = mpd->mirrorPinIndex[i];
    }
    for (i=0; i< mpd->mirrorPinIndex[mpd->numCR]; i++){
      pinGID[i] = mpd->mirrorPinGID[i];
    }
  }
  else{
    for (i=0; i<mpd->numRC; i++){
      edgeGID[i] = mpd->rcGID[i];
      vtxedge_ptr[i] = mpd->pinIndex[i];
    }
    for (i=0; i< mpd->pinIndex[mpd->numRC]; i++){
      pinGID[i] = mpd->pinGID[i];
    }
  }
}

static int phg_rc_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  int *export_procs, int *export_to_part)
{
  char *yo = "phg_rc_result";
  int ierr = ZOLTAN_OK;
  int *proclist=NULL, *partlist=NULL;
  obj_lookup *lu=NULL;
  IJTYPE nobj=0;
  IJTYPE *objList=NULL;

  if (mpd->approach == PHG_ROWS){
    if (mpd->input_type == ROW_TYPE){
      nobj = mpd->numRC;
      objList = mpd->rcGID;
    }
    else{
      nobj = mpd->numCR;
      objList = mpd->rcGID;
    }
    lu = mpd->row_lookup = create_obj_lookup_table(nobj, objList);
    proclist = mpd->rowproc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
    partlist = mpd->rowpart = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  }
  else{
    if (mpd->input_type == ROW_TYPE){
      nobj = mpd->numCR;
      objList = mpd->crGID;
    }
    else{
      nobj = mpd->numRC;
      objList = mpd->rcGID;
    }
    lu = mpd->col_lookup = create_obj_lookup_table(nobj, objList);
    proclist = mpd->colproc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
    partlist = mpd->colpart = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  }

  if (nobj && (!lu || !proclist || !partlist)){
    free_obj_lookup_table(&lu);
    ZOLTAN_FREE(&proclist);
    ZOLTAN_FREE(&partlist);
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Out of memory.\n");
    ierr = ZOLTAN_MEMERR;
  }
  else{
    ierr = get_proc_part(mpd, num_export, 
         export_global_ids, export_procs, export_to_part, 
         nobj, objList, proclist, partlist);
  }

  return ierr;
}

static int get_proc_part(ZOLTAN_MP_DATA *mpd, int num_export,
  unsigned int *gids, int *procs, int *parts,
  IJTYPE nobj, IJTYPE *objList, int *objProcs, int *objParts)
{
  obj_lookup *myIDs=NULL;
  ZOLTAN_COMM_OBJ *plan=NULL;
  int *toProc=NULL, *outprocs=NULL, *outparts=NULL;
  IJTYPE *idReqs=NULL;
  int tag = 35000, nrecv=0, idx, i;
  MPI_Comm comm = mpd->zzLib->Communicator;
  int ierr = ZOLTAN_OK;

  /*
   * Create lookup table for my gids so I can reply to other procs.
   */
  myIDs = create_obj_lookup_table(num_export, gids);
  if (myIDs == NULL){
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /*
   * Create a communication plan, requesting new proc/part for
   * all IDs in my objList.
   */
  toProc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  if (nobj && !toProc){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<nobj; i++){
    toProc[i] = phg_rc_obj_to_proc(mpd, objList[i]);
  }

  ierr = Zoltan_Comm_Create(&plan, nobj, toProc, comm, tag, &nrecv);
  if (ierr != ZOLTAN_OK) goto End;
  /*
   * Send list of obj IDs I need proc/part for, get list from others.
   */
  idReqs = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * nrecv);
  outprocs = (int *)ZOLTAN_MALLOC(sizeof(int) * nrecv);
  outparts = (int *)ZOLTAN_MALLOC(sizeof(int) * nrecv);
  if (nrecv && (!idReqs || outprocs || outparts)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  ierr = Zoltan_Comm_Do(plan, --tag, (char *)objList, sizeof(IJTYPE),
                        (char *)idReqs);
  if (ierr != ZOLTAN_OK) goto End;

  /*
   * Create list of proc/part for all requests I received.
   */
  for (i=0; i<nrecv; i++){
    idx = Zoltan_Lookup_Obj(myIDs, idReqs[i], 0);
    if (idx == -1){
      ierr = ZOLTAN_FATAL;
      goto End;
    }
    outprocs[i] = procs[idx];
    outparts[i] = parts[idx];
  }

  /*
   * Send proc/parts to others and also get my requests.
   */
  ierr = Zoltan_Comm_Do_Reverse(plan, --tag, (char *)outprocs,
             sizeof(unsigned int), NULL, (char *)objProcs);
  if (ierr != ZOLTAN_OK) goto End;

  ierr = Zoltan_Comm_Do_Reverse(plan, --tag, (char *)outparts,
             sizeof(unsigned int), NULL, (char *)objParts);
  if (ierr != ZOLTAN_OK) goto End;

End:
  Zoltan_Comm_Destroy(&plan); 
  free_obj_lookup_table(&myIDs);
  ZOLTAN_FREE(&toProc);
  ZOLTAN_FREE(&idReqs);
  ZOLTAN_FREE(&outprocs);
  ZOLTAN_FREE(&outparts);

  return ierr;
}
/****************************************************************/
/****************************************************************/

/****************************************************************/
/****************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_MATRIX_PARTITION */
