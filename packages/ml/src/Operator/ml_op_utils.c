/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_struct.h"
#include "ml_op_utils.h"

/* ******************************************************************** */
/* Blow away any inter-mixing between boundary and interior points in   */
/* the matrix ml_handle->Pmat[level2].                                  */
/* -------------------------------------------------------------------- */

int oldML_Mdfy_Prolongator_DirBdry(ML *ml_handle, int level2, 
                                   double *boundary, double *cboundary)
{
   int *cols, *row_ptr, i, j;
   double *vals;
   struct ML_CSR_MSRdata *temp;

   int size;

   if (ml_handle->Pmat[level2].getrow->external != CSR_getrows)
     perror("ML_Mdfy_Prolongator_DirBdry can only be used with CSR matrices\n");

   temp    = (struct ML_CSR_MSRdata *) ml_handle->Pmat[level2].data;
   cols    = temp->columns;
   vals    = temp->values;
   row_ptr = temp->rowptr;
   size    = ml_handle->Pmat[level2].outvec_leng;

   /* zero out any inter-mixing between boundary  */
   /* and interior in the interpolation operator. */

   for (i = 0; i < size; i++) {
      if (boundary[i] == 1.) {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( cboundary[cols[j]] == 0.) vals[j] = 0.0;
      }
      else {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( cboundary[cols[j]] == 1.) vals[j] = 0.0;
      }
   }
   return(1);
}

/* ******************************************************************** */
/* Blow away any inter-mixing between boundary and interior points in   */
/* the matrix ml_handle->Pmat[level2].                                  */
/* -------------------------------------------------------------------- */

int ML_Mdfy_Prolongator_DirBdry(ML *ml_handle, int level2, int size,
     int fine_size )
{
   int fboundary_length,cboundary_length, *fboundary_list, *cboundary_list;
   char *f_bdry;
   int *cols, *row_ptr, i, j, jk;
   double *vals;
   struct ML_CSR_MSRdata *temp;
   double *dtemp;
   ML_CommInfoOP *comm_info;

   comm_info        = ml_handle->Pmat[level2].getrow->pre_comm;
   fboundary_length = ml_handle->Pmat[level2].to->BCs->Dirichlet_grid_length;
   fboundary_list   = ml_handle->Pmat[level2].to->BCs->Dirichlet_grid_list;
   cboundary_length = ml_handle->Pmat[level2].from->BCs->Dirichlet_grid_length;
   cboundary_list   = ml_handle->Pmat[level2].from->BCs->Dirichlet_grid_list;
   dtemp       = (double *) malloc((size+1)*sizeof(double));
   f_bdry      = (char *) malloc((fine_size+1)*sizeof(char));
   if (f_bdry == NULL) {
        printf("No space to compute coarse boundary\n");
        exit(1);
   }
   for (jk = 0; jk < fine_size; jk++)       f_bdry[jk] = 'i';
   for (jk = 0; jk < fboundary_length; jk++) f_bdry[fboundary_list[jk]] = 'b';
   for (jk = 0; jk < size; jk++)       dtemp[jk] = 0.0;
   for (jk = 0; jk < cboundary_length; jk++) dtemp[cboundary_list[jk]] = 1.;
   if ( comm_info != NULL)
      ML_exchange_bdry(dtemp, comm_info, size, ml_handle->comm, ML_OVERWRITE);



   if (ml_handle->Pmat[level2].getrow->external != CSR_getrows)
     perror("ML_Mdfy_Prolongator_DirBdry can only be used with CSR matrices\n");

   temp    = (struct ML_CSR_MSRdata *) ml_handle->Pmat[level2].data;
   cols    = temp->columns;
   vals    = temp->values;
   row_ptr = temp->rowptr;

   /* zero out any inter-mixing between boundary  */
   /* and interior in the interpolation operator. */

   for (i = 0; i < fine_size; i++) {
      if (f_bdry[i] == 'b') {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( dtemp[cols[j]] == 0.0) vals[j] = 0.0;
      }
      else {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( dtemp[cols[j]] == 1.0) vals[j] = 0.0;
      }
   }
   return(1);
}

/* ******************************************************************** */
/* -------------------------------------------------------------------- */

int ML_Compute_Coarse_Bdry(ML *ml_handle, int level, int size, int fine_size)
{
   int     *cols, *row_ptr, jk, Ncoarse;
   struct  ML_CSR_MSRdata *temp;
   int     *boundary_list, boundary_length, *cboundary_list, cboundary_length = 0;
   char    *f_bdry, *c_bdry;

   Ncoarse     = ml_handle->Pmat[level].invec_leng;
   c_bdry      = (char *) malloc((size+1)*sizeof(char));
   f_bdry      = (char *) malloc((fine_size+1)*sizeof(char));
   if (f_bdry == NULL) {
      printf("No space to compute coarse boundary\n");
      exit(1);
   }
   boundary_length = ml_handle->Pmat[level].to->BCs->Dirichlet_grid_length;
   boundary_list   = ml_handle->Pmat[level].to->BCs->Dirichlet_grid_list;
   for (jk = 0; jk < fine_size; jk++) f_bdry[jk] = 'i';
   for (jk = 0; jk < boundary_length; jk++) f_bdry[boundary_list[jk]] = 'b';

   /* Mark the coarse grid boundary */

   temp        = (struct ML_CSR_MSRdata*) ml_handle->Pmat[level].data;
   cols        = temp->columns;
   row_ptr     = temp->rowptr;

   for (jk = 0; jk < size; jk++) c_bdry[jk] = 'i';
   for (jk = 0; jk < fine_size; jk++) {
      if ((row_ptr[jk+1] - row_ptr[jk] == 1) && ( f_bdry[jk] == 'b')) {
         c_bdry[ cols[row_ptr[jk]]] = 'b';
      }
   }

   /* stuff the coarse boundary information into ML */

   for (jk = 0; jk < Ncoarse; jk++) {
      if (c_bdry[jk] == 'b') cboundary_length++;
   }
   cboundary_list = (int *) malloc((cboundary_length+1)*sizeof(int));
   if (cboundary_list == NULL) {
      printf("No space to compute coarse boundary\n");
      exit(1);
   }
   cboundary_length = 0;
   for (jk = 0; jk < Ncoarse; jk++) {
      if (c_bdry[jk] == 'b') cboundary_list[cboundary_length++] = jk;
   }
   ML_Set_BoundaryTypes(ml_handle, ml_handle->Pmat[level].from->levelnum,
                        ML_BDRY_DIRICHLET,cboundary_length,cboundary_list);

   free(c_bdry);
   free(f_bdry);
   free(cboundary_list);
   return(1);
}

/* ******************************************************************** */
/*
 * Take the 'matrix' defined by getrow() and vec_comm() and make a CSR
 * copy of it which will be stored in the matrix ml_handle->Pmat[level2].
 *
 *
 * Parameters
 * ==========
 *
 *    isize      On input, number of local columns in matrix to be copied.
 *
 *    osize      On input, number of local rows in matrix to be copied.
 *
 *    getrow()   On input, user's function to get rows of the matrix.
 *
 *    vec_comm() On input, user's function to update a vector via communication
 *               so that a local matrix-vector product can occur.
 *
 * -------------------------------------------------------------------- */

int ML_Gen_Prolongator_Getrow(ML *ml_handle, int level2, int level, int isize,
   int osize, int (*getrow)(void* , int , int *, int , int *, double *, int *),
   int (*vec_comm)(double *, void*), void *data, int Nghost)
{
   int *cols, *row_ptr, space, flag, nz_ptr, i, length;
   double *vals;
   double dsize, di;
   struct ML_CSR_MSRdata *temp;

   row_ptr = (int *) malloc(sizeof(int)*(osize+1));

   space = osize*5+30;

   flag = 0;

   while (flag == 0) {
      cols    = (int    *) malloc(sizeof(int)*space);
      vals    = (double *) malloc(sizeof(double)*space);

      nz_ptr = 0;
      row_ptr[0] = nz_ptr;
      for (i = 0; i < osize; i++) {
         flag = getrow(data, 1, &i, space-nz_ptr, &(cols[nz_ptr]),
                       &(vals[nz_ptr]), &length);
         if (flag == 0) break;
         nz_ptr += length;
         row_ptr[i+1] = nz_ptr;
      }
      if (flag == 0) {
         dsize = (double) osize;
         di    = (double) (i+1);
         dsize = 1.2*dsize/di;
         space = (int) ( ((double) space)*dsize);
         space++;
         ML_free(vals);
         ML_free(cols);
      }
   }

   /* store the matrix into ML */

   temp = (struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
   temp->columns = cols;
   temp->values  = vals;
   temp->rowptr = row_ptr;

   ml_handle->Pmat[level2].data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Init_Prolongator(ml_handle, level2, level, isize, osize, (void *) temp);
   ML_Operator_Set_ApplyFunc(&(ml_handle->Pmat[level2]),ML_INTERNAL,CSR_matvec);
/*
printf("we've changed the data pointer ? ....\n");
*/
   if (vec_comm != NULL) {
      ML_CommInfoOP_Generate(&(ml_handle->Pmat[level2].getrow->pre_comm), 
			  vec_comm, data, ml_handle->comm, 
			  ml_handle->Pmat[level2].invec_leng, Nghost);
   }
   else ml_handle->Pmat[level2].getrow->pre_comm = NULL;

   ML_Operator_Set_Getrow(&(ml_handle->Pmat[level2]), ML_EXTERNAL, 
			  ml_handle->Pmat[level2].outvec_leng, CSR_getrows);

/*
   ML_CommInfoOP_Generate(&(ml_handle->Pmat[level2].getrow->pre_comm),
                   vec_comm, data, ml_handle->comm, 
		   ml_handle->Pmat[level2].invec_leng, Nghost);
*/

   return(1);
}

/* ******************************************************************** */
/* Take the Prolongator given by ml_handle->Pmat[level2], transpose
 * it and store it into the matrix given by ml_handle->Rmat[level].
 *
 * -------------------------------------------------------------------- */

int ML_Gen_Restrictor_TransP(ML *ml_handle, int level, int level2)
{

   ML_Operator *Pmat, *Rmat;
   int  *row_ptr, *colbuf, *cols;
   int isize, osize, i, j, N_nzs, flag, length, sum, new_sum;
   int Nneighbors, *neigh_list, *send_list, *rcv_list, Nsend, Nrcv;
   void *data = NULL;
   double *valbuf, *vals;
   int (*getrow)(void* , int , int *, int , int *, double *, int *) = NULL;
   struct ML_CSR_MSRdata *temp;
   int Nghost = 0, Nghost2 = 0;
   int *remap, remap_leng;
   ML_CommInfoOP *c_info, **c2_info;
   char procname[128];

   sprintf(procname,"p%d",ml_handle->comm->ML_mypid);

   /* pull out things from ml_handle */


   Pmat  = &(ml_handle->Pmat[level2]);
   temp = (struct ML_CSR_MSRdata *) Pmat->data;
   Rmat  = &(ml_handle->Rmat[level]);
   isize = Pmat->outvec_leng;
   osize = Pmat->invec_leng;
   if (Pmat->getrow->ML_id == ML_EXTERNAL) {
       data   = Pmat->data;
       getrow = Pmat->getrow->external;
   }
   else if (Pmat->getrow->ML_id == ML_INTERNAL) {
       data   = (void *) Pmat;
       getrow = Pmat->getrow->internal;
   }
   else perror("ML_Gen_Restrictor_TransP: Getrow not defined for P!\n");

   /* transpose Pmat's communication list. This means that PRE communication */
   /* is replaced by POST, ML_OVERWRITE is replaced by ML_ADD, and the send  */
   /* send and receive lists are swapped.                                    */

   c_info     = Pmat->getrow->pre_comm;
   Nneighbors = ML_CommInfoOP_Get_Nneighbors(c_info);
   neigh_list = ML_CommInfoOP_Get_neighbors(c_info);
   remap_leng = osize;
   Nrcv = 0;
   Nsend = 0;
   for (i = 0; i < Nneighbors; i++) {
      Nrcv  += ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nsend += ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
   }
   remap_leng = osize + Nrcv + Nsend;
   remap = (int *) malloc( remap_leng*sizeof(int));
   for (i = 0; i < osize; i++) remap[i] = i;
   for (i = osize; i < osize+Nrcv+Nsend; i++) 
      remap[i] = -1;
 
   c2_info     = &(Rmat->getrow->post_comm);
   ML_CommInfoOP_Set_neighbors(c2_info, Nneighbors,
 			      neigh_list,ML_ADD,remap,remap_leng);
   ML_free(remap);
   for (i = 0; i < Nneighbors; i++) {
      Nsend      = ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
      send_list  = ML_CommInfoOP_Get_sendlist (c_info, neigh_list[i]);
      Nrcv       = ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nghost    += Nrcv;
      rcv_list   = ML_CommInfoOP_Get_rcvlist(c_info, neigh_list[i]);
      /* handle empty rows ... i.e. ghost variables not used */
      if (rcv_list != NULL) {
         for (j = 0; j < Nrcv; j++) {
            if (rcv_list[j] > Nghost2 + osize - 1)
               Nghost2 = rcv_list[j] - osize + 1;
         }
      }
 
      ML_CommInfoOP_Set_exch_info(*c2_info, neigh_list[i], Nsend, send_list,
 				 Nrcv,rcv_list);
      if (send_list != NULL) ML_free(send_list);
      if ( rcv_list != NULL) ML_free( rcv_list);
   }
   if (Nghost2 > Nghost) Nghost = Nghost2;
   if (neigh_list != NULL) ML_free(neigh_list);

   row_ptr = (int    *) malloc(sizeof(int)*(Nghost+osize+1));
   colbuf  = (int    *) malloc(sizeof(int)*(Nghost+osize+1));
   valbuf  = (double *) malloc(sizeof(double)*(Nghost+osize+1));

   /* count the total number of nonzeros and compute */
   /* the length of each row in the transpose.       */
 
   for (i = 0; i < Nghost+osize; i++) row_ptr[i] = 0;

   N_nzs = 0;
   for (i = 0; i < isize; i++) {
      flag = getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      if (flag == 0) pr_error("ML_Transpose_Prolongator: sizes don't work\n");
      N_nzs += length;
      for (j = 0; j < length; j++)
         row_ptr[  colbuf[j] ]++;
   }

   cols    = (int    *) malloc(sizeof(int   )*(N_nzs+1));
   vals    = (double *) malloc(sizeof(double)*(N_nzs+1));
   if (vals == NULL) 
      pr_error("ML_Gen_Restrictor_TransP: Out of space\n");

   /* set 'row_ptr' so it points to the beginning of each row */

   sum = 0;
   for (i = 0; i < Nghost+osize; i++) {
      new_sum = sum + row_ptr[i];
      row_ptr[i] = sum;
      sum = new_sum;
   }
   row_ptr[osize+Nghost] = sum;

   /* read in the prolongator matrix and store transpose in Rmat */

   for (i = 0; i < isize; i++) {
      getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      for (j = 0; j < length; j++) {
         cols[ row_ptr[ colbuf[j] ]   ] = i;
         vals[ row_ptr[ colbuf[j] ]++ ] = valbuf[j];
      }
   }

   /* row_ptr[i] now points to the i+1th row.    */
   /* Reset it so that it points to the ith row. */

   for (i = Nghost+osize; i > 0; i--)
      row_ptr[i] = row_ptr[i-1];
   row_ptr[0] = 0;

   ML_free(valbuf);
   ML_free(colbuf);

   /* store the matrix into ML */

   temp = (struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
   temp->columns = cols;
   temp->values  = vals;
   temp->rowptr  = row_ptr;
 
   ml_handle->Rmat[level].data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Init_Restrictor(ml_handle, level, level2, isize, osize, (void *) temp);
   ML_Operator_Set_ApplyFunc(Rmat,ML_INTERNAL, CSR_matvec);
   ML_Operator_Set_Getrow(&(ml_handle->Rmat[level]), ML_EXTERNAL,
                                 Nghost+osize, CSR_getrows);
  return(1);
}
/*  SYMB_GRID::partitionBlocksNodes ****************************************
 *
 *   - do local partition
 *
 *   INPUT:
 *    - nblk: number of blocks
 *    - pnode_part[nLocalNd]: local nodes to block map (out)
 *
 */
#ifdef METIS
#include "metis.h"
#endif
int ML_Operator_BlockPartition(ML_Operator *matrix, int nLocalNd, int *nblk,
                         int *pnode_part,
                         int *ndwts /*=NULL*/, int *egwts/*=NULL*/,
                         int nedges /*= 0*/ )
{
#ifdef METIS
  int locid, ii, numadjac, *bindx = NULL;
  idxtype *xadj, *adjncy, *blks;
  int options[5]={0,3,1,1,0};
  int weightflag = ndwts ? 2 : 0;
  int nmbng = 0, edgecut = -1, n = nLocalNd, np = *nblk;
  double *val = NULL;
  int allocated = 0, row_length, j;

  if( egwts ) weightflag++;

  if( *nblk == 1 || nLocalNd < 1 ){
    for( ii = 0 ; ii < nLocalNd ; ii++ ) pnode_part[ii] = 0;
    return 0;
  }

  /* set 'xadj' & 'adjncy' adjacentcy data */

  xadj = (idxtype *) ML_allocate( (nLocalNd+1) * sizeof(idxtype) );
  if (xadj == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");

  ii = 0;
  for( locid = 0 ; locid < nLocalNd ; locid++ ) {
    ML_get_matrix_row(matrix, 1, &locid, &allocated, &bindx, &val,
                      &row_length, 0);

    for (j = 0; j < row_length; j++) {
       if ( bindx[j] < nLocalNd) ii++;
    }
  }
  numadjac = ii;
  adjncy = (idxtype *) ML_allocate( (numadjac+1) * sizeof(idxtype) );
  if (adjncy == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");

  ii = 0;
  for( locid = 0 ; locid < nLocalNd ; locid++ ) {
    xadj[locid] = ii;
    ML_get_matrix_row(matrix, 1, &locid, &allocated, &bindx,
                      &val, &row_length, 0);

    for (j = 0; j < row_length; j++) {
       if ( bindx[j] < nLocalNd) adjncy[ii++] = bindx[j];
    }
  }
  xadj[nLocalNd] = ii;
  ML_free(val);
  ML_free(bindx);

  /* get local partition */

  METIS_PartGraphKway( &n, xadj, adjncy, ndwts, egwts, &weightflag, &nmbng,
                       &np, options, &edgecut, pnode_part );
  ML_free(xadj);
  ML_free(adjncy);

  blks = (int *) ML_allocate((np+1)*sizeof(int));
  if (blks == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");
  for (j = 0; j < *nblk; j++) blks[j] = -1;
  for (j = 0; j < n; j++) blks[pnode_part[j]] = -2;
  ii = 0;
  for (j = 0; j < *nblk; j++) {
    if ( blks[j] == -2) {
       blks[j] = ii;
       ii++;
    }
  }
  for (j = 0; j < n; j++) pnode_part[j] = blks[pnode_part[j]];
  *nblk = ii;
  
  ML_free(blks);

#else
  printf("ML_partitionBlocksNodes: Metis not linked\n");
#endif

  return 0;
}
/* ******************************************************************** */
/* Take the Transpose of an ML_Operator and stick it in a new matrix    */
/* -------------------------------------------------------------------- */

int ML_Operator_Transpose(ML_Operator *Amat, ML_Operator *Amat_trans )
{

   int  *row_ptr, *colbuf, *cols;
   int isize, osize, i, j, N_nzs, flag, length, sum, new_sum;
   int Nneighbors, *neigh_list, *send_list, *rcv_list, Nsend, Nrcv;
   void *data = NULL;
   double *valbuf, *vals;
   int (*getrow)(void* , int , int *, int , int *, double *, int *) = NULL;
   struct ML_CSR_MSRdata *temp;
   int Nghost = 0, Nghost2 = 0;
   int *remap, remap_leng;
   ML_CommInfoOP *c_info, **c2_info;

   /* pull out things from ml_handle */


   temp = (struct ML_CSR_MSRdata *) Amat->data;
   isize = Amat->outvec_leng;
   osize = Amat->invec_leng;
   if (Amat->getrow->ML_id == ML_EXTERNAL) {
       data   = Amat->data;
       getrow = Amat->getrow->external;
   }
   else if (Amat->getrow->ML_id == ML_INTERNAL) {
       data   = (void *) Amat;
       getrow = Amat->getrow->internal;
   }
   else perror("ML_Operator_Transpose: Getrow not defined for P!\n");

   /* transpose Amat's communication list. This means that PRE communication */
   /* is replaced by POST, ML_OVERWRITE is replaced by ML_ADD, and the send  */
   /* send and receive lists are swapped.                                    */

   c_info     = Amat->getrow->pre_comm;
   Nneighbors = ML_CommInfoOP_Get_Nneighbors(c_info);
   neigh_list = ML_CommInfoOP_Get_neighbors(c_info);
   remap_leng = osize;
   Nrcv = 0;
   Nsend = 0;
   for (i = 0; i < Nneighbors; i++) {
      Nrcv  += ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nsend += ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
   }
   remap_leng = osize + Nrcv + Nsend;
   remap = (int *) malloc( remap_leng*sizeof(int));
   for (i = 0; i < osize; i++) remap[i] = i;
   for (i = osize; i < osize+Nrcv+Nsend; i++) 
      remap[i] = -1;
 
   c2_info     = &(Amat_trans->getrow->post_comm);
   ML_CommInfoOP_Set_neighbors(c2_info, Nneighbors,
 			      neigh_list,ML_ADD,remap,remap_leng);
   ML_free(remap);
   for (i = 0; i < Nneighbors; i++) {
      Nsend      = ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
      send_list  = ML_CommInfoOP_Get_sendlist (c_info, neigh_list[i]);
      Nrcv       = ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nghost    += Nrcv;
      rcv_list   = ML_CommInfoOP_Get_rcvlist(c_info, neigh_list[i]);
      /* handle empty rows ... i.e. ghost variables not used */
      if (rcv_list != NULL) {
         for (j = 0; j < Nrcv; j++) {
            if (rcv_list[j] > Nghost2 + osize - 1)
               Nghost2 = rcv_list[j] - osize + 1;
         }
      }
 
      ML_CommInfoOP_Set_exch_info(*c2_info, neigh_list[i], Nsend, send_list,
 				 Nrcv,rcv_list);
      if (send_list != NULL) ML_free(send_list);
      if ( rcv_list != NULL) ML_free( rcv_list);
   }
   if (Nghost2 > Nghost) Nghost = Nghost2;
   if (neigh_list != NULL) ML_free(neigh_list);

   row_ptr = (int    *) malloc(sizeof(int)*(Nghost+osize+1));
   colbuf  = (int    *) malloc(sizeof(int)*(Nghost+osize+1));
   valbuf  = (double *) malloc(sizeof(double)*(Nghost+osize+1));

   /* count the total number of nonzeros and compute */
   /* the length of each row in the transpose.       */
 
   for (i = 0; i < Nghost+osize; i++) row_ptr[i] = 0;

   N_nzs = 0;
   for (i = 0; i < isize; i++) {
      flag = getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      if (flag == 0) perror("ML_Transpose_Prolongator: sizes don't work\n");
      N_nzs += length;
      for (j = 0; j < length; j++)
         row_ptr[  colbuf[j] ]++;
   }

   cols    = (int    *) malloc(sizeof(int   )*(N_nzs+1));
   vals    = (double *) malloc(sizeof(double)*(N_nzs+1));
   if (vals == NULL) 
      pr_error("ML_Gen_Restrictor_TransP: Out of space\n");

   /* set 'row_ptr' so it points to the beginning of each row */

   sum = 0;
   for (i = 0; i < Nghost+osize; i++) {
      new_sum = sum + row_ptr[i];
      row_ptr[i] = sum;
      sum = new_sum;
   }
   row_ptr[osize+Nghost] = sum;

   /* read in the prolongator matrix and store transpose in Amat_trans */

   for (i = 0; i < isize; i++) {
      getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      for (j = 0; j < length; j++) {
         cols[ row_ptr[ colbuf[j] ]   ] = i;
         vals[ row_ptr[ colbuf[j] ]++ ] = valbuf[j];
      }
   }

   /* row_ptr[i] now points to the i+1th row.    */
   /* Reset it so that it points to the ith row. */

   for (i = Nghost+osize; i > 0; i--)
      row_ptr[i] = row_ptr[i-1];
   row_ptr[0] = 0;

   ML_free(valbuf);
   ML_free(colbuf);

   /* store the matrix into ML */

   temp = (struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
   temp->columns = cols;
   temp->values  = vals;
   temp->rowptr  = row_ptr;
   Amat_trans->data_destroy = ML_CSR_MSRdata_Destroy;
   
   ML_Operator_Set_ApplyFuncData(Amat_trans, isize, osize, 
                                  ML_EMPTY, temp, osize, NULL, 0);
   ML_Operator_Set_ApplyFunc(Amat_trans,ML_INTERNAL, CSR_matvec);
   ML_Operator_Set_Getrow(Amat_trans, ML_EXTERNAL,
                                 Nghost+osize, CSR_getrows);
  return(1);
}

/************************************************************************/
/* Take a matrix that is effectively partitioned by columns and         */
/* transform it into one that is partitioned by rows. The original      */
/* matrix was most likely created by transposing a matrix partitioned   */
/* by row.                                                              */
/*----------------------------------------------------------------------*/

int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans)
{ 
 
  ML_Operator *eye1, *eye2;
 
  eye1 = ML_Operator_Create(A->comm);
  eye2 = ML_Operator_Create(A->comm);
 
  ML_Operator_Set_ApplyFuncData(eye1, A->invec_leng, A->invec_leng,
            ML_EXTERNAL,NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye1, ML_EXTERNAL, A->invec_leng, eye_getrows);
 
  ML_Operator_Set_ApplyFuncData(eye2, A->invec_leng, A->invec_leng,
            ML_EXTERNAL,NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye2, ML_EXTERNAL, A->invec_leng, eye_getrows);
  ML_2matmult(A, eye1, Atrans);
 
  return 1;
}  


#ifdef new2row
int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans)
{
  int         max_per_proc;
  ML_Operator *Acomm, *tptr;
 
  if (A->getrow->use_loc_glob_map == ML_YES)
     pr_error("ML_Operator_ColPartition2RowPartition: Matrix already has local"
              "column indices mapped to global indices\n");
  if (A->getrow->pre_comm != NULL)
     pr_error("ML_Operator_ColPartition2RowPartiion: Matrix has a"
              "pre-communication structure?\n");
 
  ML_create_unique_col_id(A->invec_leng, &(A->getrow->loc_glob_map),
                           NULL, &max_per_proc, A->comm);
 
  A->getrow->use_loc_glob_map = ML_YES;
 
  if (A->getrow->post_comm != NULL)
      ML_exchange_rows( A, &Acomm, A->getrow->post_comm);
  else Acomm = A;
 
  ML_back_to_csrlocal(Acomm, Atrans, max_per_proc);
 
  ML_free(A->getrow->loc_glob_map); A->getrow->loc_glob_map = NULL;
  A->getrow->use_loc_glob_map = ML_NO;
 
  if (A->getrow->post_comm != NULL) {
      tptr = Acomm;
      while ( (tptr!= NULL) && (tptr->sub_matrix != A))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(Acomm);
      ML_Operator_Destroy(Acomm);
   }
 
  return 1;
} 
#endif


/************************************************************************/
/* Getrow function for the identity matrix.                             */
/*----------------------------------------------------------------------*/

int eye_getrows(void *data, int N_requested_rows, int requested_rows[],
                int allocated_space, int columns[], double values[],
				int row_lengths[])
{
   double *temp;
   int    i;

   temp = (double *) data;

   if (allocated_space < N_requested_rows) return(0);

   for (i = 0; i < N_requested_rows; i++) {
      row_lengths[i] = 1;
      columns[i]     = requested_rows[i];
      values[i]      = 1.;
   }
   return(1);
}

/************************************************************************/
/* Matvec function for the identity matrix.                             */
/*----------------------------------------------------------------------*/

int eye_matvec(void *Amat_in, int ilen, double p[], int olen, double ap[])
{
  int i;

  for (i = 0; i < olen; i++) ap[i] = p[i];

  return(1);
}

/************************************************************************/
/* Take the transpose of an ML_Operator and realign resulting matrix    */
/* so that it is partitioned by rows.                                   */
/*----------------------------------------------------------------------*/

int ML_Operator_Transpose_byrow(ML_Operator *A, ML_Operator *Atrans)
{
  ML_Operator *temp;

  temp = ML_Operator_Create(A->comm);
  ML_Operator_Transpose(A, temp);
  ML_Operator_ColPartition2RowPartition(temp, Atrans);
  return 1;
}
