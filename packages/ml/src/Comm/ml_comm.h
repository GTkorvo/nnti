/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML communicator structure                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : September, 1998                                      */
/* ******************************************************************** */

#ifndef _MLCOMM_
#define _MLCOMM_

#include <stdio.h>
#include "ml_defs.h"
#include "ml_memory.h"

/* ******************************************************************** */
/* communicator can be loaded by users or use the resident mpi          */ 
/* -------------------------------------------------------------------- */

#ifdef ML_MPI
#include <mpi.h>
#define USR_COMM MPI_Comm
#define USR_REQ  MPI_Request
#else
#define USR_COMM int
#define USR_REQ  int
#endif

/* ******************************************************************** */
/* components of the ML communicator                                    */
/* -------------------------------------------------------------------- */

typedef struct ML_Comm_Struct
{
   int      ML_id;
   int      ML_mypid;
   int      ML_nprocs;
   USR_COMM USR_comm;
   int      (*USR_sendbytes)(void*,unsigned int,int,int,USR_COMM);
   int      (*USR_irecvbytes)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*);
   int      (*USR_waitbytes)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*);

} ML_Comm;

extern ML_Comm *global_comm; /* should be made obsolete */

/* ******************************************************************** */
/* functions for the ML communicator                                    */
/* -------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C"
{
#endif

extern int    ML_Comm_Create( ML_Comm ** comm );
extern int    ML_Comm_Destroy( ML_Comm ** comm );
extern int    ML_Comm_Check( ML_Comm *comm );

extern int    ML_Comm_Set_UsrComm( ML_Comm *comm, USR_COMM com );
extern int    ML_Comm_Set_Mypid( ML_Comm *comm, int mypid );
extern int    ML_Comm_Set_Nprocs( ML_Comm *comm, int nprocs);
extern int    ML_Comm_Set_SendFcn( ML_Comm *comm, int (*SendFcn)(void*,unsigned int,int,int,USR_COMM));
extern int    ML_Comm_Set_RecvFcn( ML_Comm *comm, int (*RecvFcn)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*));
extern int    ML_Comm_Set_WaitFcn( ML_Comm *comm, int (*WaitFcn)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*));

extern int    ML_Comm_GmaxInt( ML_Comm *comm, int intdata );
extern double ML_Comm_GmaxDouble(ML_Comm *comm, double ddata );
extern int    ML_Comm_GsumInt( ML_Comm *comm, int intdata );
extern double ML_Comm_GsumDouble(ML_Comm *comm, double ddata );
extern int    ML_Comm_GsumVecInt(ML_Comm *,int *invec,int *tmpvec,int vleng);
extern int    ML_Comm_GappendInt(ML_Comm *,int *invec,int *lleng, int tleng);
extern int    ML_Comm_GappendDouble(ML_Comm*,double *dvec,int *lleng,int tleng);

extern int    ML_Comm_Irecv(void*,unsigned int,int *,int *,USR_COMM,USR_REQ*);
extern int    ML_Comm_Wait (void*,unsigned int,int *,int *,USR_COMM,USR_REQ*);
extern int    ML_Comm_Send (void*,unsigned int,int,  int,  USR_COMM );

#ifdef __cplusplus
}
#endif

#endif

