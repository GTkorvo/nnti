/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  Thread Pool Interface (TPI).
 *
 *  A simple and miminalistic interface for executing subprograms
 *  in a thread parallel, shared memory mode.
 *
 *  States: the underlying thread pool has three states.
 *    1) Uninitialized: no extra threads exist, this is the initial state.
 *    2) Paused: extra threads exist and are blocked.
 *    3) Active: extra threads are calling the subprogram.
 *
 *  Threads are created and blocked ready-for-use so that the thread
 *  creation cost is paid only at initialization.
 *  The extra threads are blocked, as opposed to spinning, in the paused
 *  state so that they do not compete for compute cycles with other threads
 *  created external to the TPI interface.  For example, threads created
 *  by OpenMP.
 */

#ifndef ThreadPoolInterface_h
#define ThreadPoolInterface_h

#if defined( __cplusplus )
extern "C" {
#endif

/*--------------------------------------------------------------------*/
/** \brief  A utility to measure wall-clock time, which is frequently
 *          needed when performance testing HPC algorithms.
 */
double TPI_Walltime();

/*--------------------------------------------------------------------*/
/* All functions return zero for success. */

#define TPI_LOCK_BUSY      ((int)  1)  /**<  trylock or unlock failed */
#define TPI_ERROR_NULL     ((int) -1)  /**<  NULL input */
#define TPI_ERROR_SIZE     ((int) -2)  /**<  BAD input: size or index */
#define TPI_ERROR_LOCK     ((int) -3)  /**<  BAD lock or unlock */
#define TPI_ERROR_ACTIVE   ((int) -4)  /**<  BAD input: the pool is active  */
#define TPI_ERROR_INTERNAL ((int) -5)  /**< internal resource error */

/*--------------------------------------------------------------------*/

/** \brief  Work information passed to a work subprogram. */
struct TPI_Work_Struct {
  void * shared ;     /**<  Shared data input to TPI_Run */
  void * reduce ;     /**<  Data for reduce operation, if any */
  int    lock_count ; /**<  Count of locks requested via TPI_Run */
  int    work_count ; /**<  Count of work  requested via TPI_Run */
  int    work_rank ;  /**<  Rank of work fot the current call */
};

/** \brief  Typedef for work subprogram argument */
typedef const struct TPI_Work_Struct TPI_Work ;

/**  The interface for a parallel task */
typedef void (*TPI_work_subprogram)( TPI_Work * );

/**  The interface for a parallel reduction operation.
 *   Reduce the 'src' data into the 'dest' data.
 */
typedef
void (*TPI_reduce_subprogram)( void * dest , const void * src );

/** Run a work subprogram in thread parallel.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus a recursive call to TPI_Run is illegal.
 */
int TPI_Run( TPI_work_subprogram subprogram  ,
             void *              shared ,
             int                 work_count  ,
             int                 lock_count  );

/** Run a work and reduction subprograms in thread parallel.
 *  Each call to the work_subprogram has exclusive (thread safe)
 *  access to its work->reduce data.
 *  The reduce_subprogram has exclusive access to its input data.
 *  The reduce data given to the work_subprogram and reduce_subprogram
 *  may be the originally input reduce_data.  However, this data will
 *  most often be a memory-copy (memcpy) initialized version of the
 *  originally input reduce_data.
 */
int TPI_Run_reduce( TPI_work_subprogram   work_subprogram  ,
                    void *                work_shared ,
                    int                   work_count  ,
                    TPI_reduce_subprogram reduce_subprogram ,
                    void *                reduce_data ,
                    int                   reduce_size );

/** Run a work subprogram in thread parallel;
 *  run it exactly once on each thread.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus a recursive call to TPI_Run is illegal.
 */
int TPI_Run_threads( TPI_work_subprogram subprogram ,
                     void *              shared ,
                     int                 lock_count );

int TPI_Run_threads_reduce( TPI_work_subprogram   subprogram ,
                            void *                shared ,
                            TPI_reduce_subprogram reduce_subprogram ,
                            void *                reduce_data ,
                            int                   reduce_size );

/*--------------------------------------------------------------------*/
/** \brief  Blocks until lock lock_rank is obtained.
 *          The thread pool must be in the 'active' state.
 */
int TPI_Lock( int lock_rank );

/** \brief  Tries to lock lock_rank, returns 0 if successful. 
 *          The thread pool must be in the 'active' state.
 */
int TPI_Trylock( int lock_rank );

/** \brief  Unlocks lock lock_rank.
 *          The thread pool must be in the 'active' state.
 */
int TPI_Unlock( int lock_rank );

/*--------------------------------------------------------------------*/
/** Start up the requested number of threads, less the calling thread.
 *  Return the actual number of threads, including the calling thread,
 *  otherwise return an error.
 */
int TPI_Init( int thread_count );

/** Shut down all started threads. */
int TPI_Finalize();

/*--------------------------------------------------------------------*/

#if defined( __cplusplus )
}
#endif

#endif

