// *****************************************************************************
// * Zoltan Library for Parallel Applications                                  *
// * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
// * This software is distributed under the GNU Lesser General Public License. *
// * For more info, see the README file in the top-level Zoltan directory.     *
// *****************************************************************************
// *****************************************************************************
// * CVS File Information :
// *    $RCSfile$
// *    $Author$
// *    $Date$
// *    $Revision$
// *****************************************************************************
//
// ************************************************************************
// 
// C++ wrappers for Zoltan communication library.  
//
// Two styles of initialization: 
//
//   C++ style:  Zoltan_Comm comm(nvals, assign, comm, tag, pnvals_recv);
//
//   C style: Zoltan_Comm comm;
//            comm.Create(nvals, assign, comm, tag, pnvals_recv);
// 
// ************************************************************************

#ifndef ZOLTAN_COMM_CPP_H_
#define ZOLTAN_COMM_CPP_H_

#include "comm.h"

class Zoltan_Comm {

public:
  
  Zoltan_Comm(int nvals, int *assign, MPI_Comm comm, int tag, int *pnvals_recv)
    {
    // Assumption: MPI has been initialized prior to this call.
    Zoltan_Comm_Create(&this->Plan, nvals, assign, comm, tag, pnvals_recv);
    }

  Zoltan_Comm()
    {
    this->Plan = NULL;
    // Caller will have to call Create to finish initialization of object
    }

  int Create(int nvals, int *assign, MPI_Comm comm, int tag, int *pnvals_recv)
    {
    if (this->Plan)
      {
      Zoltan_Comm_Destroy(&this->Plan);
      this->Plan = NULL;
      }

    int rc = Zoltan_Comm_Create(&this->Plan, nvals, assign, comm, tag, pnvals_recv);

    return rc;
    }

  ~Zoltan_Comm()
    {
    Zoltan_Comm_Destroy(&this->Plan);
    }
      
  int Resize(int *sizes, int tag, int *sum_recv_sizes)
    {
    return Zoltan_Comm_Resize( this->Plan, sizes, tag, sum_recv_sizes);
    }
  
  int Do(int tag, char *send_data, int nbytes, char *recv_data)
    {
    return Zoltan_Comm_Do(this->Plan, tag, send_data, nbytes, recv_data);
    }

  int Do_Post( int tag, char *send_data, int nbytes, char *recv_data)
    {
    return Zoltan_Comm_Do_Post(this->Plan, tag, send_data, nbytes, recv_data);
    }

  int Do_Wait(int tag, char *send_data, int nbytes, char *recv_data)
    {
    return Zoltan_Comm_Do_Wait(this->Plan, tag, send_data, nbytes, recv_data);
    }
  
  int Do_Reverse(int tag, char *send_data, int nbytes, int *sizes, char *recv_data)
    {
    return Zoltan_Comm_Do_Reverse(this->Plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }

  int Do_Reverse_Post(int tag, char *send_data, int nbytes, int *sizes, 
    char *recv_data)
    {
    return Zoltan_Comm_Do_Reverse_Post(this->Plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }

  int Do_Reverse_Wait(int tag, char *send_data, int nbytes, int *sizes, 
    char *recv_data)
    {
    return Zoltan_Comm_Do_Reverse_Wait(this->Plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }
  
  int Info( int *nsends, int *send_procs,
    int *send_lengths, int *send_nvals, int *send_max_size, int *send_list,
    int *nrecvs, int *recv_procs, int *recv_lengths, int *recv_nvals,
    int *recv_total_size, int *recv_list, int *self_msg)
    {
      return Zoltan_Comm_Info( this->Plan, nsends, send_procs, send_lengths, 
        send_nvals, send_max_size, send_list, nrecvs, recv_procs, recv_lengths, 
        recv_nvals, recv_total_size, recv_list, self_msg);
    }
  
  int Invert_Plan()
    {
    return Zoltan_Comm_Invert_Plan(&this->Plan);
    }

  // For debugging and testing purposes

  ZOLTAN_COMM_OBJ *Get_ZOLTAN_COMM_OBJ(){return this->Plan;} 

  // Static methods

  static int Invert_Map( int *lengths_to, int *procs_to, int nsends, int self_msg,
    int **plengths_from, int **pprocs_from, int *pnrecvs, int my_proc,
    int nprocs, int out_of_mem, int tag, MPI_Comm  comm)
    {
    return Zoltan_Comm_Invert_Map( lengths_to, procs_to, nsends, self_msg,
        plengths_from, pprocs_from, pnrecvs, my_proc, nprocs, out_of_mem, 
        tag, comm);
    }
      
  static int Exchange_Sizes( int *sizes_to, int *procs_to, int  nsends, int  self_msg,
    int *sizes_from, int *procs_from, int  nrecvs, int *total_recv_size,
    int  my_proc, int  tag, MPI_Comm  comm)
    {
    return Zoltan_Comm_Exchange_Sizes(sizes_to, procs_to, nsends, self_msg,
        sizes_from, procs_from, nrecvs, total_recv_size, my_proc, tag, 
        comm);
    }

private:

  ZOLTAN_COMM_OBJ *Plan;
};

#endif
