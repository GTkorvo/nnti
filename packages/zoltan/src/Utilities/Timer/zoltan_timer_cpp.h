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
// C++ wrappers for Zoltan Timer library.  
//
// Initialization:
//
//   C++ style:  Zoltan_Timer_Object ztimer(timer_flag);
//                                   ztimer();
//
// ************************************************************************

#ifndef ZOLTAN_TIMER_CPP_H_
#define ZOLTAN_TIMER_CPP_H_

#include "zoltan_timer.h"
#include <string>

class Zoltan_Timer_Object {

public:

  Zoltan_Timer_Object(int flag = ZOLTAN_TIME_WALL) {
    // Assumption: MPI has been initialized prior to this call.
    ZTStruct = Zoltan_Timer_Create(flag);
  }

  ~Zoltan_Timer_Object() {
    Zoltan_Timer_Destroy(&ZTStruct);
  }

  int Init(const int proc, const int use_barrier, const std::string & name) {
    return Zoltan_Timer_Init(this->ZTStruct, proc, use_barrier,
                             const_cast<char*>(name.c_str()));
  }

  int Reset(const int ts_idx, const int proc, const int use_barrier, 
            const std::string & name) {
    return Zoltan_Timer_Reset(this->ZTStruct, ts_idx, proc, use_barrier,
                              const_cast<char*>(name.c_str()));
  }

  int Start(const int ts_idx, const int proc, const MPI_Comm comm) {
    return Zoltan_Timer_Start(this->ZTStruct, ts_idx, proc, comm);
  }

  int Stop(const int ts_idx, const int proc, const MPI_Comm comm) {
    return Zoltan_Timer_Stop(this->ZTStruct, ts_idx, proc, comm);
  }

  int Print(const int ts_idx, const int proc) {
    return Zoltan_Timer_Print(this->ZTStruct, ts_idx, proc);
  }

  int PrintAll(const int proc) {
    return Zoltan_Timer_PrintAll(this->ZTStruct, proc);
  }

private:

  struct Zoltan_Timer *ZTStruct;
};

#endif
