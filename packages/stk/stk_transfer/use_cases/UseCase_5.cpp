/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */      
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */      
/*  United States Government.                                             */      
/*------------------------------------------------------------------------*/

#include <Intrepid_FieldContainer.hpp>


#include <stk_mesh/base/Comm.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include <stk_transfer/Transfer.hpp>


int main(int argc, char **argv)
{
  stk::diag::Timer timer("Transfer Use Case 5", 
                          use_case::TIMER_TRANSFER, 
                          use_case::timer());
  stk::diag::Timer timer_point_to_point(" Point To Point", timer);
  use_case::timerSet().setEnabledTimerMask(use_case::TIMER_ALL);

  bool status = true;

  const double TOLERANCE = 0.000001;
  const double  rand_max = RAND_MAX;
  enum {           DIM = 3  };
  enum { FROMNUMPOINTS = 100  };
  enum {   TONUMPOINTS = 100  };

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  stk::ParallelMachine comm = use_case_environment.m_comm;
 
  typedef stk::transfer::Transfer::MDArray MDArray;

  stk::transfer::Transfer transfer;

  MDArray FromPoints (FROMNUMPOINTS,DIM), 
          ToPoints   (  TONUMPOINTS,DIM), 
          FromValues (FROMNUMPOINTS,  2), 
          ToValues   (  TONUMPOINTS,  2);
  for (unsigned i=0 ; i<FROMNUMPOINTS; ++i) {
    double l=0, q=0; 
    for (unsigned j=0 ; j<DIM; ++j) {
      FromPoints(i,j) = rand()/rand_max;
      l +=   FromPoints(i,j);
      q += j*FromPoints(i,j);
    }
    FromValues(i,0) = l;
    FromValues(i,1) = q;
  }
  for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
    for (unsigned j=0 ; j<DIM; ++j) {
      ToPoints(i,j) = rand()/rand_max;
    }
  }

  {
    stk::diag::TimeBlock __timer_point_to_point(timer_point_to_point);
    try {
      transfer.PointToPoint(ToValues, ToPoints, FromValues, FromPoints, comm);
    } catch (...) {
      status = status && false;
    } 
  }

  bool success = true;
  for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
    double check_l = 0;
    double check_q = 0;
    for (unsigned j=0 ; j<DIM; ++j) {
      check_l +=   ToPoints(i,j);
      check_q += j*ToPoints(i,j);
    }
    if (TOLERANCE < fabs(check_l-ToValues(i,0)) ) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" ToValues:"<<ToValues(i,0)
                <<" check:"<<check_l
                <<" error:"<<fabs(check_l-ToValues(i,0))
                <<std::endl;
      success = false;
    }
    if (TOLERANCE < fabs(check_q-ToValues(i,1)) ) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" ToValues:"<<ToValues(i,1)
                <<" check:"<<check_q
                <<" error:"<<fabs(check_q-ToValues(i,1))
                <<std::endl;
      success = false;
    }
  }
  status = status && success;
  timer.stop();
//stk::diag::printTimersTable(std::cout, timer, 
//      stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME, false, comm);


  const bool collective_result = use_case::print_status(comm, status);
  const int return_code = collective_result ? 0 : -1;
  return return_code;
}
