// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
// ***********************************************************************
//
// Testing Zoltan2::Environment
/*! \todo test timer
*/


#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>

using std::string;
using Teuchos::ParameterEntry;
using Teuchos::RCP;
using Teuchos::Comm;
using Zoltan2::Environment;

int checkErrorCode(Teuchos::RCP<const Teuchos::Comm<int> > &comm, int code)
{
  int rank = comm->getRank();
  if (code > 0) 
    std::cerr << "Proc " << rank << " error: " << code << std::endl;
  comm->barrier();

  // Will return 1 if any process has a non-zero code
  TEST_FAIL_AND_RETURN_VALUE(*comm, code==0, "test failure", 1);

  return 0;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0;

  ////////////////////////////////////////////////////////////
  // Create a default Environment and test it

  Environment *defEnv = NULL;

  try{
    defEnv = new Environment;
  }
  catch(std::exception &e){
    std::cerr << e.what() << std::endl;
    fail=1000;
  }

  if (checkErrorCode(comm, fail))
    return 1;

  if (!fail && defEnv->myRank_ != rank)
    fail = 1001;

  if (!fail && defEnv->numProcs_ != nprocs)
    fail = 1002;

  if (!fail && defEnv->comm_->getSize() != nprocs)
    fail = 1003;

  if (!fail && defEnv->doStatus() != true)
    fail = 1005;
  
  if (!fail && defEnv->doTiming() != false)
    fail = 1006;

  if (!fail && defEnv->doMemoryProfiling() != false)
    fail = 1007;

  if (!fail && defEnv->errorCheckLevel_ != Zoltan2::BASIC_ASSERTION)
    fail = 1008;

  if (!fail && defEnv->hasPartitioningParameters() != false)
    fail = 1009;

  if (!fail && defEnv->hasOrderingParameters() != false)
    fail = 1010;

  if (!fail && defEnv->hasColoringParameters() != false)
    fail = 1011;

  if (!fail && defEnv->hasMatchingParameters() != false)
    fail = 1012;

  if (checkErrorCode(comm, fail))
    return 1;

  delete defEnv;

  ////////////////////////////////////////////////////////////
  // Set a few parameters and create an Environment

  Teuchos::ParameterList myParams("testParameterList");

  myParams.set("debug_level", "detailed_status");        
  myParams.set("debug_procs", "all");   
  myParams.set("debug_output_stream", "std::cout");

  if (nprocs > 3)
    myParams.set("memory_procs", "0-1,3"); 
  else
    myParams.set("memory_procs", "0"); 

  myParams.set("memory_output_file", "memInfo.txt");

  myParams.set("speed_versus_quality", "speed");
  myParams.set("memory_versus_speed", "memory");

  Teuchos::ParameterList &parParams = myParams.sublist("partitioning");

  parParams.set("topology", "2,6,6");
  parParams.set("randomize_input", "true");
  parParams.set("objective", "minimize_cut_edge_weight");
  parParams.set("imbalance_tolerance", 1.2);

  Environment *env = NULL;

  try{
    env = new Environment(myParams, comm);
  }
  catch(std::exception &e){
    std::cerr << e.what() << std::endl;
    fail=2000;
  }

  if (!fail){
     try{
       env->debug(Zoltan2::BASIC_STATUS, "A basic debugging message.");
     }
     catch(std::exception &e){
       std::cerr << e.what() << std::endl;
       fail=3000;
     }
  }

  if (!fail){
     try{
       env->memory("Memory info");
       env->memory("Memory info next");
       env->memory("Memory info after");
     }
     catch(std::exception &e){
       std::cerr << e.what() << std::endl;
       fail=3002;
     }
  }

  if (checkErrorCode(comm, fail))
    return 1;

  if (!fail && env->myRank_ != rank)
    fail = 2001;

  if (!fail && env->numProcs_ != nprocs)
    fail = 2002;

  if (!fail && env->comm_->getSize() != nprocs)
    fail = 2003;

  if (!fail){
    const Teuchos::ParameterList &pl1 = env->getParameters();
    const ParameterEntry *dl = pl1.getEntryPtr("debug_level");
  
    if (!dl){
      fail = 2004;
    }
    else if (!(dl->isType<int>())){
      fail = 2013;
    }
    else{
      int value;
      int &val = dl->getValue<int>(&value);
      if (val != Zoltan2::DETAILED_STATUS) 
        fail = 2005;
    }
  }

  if (!fail && env->errorCheckLevel_ != Zoltan2::BASIC_ASSERTION)
    fail = 2008;

  if (!fail && env->hasPartitioningParameters() != true)
    fail = 2009;

  if (!fail && env->hasOrderingParameters() != false)
    fail = 2010;

  if (!fail && env->hasColoringParameters() != false)
    fail = 2011;

  if (!fail && env->hasMatchingParameters() != false)
    fail = 2012;

  if (checkErrorCode(comm, fail))
    return 1;

  if (rank==0){
    std::cout << "\nA test parameter list" << std::endl;
    const Teuchos::ParameterList &envParams = env->getParameters();
    try{
      envParams.print();
    }
    catch(std::exception &e){
      std::cerr << e.what() << std::endl;
      fail=2013;
    }
  }

  if (checkErrorCode(comm, fail))
    return 1;

  ////////////////////////////////////////////////////////////
  // Given an existing Environment, get its parameters and
  // add some new parameters and create a new Environment.

  RCP<const Comm<int> > oldComm = env->comm_;
  const Teuchos::ParameterList &oldParams = env->getUnvalidatedParameters();
  
  Teuchos::ParameterList newParams = oldParams;
  newParams.set("error_check_level", "debug_mode_assertions");
  newParams.set("memory_versus_speed", "speed");
  newParams.remove("memory_output_file");
  
  Teuchos::ParameterList &newPartParams = newParams.sublist("partitioning");
  newPartParams.set("imbalance_tolerance", "1.05");
  newPartParams.set("algorithm", "phg");
  newPartParams.set("objective", "minimize_cut_edge_weight");
  
  RCP<Environment> newEnv;

  try{
    newEnv = Teuchos::rcp(new Environment(newParams, oldComm));
  }
  catch(std::exception &e){
    std::cerr << e.what() << std::endl;
    fail=3000;
  }

  if (checkErrorCode(comm, fail))
    return 1;

  if (!fail && newEnv->errorCheckLevel_ != Zoltan2::DEBUG_MODE_ASSERTION)
    fail = 3001;

  if (!fail && !newEnv->hasPartitioningParameters())
    fail = 3002;

  if (!fail && rank==0){
    std::cout << "\nA few changes/additions to the list" << std::endl;
    const Teuchos::ParameterList &envParams = newEnv->getParameters();
    try{
      envParams.print();
    }
    catch(std::exception &e){
      std::cerr << e.what() << std::endl;
      fail=3003;
    }
  }

  if (checkErrorCode(comm, fail))
    return 1;

  delete env;

  if (rank==0)
    std::cout << "PASS" << std::endl;

  return 0;
}
