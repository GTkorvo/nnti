
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#include <string>
#include <sstream>
#include <iostream>

typedef int LO;
typedef int GO;
typedef double Scalar;

/////////////////////////////////////////////////////////////////////////
/* On a linux node, find the total memory currently allocated
 * to this process.
 *   Return the number of kilobytes allocated to this process.
 *   Return 0 if it is not possible to determine this.
 */
static long getProcessKilobytes()
{
long pageSize;

#ifdef _SC_PAGESIZE
  pageSize = sysconf(_SC_PAGESIZE);
#else
#warning "Page size query is not possible.  No per-process memory stats."
  return 0;
#endif

  pid_t pid = getpid();
  std::ostringstream fname;
  fname << "/proc/" << pid << "/statm";
  std::ifstream memFile;

  try{
    memFile.open(fname.str().c_str());
  }
  catch (...){
    return 0;
  }

  char buf[128];
  memset(buf, 0, 128);
  while (memFile.good()){
    memFile.getline(buf, 128);
    break;
  }

  memFile.close();

  std::istringstream sbuf(buf);
  long totalPages;
  sbuf >> totalPages;

  long pageKBytes = pageSize / 1024;
  totalPages = atol(buf);

  return totalPages * pageKBytes;
}

/////////////////////////////////////////////////////////////////////

int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);

  Teuchos::RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();
  int me = comm->getRank();
  int nprocs = comm->getSize();

  if (nprocs != 4)
      std::cout << "Run with 4 MPI ranks " << std::endl;

  typedef Tpetra::Map<LO, GO> map_t;
  GO numGlobalCoords = 4000000;
  LO numLocalCoords = 1000000;
  Teuchos::ParameterList myParams("testParameterList");
  myParams.set("memory_procs", "0");
  myParams.set("memory_output_stream", "std::cout");

  LO newnumLocalCoords = 1000000;
  if (me == 0)
      newnumLocalCoords = 999999;
  else if (me == 1)
      newnumLocalCoords = 1000001;
  else
      newnumLocalCoords = 1000000;

  typedef Tpetra::MultiVector<Scalar, LO, GO> mvector_t;

  long before = getProcessKilobytes();
  if (me == 0)
    std::cout << me << " "
              << getProcessKilobytes()
              << "   Before map construction " 
              << std::endl;

  for (int i = 0 ; i < 100; i++)
  {
      if (me == 0) 
        std::cout << me << " "
                  << getProcessKilobytes()
                  << "   Inside the loop " << i
                  << std::endl;
      Teuchos::RCP<const map_t> tmap = rcp(new map_t(numGlobalCoords, 
        numLocalCoords, 0, comm));
      Teuchos::RCP<const map_t> newTmap = rcp(new map_t(numGlobalCoords, 
        newnumLocalCoords, 0, comm));
      Teuchos::RCP<mvector_t> newMvector = rcp(new mvector_t(tmap, 3, true));
      Teuchos::RCP<Tpetra::Import<LO, GO> > importer = rcp(
        new Tpetra::Import<LO, GO>(tmap, newTmap));
      //defEnv->memory("Inside the loop after i = 0");
  }

  long after = getProcessKilobytes();
  if (me == 0)
    std::cout << me << " "
              << getProcessKilobytes()
              << "   After map construction "
              << std::endl;

  int iAmOK = (before == after);
  int weAreOK;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, iAmOK, &weAreOK);

  if (me == 0) {
    if (weAreOK) std::cout << "PASS" << std::endl;
    else         std::cout << "FAIL" << std::endl;
  }
}
