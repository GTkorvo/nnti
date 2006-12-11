//#include "Epetra_config.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"

#include "Aristos_YUEpetraVector.hpp"
#include "Aristos_SQPAlgo.hpp"

#include "GLdistApp_SchurOp.hpp"
#include "GLdistApp_GLdistYUEpetraObjective.hpp"
#include "GLdistApp_GLdistYUEpetraHessVec.hpp"
#include "GLdistApp_GLdistYUEpetraLagMult.hpp"
#include "GLdistApp_GLdistYUEpetraFeasStep.hpp"
#include "GLdistApp_GLdistYUEpetraConstraints.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"


int main(int argc, char *argv[])
{

  // This is a standard communicator declaration.
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  double beta = 1.0;
  
  // Want derivative check?
  bool derchk = false;

  bool wantstats = false;   // choose true if output of solver and timing info in file stats.txt is desired

  int mypid = Comm.MyPID();

  ofstream outfile("stats.txt", ios_base::app);

  // Now we can build the DataPool object ...
  GLdistApp::GLdistYUEpetraDataPool dat ( &Comm, beta, argv[1] );
  
  Epetra_Map statemap((dat.getA())->DomainMap());
  Epetra_Map controlmap((dat.getB())->DomainMap());

  Teuchos::RefCountPtr<Epetra_Vector> exy = Teuchos::rcp( new Epetra_Vector(statemap) );
  Teuchos::RefCountPtr<Epetra_Vector> exu = Teuchos::rcp( new Epetra_Vector(controlmap) );
  Teuchos::RefCountPtr<Epetra_Vector> ely = Teuchos::rcp( new Epetra_Vector(statemap) );
  Teuchos::RefCountPtr<Epetra_Vector> elu = Teuchos::null;
  Teuchos::RefCountPtr<Epetra_Vector> ecy = Teuchos::rcp( new Epetra_Vector(statemap) );
  Teuchos::RefCountPtr<Epetra_Vector> ecu = Teuchos::null;
  Teuchos::RefCountPtr<Epetra_Vector> ediry = Teuchos::rcp( new Epetra_Vector(statemap) );
  Teuchos::RefCountPtr<Epetra_Vector> ediru = Teuchos::rcp( new Epetra_Vector(controlmap) );

  Teuchos::RefCountPtr<Aristos::Vector> x   = Teuchos::rcp( new Aristos::YUEpetraVector( exy, exu ) );
  Teuchos::RefCountPtr<Aristos::Vector> c   = Teuchos::rcp( new Aristos::YUEpetraVector( ecy, ecu ) );
  Teuchos::RefCountPtr<Aristos::Vector> l   = Teuchos::rcp( new Aristos::YUEpetraVector( ely, elu ) );
  Teuchos::RefCountPtr<Aristos::Vector> dir = Teuchos::rcp( new Aristos::YUEpetraVector( ediry, ediru ) );

  Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool> rcpdat = Teuchos::rcp( &dat, false );
  
  // construct Objective object
  GLdistApp::GLdistYUEpetraObjective obj(rcpdat);

  // construct Constraints object
  GLdistApp::GLdistYUEpetraConstraints constr(rcpdat);

  // construct Feasibility Step object
  GLdistApp::GLdistYUEpetraFeasStep feasstep(rcpdat);

  // construct Lagrange Multiplier object
  GLdistApp::GLdistYUEpetraLagMult lagmult(rcpdat);

  // construct HessVec object
  GLdistApp::GLdistYUEpetraHessVec hessvec(rcpdat);

  // construct SQP algorithm
  Aristos::SQPAlgo sqp(dat, obj, constr, hessvec, lagmult, feasstep);

    
  // Derivative check.
  exy->Random();
  exu->Random();
  ely->Random();
  ediry->Random();
  ediru->Random();
  if (derchk)
    sqp.runDerivativeCheck(*x, *l, *dir);

  // SQP run.
  Teuchos::ParameterList parlist;
  parlist.set("Max Number of SQP Iterations", 50);

  int iter, iflag;
  exy->PutScalar(10.0);
  exu->PutScalar(10.0);

  if ((mypid==0) && wantstats)
    outfile << "Geometry: " << argv[1] << " , Subdomains: " << argv[2] << "   *************************" << endl;
  Epetra_Time timer(Comm);
  sqp.run(*x, *c, *l, iter, iflag, parlist);
  if ((mypid==0) && wantstats)
    outfile << endl << "Elapsed Time: " << timer.ElapsedTime() << "   *************************" << "\n\n";

  outfile.close();
  
  //dat.PrintVec(exy);

}
