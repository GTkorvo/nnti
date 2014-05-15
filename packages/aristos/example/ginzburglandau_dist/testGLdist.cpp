//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

//#include "Epetra_config.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
//#include "mpi.h"
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

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

int main(int argc, char *argv[])
{

  // This is a standard communicator declaration.
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

//#ifdef HAVE_MPI
//  MPI_Init(&argc,&argv);
//  Epetra_MpiComm Comm(MPI_COMM_WORLD);
//#else
//  Epetra_SerialComm Comm;
//#endif

  Teuchos::CommandLineProcessor clp;

  double beta = 1.0;

  std::string geomfile;
  clp.setOption("geomfile", &geomfile, "The geometry file base name.");

  clp.setOption("use-stratimikos", "use-aztecoo",
    &GLdistApp::GLdistYUEpetraDataPool::useStratimikos,
    "Use Stratimikos or AztecOO");

  clp.setOption("stratimikos-params-file",
    &GLdistApp::GLdistYUEpetraDataPool::stratimikosXmlFile,
    "Stratimikos input file.");
  
  // Want derivative check?
  bool derchk = false;

  bool wantstats = true;   // choose true if output of solver and timing info in file stats.txt is desired

  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn= clp.parse(argc, argv);
  if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
    return 0;
  }
  if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
    return 1; // Error!
  }

  TEUCHOS_ASSERT(geomfile.length());

  int mypid = Comm.MyPID();

  ofstream outfile("stats.txt", ios_base::app);

  // Now we can build the DataPool object ...
  GLdistApp::GLdistYUEpetraDataPool dat ( &Comm, beta, const_cast<char*>(geomfile.c_str()) );

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
  parlist.set("Gradient of Lagrangian Tolerance", 1e-6);
  parlist.set("Constraints Tolerance", 1e-6);

  int iter, iflag;
  exy->PutScalar(1.0);
  exu->PutScalar(1.0);

  if ((mypid==0) && wantstats)
    outfile << "Geometry: " << geomfile << " , Subdomains: " << Teuchos::GlobalMPISession::getNProc() << "   *************************" << endl;
  Epetra_Time timer(Comm);
  sqp.run(*x, *c, *l, iter, iflag, parlist);
  if ((mypid==0) && wantstats)
    outfile << endl << "Elapsed Time: " << timer.ElapsedTime() << "   *************************" << "\n\n";

  outfile.close();

  std::cout << "End Result: TEST PASSED\n";

  return 0;

  //dat.PrintSolutionVTK(exy);

}
