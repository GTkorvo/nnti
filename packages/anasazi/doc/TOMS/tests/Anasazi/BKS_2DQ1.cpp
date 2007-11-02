#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"


using namespace Teuchos;
using namespace Anasazi;

Teuchos::RefCountPtr<Epetra_CrsMatrix> createMatrix(Epetra_SerialComm &Comm, const int nx);

int main(int argc, char *argv[])
{

  Epetra_SerialComm Comm;

  const int blockSize=1;
  const std::string which("LM");
  bool testFailed;
  bool debug = false;
  bool verbose = false;
  bool herm = true;

  int nx = 50;
  int nev = 10;
  int ncv = 100;
  int maxRestarts = 8;
  double tol = 1.0e-10;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Debugging details.");
  cmdp.setOption("herm","nonherm",&herm,"Specify Hermitian?");
  cmdp.setOption("nx",&nx,"Number of elements in the Laplacian matrix.");
  cmdp.setOption("nev",&nev,"Number of desired eigenvalues.");
  cmdp.setOption("ncv",&ncv,"Number of columns in the Krylov basis.");
  cmdp.setOption("maxRestarts",&maxRestarts,"Maximum number of restarts.");
  cmdp.setOption("tol",&tol,"Tolerance of desired eigenvalues.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;

  // Create default output manager 
  RefCountPtr<OutputManager<double> > MyOM = rcp( new BasicOutputManager<double>() );

  // local output manager
  MyOM->setVerbosity( Warnings );
  // Set verbosity level for solver
  int verbosity = Errors | FinalSummary | TimingDetails;
  if (debug) {
    verbosity |= IterationDetails | Warnings | FinalSummary | TimingDetails;
  }
  if (verbose) {
    // solver verb level
    verbosity |= Warnings | FinalSummary | TimingDetails | IterationDetails;
  }
  if (debug) {
    verbosity |= Debug;
  }

  MyOM->stream(Warnings) 
       << Anasazi_Version() << endl
       << "SERIAL" << endl
       << "BlockKrylovSchur, Epetra" << endl
       << "Num elements: " << nx << endl << endl;

  typedef double ScalarType;
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef ScalarTraits<ScalarType>                   SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef MultiVecTraits<ScalarType,MV>     MVT;
  typedef OperatorTraits<ScalarType,MV,OP>  OPT;
  const ScalarType ONE  = SCT::one();


  // 
  // Create the operator
  RefCountPtr<Epetra_CrsMatrix> A = createMatrix(Comm,nx);
  //
  // Create initial vector
  RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(A->OperatorDomainMap(), blockSize) );
  // ivec->Random();
  // set to e_1
  ivec->Scale(0.0);
  (*ivec)[0][0] = 1.0;
  cout << "ivec.dim: " << MVT::GetVecLength(*ivec) << endl;

  // Create eigenproblem
  RefCountPtr<BasicEigenproblem<ScalarType,MV,OP> > MyProblem =
    rcp( new BasicEigenproblem<ScalarType, MV, OP>(A, ivec) );
  //
  // Inform the eigenproblem of Hermitian structure
  MyProblem->setHermitian(herm);
  //   
  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->setNEV( nev );
  // 
  // Inform the eigenproblem that you are finishing passing it information
  if ( MyProblem->setProblem() == false ) {
    MyOM->print(Errors,"BasicEigenproblem::setProblem() failed.");
    return -1;
  }

  //
  // Create parameter list to pass into solver manager
  //
  ParameterList MyPL;
  MyPL.set( "Which", which );
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Num Blocks", ncv );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Orthogonalization", "DGKS" );
  // set this to match ARPACK: 1/.717^2
  MyPL.set( "Orthogonalization Constant", 1.0/(.717*.717) );
  //
  // Create the solver manager
  BlockKrylovSchurSolMgr<ScalarType,MV,OP> MySolverMan(MyProblem, MyPL);
  
  // Solve the problem to the specified tolerances or length
  ReturnType returnCode = MySolverMan.solve();
  testFailed = false;
  if (returnCode != Converged) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  RefCountPtr<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    ostringstream os;
    os.setf(ios::scientific, ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<ScalarType> normV( numev );
    SerialDenseMatrix<int,ScalarType> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = sol.Evals[i].realpart;
    }
    RefCountPtr<MV> Avecs = MVT::Clone( *evecs, numev );
    OPT::Apply( *A, *evecs, *Avecs );
    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Avecs );
    // compute norm of residuals
    MVT::MvDot( *Avecs, *Avecs, &normV );

    os << "Direct residual norms computed in BlockDavidson2DQ1_e.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(2)" << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) / sol.Evals[i].realpart );
      }
      else {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) );
      }
      os << setw(20) << sol.Evals[i].realpart << setw(20) << normV[i] << endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    MyOM->stream(Warnings) << endl << os.str() << endl;
  }
  else {
    MyOM->stream(Warnings) << "solve() return no eigenpairs." << endl;
  }

  int ret;
  string msg;
  if (testFailed) {
    msg = "End Result: TEST FAILED\n";
    ret = -1;
  }
  else {
    msg = "End Result: TEST PASSED\n";
    ret = 0;
  }

  MyOM->print(Warnings,msg);
  return ret;
}



Teuchos::RefCountPtr<Epetra_CrsMatrix> createMatrix(Epetra_SerialComm &Comm, const int nx) {
  //  Dimension of the matrix
  int NumGlobalElements = nx*nx;  // Size of matrix nx*nx

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(NumGlobalElements, 0, Comm);

  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor
  std::vector<int> NumNz(NumMyElements);

  /* We are building a matrix of block structure:
  
      | T -I          |
      |-I  T -I       |
      |   -I  T       |
      |        ...  -I|
      |           -I T|

   where each block is dimension nx by nx and the matrix is on the order of
   nx*nx.  The block T is a tridiagonal matrix. 
  */

  for (int i=0; i<NumMyElements; i++) {
    if (MyGlobalElements[i] == 0 || MyGlobalElements[i] == NumGlobalElements-1 || 
        MyGlobalElements[i] == nx-1 || MyGlobalElements[i] == nx*(nx-1) ) {
      NumNz[i] = 3;
    }
    else if (MyGlobalElements[i] < nx || MyGlobalElements[i] > nx*(nx-1) || 
             MyGlobalElements[i]%nx == 0 || (MyGlobalElements[i]+1)%nx == 0) {
      NumNz[i] = 4;
    }
    else {
      NumNz[i] = 5;
    }
  }

  // Create an Epetra_Matrix

  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(::Copy, Map, &NumNz[0]) );

  // Compute coefficients for discrete convection-diffution operator
  const double one = 1.0;
  std::vector<double> Values(4);
  std::vector<int> Indices(4);
  double h = one /(nx+1);
  double h2 = h*h;
  Values[0] = -one/h2; Values[1] = -one/h2; Values[2] = -one/h2; Values[3]= -one/h2;
  // double diag = 4.0 / h2;
  // RBL: added diagonal element to remove multiplicities to discourage need for
  // re-ortho. diag is different for each diagonal element: see below
  int NumEntries, info;

  for (int i=0; i<NumMyElements; i++)
  {
    if (MyGlobalElements[i]==0)
    {
      Indices[0] = 1;
      Indices[1] = nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == nx*(nx-1))
    {
      Indices[0] = nx*(nx-1)+1;
      Indices[1] = nx*(nx-2);
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == nx-1)
    {
      Indices[0] = nx-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
      Indices[0] = 2*nx-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
    {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
      Indices[0] = nx*(nx-1)-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] < nx)
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] > nx*(nx-1))
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i]%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]+1;
      Indices[1] = MyGlobalElements[i]-nx;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if ((MyGlobalElements[i]+1)%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]-nx;
      Indices[1] = MyGlobalElements[i]+nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
      Indices[0] = MyGlobalElements[i]-1;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      Indices[3] = MyGlobalElements[i]+nx;
      NumEntries = 4;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    // Put in the diagonal entry
    // diag from before: double diag = 4.0 / h2;
    // cout << MyGlobalElements[i] << '\t' << (MyGlobalElements[i]%nx)+1 << endl;
    double diag = 4.0 / h2 + ((MyGlobalElements[i]%nx)+1)*h;
    info = A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i]);
    assert( info==0 );
  }

  // Finish up
  info = A->FillComplete();
  assert( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  info = A->OptimizeStorage();
  assert( info==0 );
  return A;
}
