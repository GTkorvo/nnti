
#include "RBGen_AnasaziPOD.h"
#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace RBGen {
  
  AnasaziPOD::AnasaziPOD() :
    isInitialized_( false ),
    isInner_( true ),
    basis_size_( 16 ),
    comp_time_( 0.0 ) 
  {
  }

  void AnasaziPOD::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                               const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
			       const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio )
  {

   // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList& rbmethod_params = params->sublist( "Reduced Basis Method" );

    if ( rbmethod_params.get("Basis Size", 16) < ss->NumVectors() ) {
      basis_size_ = rbmethod_params.get("Basis Size", 16);
    } 
    else { 
      basis_size_ = ss->NumVectors();
    }
    // Get the inner / outer product form of the operator
    isInner_ = ( rbmethod_params.get("Anasazi POD Operator Form","Inner")=="Inner"? true : false );

    // See if there is a matrix to be used for an inner-product in the orthogonal basis construction
    if (rbmethod_params.isParameter( "Inner Product Weighting Matrix" )) {
      string matFile = Teuchos::getParameter<string>( rbmethod_params, "Inner Product Weighting Matrix" );
      std::vector<std::string> filename(1,matFile);
      op_ = fileio->Read( filename );
    }
    
    // Resize the singular value vector 
    sv_.resize( basis_size_ );    

    // Set the snapshot set
    ss_ = ss;
  }

  void AnasaziPOD::computeBasis()
  {    
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm comm;
#endif    

    int MyPID = comm.MyPID();
    //
    //  Variables used for the Block Krylov Method
    //
    int step = 5;
    int num_vecs = ss_->NumVectors();
    //
    //  If the user is requesting more basis vectors than there are snapshots,
    //  compute the basis vectors using an outer product formulation.
    //
    if (basis_size_ > num_vecs) {
       step = 1;
       basis_size_ = num_vecs;
       isInner_ = false;
    }
    Epetra_Time timer( comm );
    int i, blockSize = 1;
    int nev = basis_size_;
    int maxBlocks = 2*basis_size_;
    int maxRestarts = 300;
    int verbosity = Anasazi::Warnings + Anasazi::Errors;
    double tol = 1e-14;
    string which="LM";
    //
    // Create parameter list to pass into solver
    //
    Teuchos::ParameterList MyPL;
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Num Blocks", maxBlocks );
    MyPL.set( "Maximum Restarts", maxRestarts );
    MyPL.set( "Which", which );	
    MyPL.set( "Step Size", step );
    MyPL.set( "Convergence Tolerance", tol );
    //
    //  Typedefs for Anasazi solvers
    //
    typedef Anasazi::MultiVec<double> MV;
    typedef Anasazi::Operator<double> OP;
    //
    // Create a map for the columns of the snapshot set.
    //
    Epetra_LocalMap localMap(num_vecs, 0, comm);
    //
    // Create the initial vector and randomize it.
    //
    Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> ivec;
    if (isInner_)
      ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec( localMap, blockSize ) );
    else
      ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec( ss_->Map(), blockSize ) );
    ivec->MvRandom();
    //
    Teuchos::RefCountPtr<OP> Amat;

    if (op_ != Teuchos::null) {
      // Call the constructor for the (WA)^T*WA operator
      Amat = Teuchos::rcp( new Anasazi::EpetraWSymMVOp( ss_, op_ ) );
    }
    else {
      // Call the constructor for the (A^T*A) operator
      Amat = Teuchos::rcp( new Anasazi::EpetraSymMVOp(ss_, !isInner_) );
    }
 
    // Create the eigenproblem
    Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem =
      Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Amat, ivec) );
  
    // Inform the eigenproblem that the operator A is symmetric
    MyProblem->setHermitian( true ); 
    
    // Set the number of eigenvalues requested
    MyProblem->setNEV( nev );
    
    // Inform the eigenproblem that you are finishing passing it information
    bool boolret = MyProblem->setProblem();
    if (boolret != true) {
      if (MyPID == 0) {
	cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
      }
    }
    
    // Initialize the Block Arnoldi solver
    Anasazi::BlockKrylovSchurSolMgr<double,MV,OP> MySolverMgr(MyProblem, MyPL);
    
    timer.ResetStartTime();

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr.solve();
    if (returnCode != Anasazi::Converged && MyPID==0) {
      cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
    }

    comp_time_ = timer.ElapsedTime();

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    int numev = sol.numVecs;
    
    if (numev > 0) {

      // Retrieve eigenvectors
      std::vector<int> index(numev);
      for (i=0; i<numev; i++) { index[i] = i; }
      Anasazi::EpetraMultiVec* evecs = dynamic_cast<Anasazi::EpetraMultiVec* >(sol.Evecs->CloneView( index )); 
      //
      // Compute singular values which are the square root of the eigenvalues
      //
      for (i=0; i<nev; i++) { 
	sv_[i] = Teuchos::ScalarTraits<double>::squareroot( Teuchos::ScalarTraits<double>::magnitude(evals[i].realpart) ); 
      }
      //
      // If we are using inner product formulation, 
      // then we must compute left singular vectors:
      //             u = Av/sigma
      //
      int info = 0;
      std::vector<double> tempnrm( nev );
      if (isInner_) {
	basis_ = Teuchos::rcp( new Epetra_MultiVector(ss_->Map(), nev) );
	Epetra_MultiVector AV( ss_->Map(),nev );
	
	/* A*V */   
	info = AV.Multiply( 'N', 'N', 1.0, *ss_, *evecs, 0.0 );
	AV.Norm2( &tempnrm[0] );
	
	/* U = A*V(i)/S(i) */
	Epetra_LocalMap localMap2( nev, 0, ss_->Map().Comm() );
	Epetra_MultiVector S( localMap2, nev );
	for( i=0; i<nev; i++ ) { S[i][i] = 1.0/tempnrm[i]; }
	info = basis_->Multiply( 'N', 'N', 1.0, AV, S, 0.0 );
	
	/* Compute direct residuals : || Av - sigma*u ||
	   for (i=0; i<nev; i++) { S[i][i] = sv_[i]; }
	   info = AV.Multiply( 'N', 'N', -1.0, *basis_, S, 1.0 );
	   AV.Norm2( &tempnrm[0] );
	   if (MyPID == 0) {
	   cout<<"Singular Value"<<"\t\t"<<"Direct Residual"<<endl;
	   cout<<"------------------------------------------------------"<<endl;
	   for (i=0; i<nev; i++) {
	   cout<< sv_[i] << "\t\t\t" << tempnrm[i] << endl;
	   }
	   cout<<"------------------------------------------------------"<<endl;
	   }
	*/
      } else {
	basis_ = Teuchos::rcp( new Epetra_MultiVector( Copy, *evecs, 0, nev ) );
	Epetra_MultiVector ATU( localMap, nev );
	Epetra_MultiVector V( localMap, nev );      
	
	/* A^T*U */   
	info = ATU.Multiply( 'T', 'N', 1.0, *ss_, *evecs, 0.0 );
	ATU.Norm2( &tempnrm[0] );
	
	/* V = A^T*U(i)/S(i) */
	Epetra_LocalMap localMap2( nev, 0, ss_->Map().Comm() );
	Epetra_MultiVector S( localMap2, nev );
	for( i=0; i<nev; i++ ) { S[i][i] = 1.0/tempnrm[i]; }
	info = V.Multiply( 'N', 'N', 1.0, ATU, S, 0.0 );
	
	/* Compute direct residuals : || (A^T)u - sigma*v ||     
	   for (i=0; i<nev; i++) { S[i][i] = sv_[i]; }
	   info = ATU.Multiply( 'N', 'N', -1.0, V, S, 1.0 );
	   ATU.Norm2( tempnrm );
	   if (comm.MyPID() == 0) {
	   cout<<"Singular Value"<<"\t\t"<<"Direct Residual"<<endl;
	   cout<<"------------------------------------------------------"<<endl;
	   for (i=0; i<nev; i++) {
	   cout<< sv_[i] << "\t\t\t" << tempnrm[i] << endl;
	   }
	   cout<<"------------------------------------------------------"<<endl;
	   }
	*/
      }
    }
  }
} // end of RBGen namespace




