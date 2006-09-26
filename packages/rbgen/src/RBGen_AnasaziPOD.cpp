
#include "RBGen_AnasaziPOD.h"
#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziBlockKrylovSchur.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace RBGen {
  
  AnasaziPOD::AnasaziPOD() :
    isInitialized_( false ),
    basis_size_( 16 ),
    comp_time_( 0.0 ) 
  {
  }

  void AnasaziPOD::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                               const Teuchos::RefCountPtr< Epetra_MultiVector >& ss )
  {
    if ( params->get("Basis Size", 16) < ss->NumVectors() ) {
      basis_size_ = params->get("Basis Size", 16);
    } 
    else { 
      basis_size_ = ss->NumVectors();
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
    bool inner_prod = true;
    //
    //  If the user is requesting more basis vectors than there are snapshots,
    //  compute the basis vectors using an outer product formulation.
    //
    if (basis_size_ > num_vecs) {
       step = 1;
       basis_size_ = num_vecs;
       inner_prod = false;
    }
    Epetra_Time timer( comm );
    int i, blockSize = 1;
    int nev = basis_size_;
    int maxBlocks = 2*basis_size_;
    double tol = 1e-14;
    string which="LM";
    int maxRestarts = 300;
    //
    // Create parameter list to pass into solver
    //
    Teuchos::ParameterList MyPL;
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Max Blocks", maxBlocks );
    MyPL.set( "Max Restarts", maxRestarts );
    MyPL.set( "Step Size", step );
    MyPL.set( "Tol", tol );
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
    if (inner_prod)
      ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec( localMap, blockSize ) );
    else
      ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec( ss_->Map(), blockSize ) );
    ivec->MvRandom();
    //
    // Call the constructor for the (A^T*A) operator
    Teuchos::RefCountPtr<Anasazi::EpetraSymMVOp> Amat = Teuchos::rcp( new Anasazi::EpetraSymMVOp(ss_, !inner_prod) );
    Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem =
      Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Amat, ivec) );
    
    // Inform the eigenproblem that the operator A is symmetric
    MyProblem->SetSymmetric( true ); 
    
    // Set the number of eigenvalues requested
    MyProblem->SetNEV( nev );
    
    // Inform the eigenproblem that you are finishing passing it information
    assert( MyProblem->SetProblem() == 0 );
    
    // Create a sorting manager to handle the sorting of eigenvalues in the solver
    Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
      Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );
    
    // Create an output manager to handle the I/O from the solver
    Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
      Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
    //MyOM->SetVerbosity( Anasazi::FinalSummary );
    
    // Initialize the Block Arnoldi solver
    Anasazi::BlockKrylovSchur<double,MV,OP> MyBlockKrylovSchur(MyProblem, MySort, MyOM, MyPL);
    
    timer.ResetStartTime();

    // Solve the problem to the specified tolerances or length
    MyBlockKrylovSchur.solve();

    comp_time_ = timer.ElapsedTime();
    //MyBlockKrylovSchur.currentStatus();

    // Obtain results directly
    Teuchos::RefCountPtr<std::vector<double> > evalr = MyProblem->GetEvals();
    
    // Retrieve eigenvectors
    Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> evecr = 
      Teuchos::rcp_dynamic_cast<Anasazi::EpetraMultiVec>( MyProblem->GetEvecs() );

    // Compute singular values/vectors and direct residuals.
    //
    // Compute singular values which are the square root of the eigenvalues
    //
    for (i=0; i<nev; i++) { sv_[i] = Teuchos::ScalarTraits<double>::squareroot( (*evalr)[i] ); }
    //
    // If we are using inner product formulation, 
    // then we must compute left singular vectors:
    //             u = Av/sigma
    //
    int info = 0;
    std::vector<double> tempnrm( nev );
    if (inner_prod) {
      basis_ = Teuchos::rcp( new Epetra_MultiVector(ss_->Map(), nev) );
      Epetra_MultiVector AV( ss_->Map(),nev );
      
      /* A*V */   
      info = AV.Multiply( 'N', 'N', 1.0, *ss_, *evecr, 0.0 );
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
      basis_ = Teuchos::rcp( new Epetra_MultiVector( Copy, *evecr, 0, nev ) );
      Epetra_MultiVector ATU( localMap, nev );
      Epetra_MultiVector V( localMap, nev );      
      
      /* A^T*U */   
      info = ATU.Multiply( 'T', 'N', 1.0, *ss_, *evecr, 0.0 );
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
  
} // end of RBGen namespace




