
#include "RBGen_LapackPOD.h"
#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace RBGen {
  
  LapackPOD::LapackPOD() :
    isInitialized_( false ),
    basis_size_( 16 ),
    comp_time_( 0.0 ) 
  {
  }

  void LapackPOD::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
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
    // Resize the singular value vector
    sv_.resize( ss->NumVectors() );

    // Set the snapshot set
    ss_ = ss;
  }

  void LapackPOD::computeBasis()
  {    
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm comm;
#endif    
    //
    // Variables for Epetra
    //
    Epetra_Time timer( comm );
    int num_vecs = ss_->NumVectors();
    int dim = ss_->GlobalLength();
    //
    // Check to see if there is more than one processor.
    // If not, just grab a copy of the multivector, else
    // export all information to processor 0 and compute basis.
    //
    if (comm.NumProc() > 1) {
      //
      // Create map putting all elements of vector on Processor 0.
      //
      Epetra_Map* Proc0Map;
      //
      if ( comm.MyPID() == 0 ) {
	Proc0Map = new Epetra_Map( dim, dim, 0, comm );
      } else {
	Proc0Map = new Epetra_Map( dim, 0, 0, comm );
      }
      Epetra_MultiVector Proc0MV( *Proc0Map, num_vecs );
      //
      // Create an exporter to get the global Epetra_MultiVector to a local Epetra_MultiVector.
      //
      Epetra_Export exporter( ss_->Map(), *Proc0Map );
      //
      // Export the Epetra_MultiVector
      //
      Proc0MV.Export(*ss_, exporter, Insert);
      //
      if ( comm.MyPID() == 0 ) {
	//
	// Now we can use the Proc0MV because it's on this processor.
	//
	// Create workspace for the SVD.
	//
	int info, lwork;
	double U[ 1 ], Vt[ 1 ];
	lwork = EPETRA_MAX( 3*num_vecs + dim, 5*num_vecs );
	std::vector<double> work( lwork );
	Epetra_LAPACK lapack;
	//
	// Compute the SVD.
	//
	timer.ResetStartTime();
	lapack.GESVD( 'O', 'N', dim, num_vecs, Proc0MV.Values(), Proc0MV.Stride(), &sv_[0], U, 1, Vt, 1, &work[0], &lwork, &info );
	comp_time_ = timer.ElapsedTime();      
	if (info != 0) { 
	  // THROW AN EXCEPTION HERE!
	  cout<< "The return value of the SVD is not 0!"<< endl;
	}
      }
      //
      // Communicate the singular values and vectors back to all processors.
      //
      comm.Broadcast( &sv_[0], basis_size_, 0 );
      //
      // Each processor needs to import the information back.
      //
      Epetra_Import importer( ss_->Map(), *Proc0Map );     
      ss_->Import(Proc0MV, importer, Insert);
      //
      // Create a view into the MultiVector for the basis.
      //
      basis_ = Teuchos::rcp( new Epetra_MultiVector( Copy, *ss_, 0, basis_size_ ) );
      //
      // Clean up
      //
      delete Proc0Map;
      //
    } else {
      //
      // Just use the multivector, it's on this processor.
      //
      // Create workspace for the SVD.
      //
      int info, lwork;
      double U[ 1 ], Vt[ 1 ];
      lwork = EPETRA_MAX( 3*num_vecs + dim, 5*num_vecs );
      std::vector<double> work( lwork );
      Epetra_LAPACK lapack;
      //
      // Compute the SVD
      //
      timer.ResetStartTime();
      lapack.GESVD( 'O', 'N', dim, num_vecs, ss_->Values(), ss_->Stride(), &sv_[0], U, 1, Vt, 1, &work[0], &lwork, &info );
      comp_time_ = timer.ElapsedTime();      
      if (info != 0) { 
	// THROW AN EXCEPTION HERE!
	cout<< "The return value of the SVD is not 0!"<< endl;
      }
      sv_.resize( basis_size_ );     
      basis_ = Teuchos::rcp( new Epetra_MultiVector( Copy, *ss_, 0, basis_size_ ) );
      //
    }
  }
  
} // end of RBGen namespace

