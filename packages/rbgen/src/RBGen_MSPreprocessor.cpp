

#include "RBGen_MSPreprocessor.h"
#include "Epetra_MultiVector.h"

namespace RBGen {

  MSPreprocessor::MSPreprocessor(): isInitialized_(false), steady_scale_(1.0)
  {
  }

  void MSPreprocessor::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params, 
                                   const Teuchos::RefCountPtr< FileIOHandler<Epetra_MultiVector> >& fileio )
  {
    // Save the file i/o handler.
    fileio_ = fileio;

    // Get the preprocessing method sublist.
    const Teuchos::ParameterList& preproc_params = params->sublist( "Preprocessing Method" );

    //
    // Try to get the steady state file from the parameter list
    //
    try {
      steady_file_ = Teuchos::getParameter< std::string >( preproc_params, "Steady State File" );
    }
    catch (std::exception &e) {
      cout<<"The steady state file has not been specified"<<endl;
    }

    //
    // Check if there is a scaling factor on the steady state vector
    //
    steady_scale_ = 1.0;
    if ( preproc_params.isParameter("Steady State Scaling Factor") ) {
      steady_scale_ = Teuchos::getParameter< double >( preproc_params, "Steady State Scaling Factor" );
    }

/*    try {
      scalings_ = Teuchos::getParameter< std::vector< double > >( *params, "Snapshot Scaling" );
    }
    catch (std::exception &e) {
      cout<<"The snapshot scaling vector has not been specified in the input!"<<endl;
    }
    //
    // Try to get the snapshot scaling indices vector
    //
    try {
      scaling_idx_ = Teuchos::getParameter< std::vector< std::pair<int,int> > >( *params, "Snapshot Scaling Indices" );
    }
    catch (std::exception &e) {
      cout<<"The snapshot scaling indices have not been specified in the input!"<<endl;
    }      
*/
    if ( scalings_.size() != scaling_idx_.size() ) {
      cout<<"The scaling vector is not the same size as the number of index pairs!"<< endl;
    }
    //
    // We're initialized now.
    // 
    isInitialized_ = true;
  }

  void MSPreprocessor::Preprocess( Teuchos::RefCountPtr<Epetra_MultiVector>& ss )
  {
    if (isInitialized_) {      
      //
      // Get the steady state vector from the file.
      //
      std::vector< std::string > steady_file( 1, steady_file_ );
      Teuchos::RefCountPtr<Epetra_MultiVector> steadyMV = fileio_->Read( steady_file );
      Epetra_MultiVector* colMV = 0;
      // 
      // Remove the scaled steady state vector.
      //
      for (int i=0; i<ss->NumVectors(); i++) {
	colMV = new Epetra_MultiVector( View, *ss, i, 1 );
	colMV->Update( -1.0*steady_scale_, *steadyMV, 1.0 );
	delete colMV;
      }
    }
    else {
      // Throw error here!
    }
  }
  
} // end of RBGen namespace

