#include "RBGen_EpetraMVPreprocessorFactory.h"

#include "Epetra_MultiVector.h"

namespace RBGen {

  EpetraMVPreprocessorFactory::EpetraMVPreprocessorFactory()
  {
    preproc_types.push_back("none");
    preproc_types.push_back("ModifiedSS");
  }

  Teuchos::RefCountPtr<Preprocessor<Epetra_MultiVector> >
  EpetraMVPreprocessorFactory::create( const Teuchos::ParameterList& params )
  {
    // See if the "Preprocessing" sublist exists
    if (!params.isSublist( "Preprocessing Method" ) ) {
      // TO DO:  THROW EXCEPTION!!!
    }

    // Get the preprocessing method sublist.
    const Teuchos::ParameterList& preproc_params = params.sublist( "Preprocessing Method" );

    // Get the preprocessing parameter.
    std::string method = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(preproc_params),
                                                             "Method" );

    // Create the preprocessor.
    Teuchos::RefCountPtr<Preprocessor< Epetra_MultiVector > > RBPreprocessor;

    // Basic preprocessor does nothing
    if ( method == "none" ) {
      RBPreprocessor = Teuchos::rcp( new NoPreprocessor<Epetra_MultiVector>() );
    } else
    // Modified snapshot preprocessor
    if ( method == "Modified Snapshot" ) {
      RBPreprocessor = Teuchos::rcp( new MSPreprocessor() );
    } else
    {
       // Throw an exception because the preprocessing method was not recognized by this factory.
    }
    //
    // Return the preprocessor created
    //
    return RBPreprocessor;

  }

} // end of RBGen namespace

