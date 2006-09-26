#include "RBGen_EpetraMVPreprocessorFactory.h"

#include "Epetra_MultiVector.h"
#include <vector>
#include <string>

namespace RBGen {

  EpetraMVPreprocessorFactory::EpetraMVPreprocessorFactory()
  {
    preproc_types.push_back("none");
    preproc_types.push_back("ModifiedSS");
  }

  Teuchos::RefCountPtr<Preprocessor<Epetra_MultiVector> >
  EpetraMVPreprocessorFactory::create( const Teuchos::ParameterList& params )
  {

    if (!params.isParameter( "Preprocessing" ) ) {
      // TO DO:  THROW EXCEPTION!!!
    }

    // Get the preprocessing parameter.
    std::string preproc = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(params),
                                                              "Preprocessing" );

    // Create the preprocessor.
    Teuchos::RefCountPtr<Preprocessor< Epetra_MultiVector > > RBPreprocessor;

    // Basic preprocessor does nothing
    if ( preproc == "none" )
      RBPreprocessor = Teuchos::rcp( new NoPreprocessor<Epetra_MultiVector>() );

    // Modified snapshot preprocessor
    if ( preproc == "ModifiedSS" )
      RBPreprocessor = Teuchos::rcp( new MSPreprocessor() );
    //
    // Return the preprocessor created
    //
    return RBPreprocessor;

  }

} // end of RBGen namespace

