#include "RBGen_EpetraMVMethodFactory.h"
#include "RBGen_LapackPOD.h"

#if HAVE_RBGEN_ANASAZI
#include "RBGen_AnasaziPOD.h"
#endif

namespace RBGen {
  
  Teuchos::RefCountPtr<Method<Epetra_MultiVector> > EpetraMVMethodFactory::create( const Teuchos::ParameterList& params )
  {
     // See if the "Reduced Basis Method" sublist exists
    if ( !params.isSublist( "Reduced Basis Method" ) ) {
      //  TO DO:  THROW EXCEPTION!!!!
    }

    // Get the "Reduced Basis Method" sublist.
    const Teuchos::ParameterList& rbmethod_params = params.sublist( "Reduced Basis Method" );

    // Get the file format type
    std::string method = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(rbmethod_params),
                                                             "Method" );

    Teuchos::RefCountPtr< Method<Epetra_MultiVector> > RBMethod;

    // POD computed using exact SVD through LAPACK
    if ( method == "LAPACK POD" ) {
      RBMethod = Teuchos::rcp( new LapackPOD() );
    } else
    // Inexact POD computed using inexact SVD through Anasazi
#if HAVE_RBGEN_ANASAZI
    if ( method == "Anasazi POD" ) {
      RBMethod = Teuchos::rcp( new AnasaziPOD() );
    } else 
#endif
    {
	// TO DO:  Throw exception because method was not valid
    }
    //
    // Return the method created
    //
    return RBMethod;
  }  

} // end of RBGen namespace

