#include "RBGen_EpetraMVMethodFactory.h"
#include "RBGen_LapackPOD.h"
#include "RBGen_AnasaziPOD.h"
#include <string>

namespace RBGen {
  
  Teuchos::RefCountPtr<Method<Epetra_MultiVector> > EpetraMVMethodFactory::create( const Teuchos::ParameterList& params )
  {
    if ( !params.isParameter( "Reduced Basis Method" ) ) {
      //  TO DO:  THROW EXCEPTION!!!!
    }

    // Get the file format type
    std::string method = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(params),
                                                             "Reduced Basis Method" );

    Teuchos::RefCountPtr< Method<Epetra_MultiVector> > RBMethod;

    // POD computed using exact SVD through LAPACK
    if ( method == "LAPACK POD" )
      RBMethod = Teuchos::rcp( new LapackPOD() );
    
    // Inexact POD computed using inexact SVD through Anasazi
    if ( method == "Anasazi POD" )
      RBMethod = Teuchos::rcp( new AnasaziPOD() );
    //
    // Return the method created
    //
    return RBMethod;
  }  

} // end of RBGen namespace

