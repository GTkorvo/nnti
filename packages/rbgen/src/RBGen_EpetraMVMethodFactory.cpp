#include "RBGen_EpetraMVMethodFactory.h"
#include "RBGen_LapackPOD.h"

#if HAVE_RBGEN_ANASAZI
#include "RBGen_AnasaziPOD.h"
#include "RBGen_ISVD_SingleUDV.h"
#include "RBGen_ISVD_MultiCDUDV.h"
#include "RBGen_ISVD_MultiSDAUDV.h"
#include "RBGen_ISVD_MultiSDBUDV.h"
#include "RBGen_StSVD_RTR.h"
#endif

namespace RBGen {
  
  Teuchos::RCP<Method<Epetra_MultiVector,Epetra_CrsMatrix> > EpetraMVMethodFactory::create( const Teuchos::ParameterList& params )
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

    Teuchos::RCP< Method<Epetra_MultiVector,Epetra_CrsMatrix> > RBMethod;

    // POD computed using exact SVD through LAPACK
    if ( method == "Lapack POD" ) {
      RBMethod = Teuchos::rcp( new LapackPOD() );
    } 
    // Inexact POD computed using inexact SVD through Anasazi
    // IncSVDPOD uses Anasazi utility classes, while AnasaziPOD uses Anasazi for the solution
#if HAVE_RBGEN_ANASAZI
    else if ( method == "IncSVD POD" ) {
      std::string incsvdmethod = rbmethod_params.get<std::string>("IncSVD Method");
      if ( incsvdmethod == "Single/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_SingleUDV() );
      }
      else if ( incsvdmethod == "MultiCD/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_MultiCDUDV() );
      }
      else if ( incsvdmethod == "MultiSDA/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_MultiSDAUDV() );
      }
      else if ( incsvdmethod == "MultiSDB/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_MultiSDBUDV() );
      }
    } 
    else if ( method == "StSVD/RTR") {
      RBMethod = Teuchos::rcp( new StSVDRTR() );
    }
    else if ( method == "Anasazi POD" ) {
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

