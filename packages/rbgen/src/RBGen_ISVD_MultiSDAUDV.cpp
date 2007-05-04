#include "RBGen_ISVD_MultiSDAUDV.h"

namespace RBGen {

  ISVD_MultiSDAUDV::ISVD_MultiSDAUDV() : IncSVDPOD(), ISVDUDV(), ISVDMultiSDA() {}

  void ISVD_MultiSDAUDV::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDMultiSDA::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
