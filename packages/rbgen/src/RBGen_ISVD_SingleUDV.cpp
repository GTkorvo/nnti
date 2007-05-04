#include "RBGen_ISVD_SingleUDV.h"

namespace RBGen {

  ISVD_SingleUDV::ISVD_SingleUDV() : ISVDUDV(), ISVDSingle() {}

  void ISVD_SingleUDV::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDSingle::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
