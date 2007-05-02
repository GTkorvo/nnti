#include "RBGen_ISVD_MultiCDUDV.h"

namespace RBGen {

  ISVD_MultiCDUDV::ISVD_MultiCDUDV() : ISVDUDV(), ISVDMultiCD() {}

  void ISVD_MultiCDUDV::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDMultiCD::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
