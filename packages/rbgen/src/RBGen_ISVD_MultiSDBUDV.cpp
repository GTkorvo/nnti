#include "RBGen_ISVD_MultiSDBUDV.h"

namespace RBGen {

  ISVD_MultiSDBUDV::ISVD_MultiSDBUDV() : IncSVDPOD(), ISVDUDV(), ISVDMultiSDB() {}

  void ISVD_MultiSDBUDV::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDMultiSDB::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
