#ifndef RBGEN_ISVD_MULTISDBUDV_H
#define RBGEN_ISVD_MULTISDBUDV_H

#include "RBGen_ISVDMultiSDB.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  class ISVD_MultiSDBUDV : public virtual ISVDUDV, public virtual ISVDMultiSDB {
    public:
      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_MultiSDBUDV();

      //! Destructor.
      virtual ~ISVD_MultiSDBUDV() {};
      //@}

      //! @name Set Methods
      //@{

      //! Initialize the method with the given parameter list and snapshot set.
      void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
          const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
          const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio = Teuchos::null );

      //@}
  };

} // end of RBGen namespace

#endif // RBGEN_ISVD_MULTISDBUDV_H
