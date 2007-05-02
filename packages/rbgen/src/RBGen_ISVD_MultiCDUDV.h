#ifndef RBGEN_ISVD_MULTICDUDV_H
#define RBGEN_ISVD_MULTICDUDV_H

#include "RBGen_ISVDMultiCD.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  class ISVD_MultiCDUDV : public virtual ISVDUDV, public virtual ISVDMultiCD {
    public:
      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_MultiCDUDV();

      //! Destructor.
      virtual ~ISVD_MultiCDUDV() {};
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

#endif // RBGEN_ISVD_MULTICDUDV_H
