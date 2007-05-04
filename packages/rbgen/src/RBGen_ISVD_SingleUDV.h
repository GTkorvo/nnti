#ifndef RBGEN_ISVD_SINGLEUDV_H
#define RBGEN_ISVD_SINGLEUDV_H

#include "RBGen_ISVDSingle.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  class ISVD_SingleUDV : public virtual ISVDUDV, public virtual ISVDSingle {
    public:

      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_SingleUDV();

      //! Destructor.
      virtual ~ISVD_SingleUDV() {};
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

#endif // RBGEN_ISVD_SINGLEUDV_H
