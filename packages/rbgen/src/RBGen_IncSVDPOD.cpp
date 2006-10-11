#include "RBGen_IncSVDPOD.h"

namespace RBGen {

  IncSVDPOD::IncSVDPOD() : 
    isInitialized_(false),
    filter_(Teuchos::null),
    max_basis_size_(0),
    cur_rank_(0),
    comp_time_(0),
    U_(Teuchos::null),
    V_(Teuchos::null),
    sigma_(0),
    twoSided_(false),
    num_proc_(0) 
  {}
                      
  void IncSVDPOD::computeBasis() {
  }
  
  void IncSVDPOD::updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss ) {
  }

  Teuchos::RefCountPtr<const Epetra_MultiVector> IncSVDPOD::getBasis() const {
    return Teuchos::null;
  }

  Teuchos::RefCountPtr<const Epetra_MultiVector> IncSVDPOD::getRightBasis() const {
    return Teuchos::null;
  }

  void IncSVDPOD::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                   const Teuchos::RefCountPtr< Epetra_MultiVector >& init        ) {
  }

  void IncSVDPOD::Reset( const Teuchos::RefCountPtr<Epetra_MultiVector>& new_ss ) {
  }
  
} // end of RBGen namespace
