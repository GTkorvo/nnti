

#ifndef RBGEN_ANASAZI_POD_H
#define RBGEN_ANASAZI_POD_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"

namespace RBGen {
  
  class AnasaziPOD : public virtual Method<Epetra_MultiVector>, public virtual PODMethod<double> {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    AnasaziPOD();

    //! Destructor.
    virtual ~AnasaziPOD() {};
    //@}

    //@{ @name Computation Methods

    void computeBasis();
    
    void updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss ) {};

    //@}

    //@{ @name Get Methods
    
    Teuchos::RefCountPtr<const Epetra_MultiVector> getBasis() const { return basis_; }

    const std::vector<double> getSingularValues() const { return sv_; }

    double getCompTime() const { return comp_time_; }
    //@}
    
    //@{ @name Set Methods
    
    //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                     const Teuchos::RefCountPtr< Epetra_MultiVector >& ss );

    void Reset( const Teuchos::RefCountPtr<Epetra_MultiVector>& new_ss ) { ss_ = new_ss; }

    //@}

    //@{ @name Status Methods

    bool isInitialized() { return isInitialized_; }
 
    //@}

  private:

    // Is this object initialized.
    bool isInitialized_;

    // Size of the basis that this method will compute.
    int basis_size_;

    // Computational time (in seconds, using wall clock).
    double comp_time_;

    // Pointers to the snapshots and reduced basis.
    Teuchos::RefCountPtr<Epetra_MultiVector> ss_, basis_;

    // Vector holding singular values.
    std::vector<double> sv_;
  };
  
} // end of RBGen namespace

#endif // RBGEN_ANASAZI_POD_H
