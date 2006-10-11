#ifndef RBGEN_INCSVD_POD_H
#define RBGEN_INCSVD_POD_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "RBGen_Filter.hpp"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"

namespace RBGen {

  class IncSVDPOD : public virtual Method<Epetra_MultiVector>, public virtual PODMethod<double> {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    IncSVDPOD();

    //! Destructor.
    virtual ~IncSVDPOD() {};
    //@}

    //@{ @name Computation Methods

    void computeBasis();
    
    void updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss );

    //@}

    //@{ @name Get Methods
    
    Teuchos::RefCountPtr<const Epetra_MultiVector> getBasis() const;

    Teuchos::RefCountPtr<const Epetra_MultiVector> getRightBasis() const;

    const std::vector<double> getSingularValues() const { return sigma_; }

    double getCompTime() const { return comp_time_; }

    int numSnapshots() const { return num_proc_; }

    //@}
    
    //@{ @name Set Methods
    
    //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                     const Teuchos::RefCountPtr< Epetra_MultiVector >& init        );

    void Reset( const Teuchos::RefCountPtr<Epetra_MultiVector>& new_ss );

    //@}

    //@{ @name Status Methods

    bool isInitialized() { return isInitialized_; }
 
    //@}

  private:

    // Is this object initialized?
    bool isInitialized_;

    // Singular value filter
    Teuchos::RefCountPtr< Filter<double> > filter_;

    // Max allocation size.
    // The maximum rank DURING any step:
    //    lup <= max_basis_size_ - cur_rank_
    int max_basis_size_;

    // Current rank of the factorization
    int cur_rank_;

    // Computational time (in seconds, using wall clock).
    double comp_time_;

    // Pointers to the snapshots and reduced basis.
    Teuchos::RefCountPtr<Epetra_MultiVector> U_, V_;

    // Vector holding singular values.
    std::vector<double> sigma_;

    // Are we two-sided or left-sided?
    bool twoSided_;

    // Number of snapshots processed thus far
    int num_proc_;

  };
  
} // end of RBGen namespace

#endif // RBGEN_INCSVD_POD_H
