#ifndef RBGEN_INCSVD_POD_H
#define RBGEN_INCSVD_POD_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "RBGen_Filter.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziOrthoManager.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"

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

    std::vector<double> getSingularValues() const;

    double getCompTime() const { return timerComp_.totalElapsedTime(); }

    int numSnapshots() const { return numProc_; }

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

    // private member for performing inc steps
    void incStep(const Teuchos::RefCountPtr<const Epetra_MultiVector> &Aplus);

    // Is this object initialized?
    bool isInitialized_;

    // Singular value filter
    Teuchos::RefCountPtr< Filter<double> > filter_;

    // Max allocation size.
    // The maximum rank DURING any step:
    //    lup <= maxBasisSize_ - curRank_
    int maxBasisSize_;

    // Current rank of the factorization
    int curRank_;

    // Pointers to the snapshots and reduced basis.
    Teuchos::RefCountPtr<Epetra_MultiVector> A_, U_, V_, workU_;

    // Vector holding singular values.
    std::vector<double> sigma_;

    // Are we two-sided or left-sided?
    bool twoSided_;

    // Number of snapshots processed thus far
    int numProc_;

    // min,max number of update vectors
    int lmin_;
    int lmax_;
    int startRank_;

    // ortho manager
    Teuchos::RefCountPtr< Anasazi::OrthoManager<double,Epetra_MultiVector> > ortho_;

    // timer
    Teuchos::Time timerComp_;
  };
  
} // end of RBGen namespace

#endif // RBGEN_INCSVD_POD_H
