#ifndef RBGEN_STSVD_RTR_H
#define RBGEN_STSVD_RTR_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"


//
// This class encapsulates a dominant SVD computation via manifold
// optimization.
// 
// We are given a snapshot matrix A of dimension m by n. We wish to find the rank-k 
// dominant SVD.
//
// The manifold is M = St(k,m) x St(k,n), where St(k,p) is the manifold of rank-k orthonormal
// bases of R^p (the compact Stiefel manifold).
//
// The objective function is
// f : M --> R 
//   : (U,V) --> trace(U^T A V N), 
// where N is a diagonal matrix with distinct, ascending (or descending) strictly positive elements.
// For example, N = diag(5,4,3,2,1)
// N serves to order the singular vectors in U and V
//
// This solver applies the Riemannian Trust-Region method (Absil, Baker and Gallivan) to maximize f 
// over M.
//

namespace RBGen {

  //! Class for producing a basis using the Incremental SVD
  class StSVDRTR : public virtual Method<Epetra_MultiVector,Epetra_CrsMatrix>, public virtual PODMethod<double> {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    StSVDRTR();

    //! Destructor.
    virtual ~StSVDRTR() {};
    //@}

    //! @name Computation Methods
    //@{

    //! Computes bases for the left and (optionally) right singular subspaces, along with singular vaues.
    void computeBasis();

    //! Update the current basis using a new set of snapshots.
    void updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss );

    //@}

    //! @name Get Methods
    //@{

    //! Return a basis for the left singular subspace.
    Teuchos::RefCountPtr<const Epetra_MultiVector> getBasis() const;

    //! Return a basis for the right singular subspace.
    Teuchos::RefCountPtr<const Epetra_MultiVector> getRightBasis() const;

    //! Return the singular values.
    std::vector<double> getSingularValues() const;

    //! Return the cummulative wall-clock time.
    double getCompTime() const { return timerComp_.totalElapsedTime(); }

    //! Return the scaled residual norms.
    const std::vector<double> & getResNorms();

    //@}

    //! @name Set Methods
    //@{

    //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                     const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
                     const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio = Teuchos::null );

    void Reset( const Teuchos::RefCountPtr<Epetra_MultiVector>& new_ss );

    //@}

    //! @name Status Methods
    //@{

    bool isInitialized() { return isInitialized_; }

    //@}

  protected:

    // subproblem return codes
    enum trRetType {
      UNINITIALIZED = 0,
      MAXIMUM_ITERATIONS,
      NEGATIVE_CURVATURE,
      EXCEEDED_TR,
      KAPPA_CONVERGENCE,
      THETA_CONVERGENCE
    };

    // private members 
    // Riemannian metric
    // g(eta,eta)
    double innerProduct( Teuchos::RefCountPtr< Epetra_MultiVector > etaU, 
                         Teuchos::RefCountPtr< Epetra_MultiVector > etaV );
    // g(eta,zeta)
    double innerProduct( Teuchos::RefCountPtr< Epetra_MultiVector > etaU, 
                         Teuchos::RefCountPtr< Epetra_MultiVector > etaV, 
                         Teuchos::RefCountPtr< Epetra_MultiVector > zetaU, 
                         Teuchos::RefCountPtr< Epetra_MultiVector > zetaV );
    // Retraction, in situ
    void retract( Teuchos::RefCountPtr< Epetra_MultiVector > xU, 
                  Teuchos::RefCountPtr< Epetra_MultiVector > xV, 
                  Teuchos::RefCountPtr< Epetra_MultiVector > etaU, 
                  Teuchos::RefCountPtr< Epetra_MultiVector > etaV );
    // Compute gradient
    void grad( Teuchos::RefCountPtr< Epetra_MultiVector > xU, 
               Teuchos::RefCountPtr< Epetra_MultiVector > xV, 
               Teuchos::RefCountPtr< Epetra_MultiVector > gradU, 
               Teuchos::RefCountPtr< Epetra_MultiVector > gradV );
    // Apply Hessian
    void Hess( Teuchos::RefCountPtr< Epetra_MultiVector > xU, 
               Teuchos::RefCountPtr< Epetra_MultiVector > xV, 
               Teuchos::RefCountPtr< Epetra_MultiVector > etaU, 
               Teuchos::RefCountPtr< Epetra_MultiVector > etaV );
    // Project back to tangent plane
    void Proj( Teuchos::RefCountPtr< Epetra_MultiVector > xU, 
               Teuchos::RefCountPtr< Epetra_MultiVector > xV, 
               Teuchos::RefCountPtr< Epetra_MultiVector > etaU, 
               Teuchos::RefCountPtr< Epetra_MultiVector > etaV );
    // solve the trust-region subproblem
    void solveTRSubproblem();


    // Is this object initialized?
    bool isInitialized_;

    // Vector holding N
    std::vector<double> N_;

    // Pointer to the snapshots
    Teuchos::RefCountPtr<Epetra_MultiVector> A_;
    // state multivecs
    Teuchos::RefCountPtr<Epetra_MultiVector> U_, V_,           // current iterate
                                             RU_, RV_,         // residual, gradient and residual of model minimization
                                             etaU_, etaV_,     // subproblem solution
                                             deltaU_, deltaV_; // search directions

    // Vector holding singular values.
    std::vector<double> sigma_;

    // Convergence tolerance
    double tol_;

    // ortho manager: need this for qf()
    Teuchos::RefCountPtr< Anasazi::OrthoManager<double,Epetra_MultiVector> > ortho_;

    // cummulative timer
    Teuchos::Time timerComp_;

    // debug flag
    bool debug_;

    // verb level
    int verbLevel_;

    // residual norms
    std::vector<double> resNorms_;

    // trust-region state
    // initial and current trust-region radius
    double Delta0_, Delta_, Delta_bar_;
    // acceptance parameter
    double rhoPrime_;
    // norm of the initial gradient
    double normGrad0_;
    // number of outer iterations
    int iter_;
    // convergence parameters
    double conv_kappa_, conv_theta_;
    // most recent rho
    double rho_;
    // current objective function
    double fx_;
    // last inner stop
    int innerStop_;
    // dimensions of problem
    int m_, n_;
  };

} // end of RBGen namespace

#endif // RBGEN_STSVD_RTR_H
