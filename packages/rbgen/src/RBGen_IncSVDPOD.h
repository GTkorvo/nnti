#ifndef RBGEN_INCSVD_POD_H
#define RBGEN_INCSVD_POD_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "RBGen_Filter.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziOrthoManager.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"


//
// computeBasis()
//       |
//   makePass()     ___
//   while () {    /    expand()
//      incStep()  ---- SVD()
//   }             \___ contract()
//
// makePass(), expand() and contract() are pure virtual
//
// expand(),contract() are implemented in a base class that 
// decides the representation: UDV, QRW, QBW
// 
// makePass() is implemented in a base class that decides 
// if a method is Multipass, and if so, in what manner
//
//                                   IncSVDPOD
//                                       |
//      -------------------------------------------------------------------
//      |       |        |           |            |            |          |
//  ISVDUDV  ISVDQRW  ISVDQBW   ISVDMultiCD ISVDMultiSDA ISVDMultiSDB ISVDSingle
//      |       |        |           |            |            |          |
//      ------------------           --------------------------------------
//              \                                       /
//               \                                     /
//                \                                   /
//                 \                                 /
//                  \                               /
//                   \---  Concrete Base Class ----/
//
// Then a concrete base class (one of 3x4==12 varieties) is formed simply through
// inheritence. This is the Template Pattern type of Design Pattern.
//

namespace RBGen {

  //! Class for producing a basis using the Incremental SVD
  class IncSVDPOD : public virtual Method<Epetra_MultiVector,Epetra_CrsMatrix>, public virtual PODMethod<double> {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    IncSVDPOD();

    //! Destructor.
    virtual ~IncSVDPOD() {};
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

    // SerialDenseMatrix holding current core matrix B
    Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> B_;

    // Vector holding singular values.
    std::vector<double> sigma_;

    // Are we two-sided or left-sided?
    bool twoSided_;

    // Number of snapshots processed thus far
    int numProc_;

    // Maximum allowable number of passes through A
    int maxNumPasses_;

    // Current number of passes through A
    int curNumPasses_;

    // Convergence tolerance
    double tol_;

    // min,max number of update vectors
    int lmin_;
    int lmax_;
    int startRank_;

    // ortho manager
    Teuchos::RefCountPtr< Anasazi::OrthoManager<double,Epetra_MultiVector> > ortho_;

    // cummulative timer
    Teuchos::Time timerComp_;
  };

} // end of RBGen namespace

#endif // RBGEN_INCSVD_POD_H
