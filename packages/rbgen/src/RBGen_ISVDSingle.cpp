#include "RBGen_ISVDSingle.h"
#include "Teuchos_ScalarTraits.hpp"

namespace RBGen {

  ISVDSingle::ISVDSingle() : IncSVDPOD() {
    maxNumPasses_ = 1;
  }
                                          
  int ISVDSingle::makePass() {
    // ISVDSingle only makes a single pass
    TEST_FOR_EXCEPTION(maxNumPasses_ != 1,std::logic_error,
        "RBGen::ISVDSingle::makePass(): Max Num Passes should be 1, but is not.");
    // did we already make our one pass?
    if (curNumPasses_ > 0) return -1;
    const int numCols = A_->NumVectors();
    while (numProc_ < numCols) {
      // determine lup
      int lup;
      if (curRank_ == 0) {
        // first step
        lup = startRank_;
      }
      else {
        // this value minimizes overall complexity, assuming fixed rank
        lup = (int)(curRank_ / Teuchos::ScalarTraits<double>::squareroot(2.0));
      }
      // now cap lup via lmin,lmax,maxBasisSize
      // want lup >= lmin
      // need lup <= numCols - numProc
      //      lup <= lmax
      //      lup <= maxBasisSize - curRank
      lup = (lup < lmin_ ? lmin_ : lup);
      lup = (lup > numCols - numProc_ ? numCols - numProc_ : lup);
      lup = (lup > lmax_ ? lmax_ : lup);
      lup = (lup > maxBasisSize_ - curRank_ ? maxBasisSize_ - curRank_ : lup);

      // get view of new vectors
      Teuchos::RefCountPtr<const Epetra_MultiVector> Aplus;
      Teuchos::RefCountPtr<Epetra_MultiVector> Unew;
      Aplus = Teuchos::rcp( new Epetra_MultiVector(::View,*A_,numProc_,lup));
      Unew = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,curRank_,lup));
      // put them in U
      *Unew = *Aplus;
      Unew = Teuchos::null;
      Aplus = Teuchos::null;

      // perform the incremental step
      incStep(lup);
  
      // increment the column pointer
      numProc_ += lup;
    }
    curNumPasses_++;
    return 0;
  }
    
} // end of RBGen namespace


