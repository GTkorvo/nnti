#include "RBGen_ISVDSingle.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_Comm.h"

namespace RBGen {

  ISVDSingle::ISVDSingle() {}
                                          
  void ISVDSingle::makePass() {
    // ISVDSingle only makes a single pass
    TEST_FOR_EXCEPTION(maxNumPasses_ != 1,std::logic_error,
        "RBGen::ISVDSingle::makePass(): Max Num Passes should be 1, but is not.");
    // did we already make our one pass? we can't afford to run this again
    if (curNumPasses_ > 0) return;
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

    //
    // compute the new residuals
    // we know that A V = U S
    // if, in addition, A^T U = V S, then have singular subspaces
    // check residuals A^T U - V S, scaling the i-th column by sigma[i]
    //
    {
      Epetra_LocalMap lclmap(A_->NumVectors(),0,A_->Comm());
      Epetra_MultiVector ATU(lclmap,maxBasisSize_,false);

      // we know that A V = U S
      // if, in addition, A^T U = V S, then have singular subspaces
      // check residuals A^T U - V S, scaling the i-th column by sigma[i]
      Epetra_MultiVector ATUlcl(::View,ATU,0,curRank_);
      Epetra_MultiVector Ulcl(::View,*U_,0,curRank_);
      Epetra_MultiVector Vlcl(::View,*V_,0,curRank_);
      // compute A^T U
      int info = ATUlcl.Multiply('T','N',1.0,*A_,Ulcl,0.0);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply for A^T U.");
      Epetra_LocalMap rankmap(curRank_,0,A_->Comm());
      Epetra_MultiVector S(rankmap,curRank_,true);
      for (int i=0; i<curRank_; i++) {
        S[i][i] = sigma_[i];
      }
      // subtract V S from A^T U
      info = ATUlcl.Multiply('N','N',-1.0,Vlcl,S,1.0);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::computeBasis(): Error calling Epetra_MultiVector::Multiply for V S.");
      resNorms_.resize(curRank_);
      ATUlcl.Norm2(&resNorms_[0]);
      // scale by sigmas
      for (int i=0; i<curRank_; i++) {
        if (sigma_[i] != 0.0) {
          resNorms_[i] /= sigma_[i];
        }
      }

    }


    // debugging checks
    std::vector<double> errnorms(curRank_);
    if (debug_) {
      int info;
      // Check that A V = U Sigma
      // get pointers to current U and V, create workspace for A V - U Sigma
      Epetra_MultiVector work(U_->Map(),curRank_,false), 
                         curU(::View,*U_,0,curRank_),
                         curV(::View,*V_,0,curRank_);
      // create local MV for sigmas
      Epetra_LocalMap lclmap(curRank_,0,A_->Comm());
      Epetra_MultiVector curS(lclmap,curRank_,true);
      for (int i=0; i<curRank_; i++) {
        curS[i][i] = sigma_[i];
      }
      info = work.Multiply('N','N',1.0,curU,curS,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for debugging U S.");
      info = work.Multiply('N','N',-1.0,*A_,curV,1.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for debugging U S - A V.");
      work.Norm2(&errnorms[0]);
      for (int i=0; i<curRank_; i++) {
        if (sigma_[i] != 0.0) {
          errnorms[i] /= sigma_[i];
        }
      }
    }


    // update pass counter
    curNumPasses_++;

    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    if (comm->MyPID() == 0 && verbLevel_ >= 1) {
      cout 
        << "------------- ISVDSingle::makePass() -----------" << endl
        << "| Number of passes: " << curNumPasses_ << endl
        << "|     Current rank: " << curRank_ << endl
        << "|   Current sigmas: " << endl;
      for (int i=0; i<curRank_; i++) {
        cout << "|             " << sigma_[i] << endl;
      }
      if (debug_) {
        cout << "|DBG   US-AV norms: " << endl;
        for (int i=0; i<curRank_; i++) {
          cout << "|DBG          " << errnorms[i] << endl;
        }
      }
    }

    return;
  }
    
  void ISVDSingle::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio
      ) 
  {
    maxNumPasses_ = 1;
  }

} // end of RBGen namespace


