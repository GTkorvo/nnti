#include "RBGen_ISVDMultiCD.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_LAPACK.h"

namespace RBGen {

  ISVDMultiCD::ISVDMultiCD() : IncSVDPOD() {}

  int ISVDMultiCD::makePass() {
    bool firstPass = (curRank_ == 0);
    if (firstPass) {
      if (curNumPasses_ + 1 > maxNumPasses_) return -1;
    }
    else {
      if (curNumPasses_ + 2 > maxNumPasses_) return -1;
    }
    const int numCols = A_->NumVectors();
    numProc_ = 0;

    // compute W = I - Z T Z^T from current V_
    if (!firstPass) {
      Epetra_LAPACK lapack;
      // copy V_ into workZ_
      Teuchos::RefCountPtr<Epetra_MultiVector> lclHV 
        = Teuchos::rcp( new Epetra_MultiVector(::View,*workZ_,0,curRank_) );
      Teuchos::RefCountPtr<Epetra_MultiVector> lclV 
        = Teuchos::rcp( new Epetra_MultiVector(::View,*V_,0,curRank_) );
      *lclHV = *lclV;
      lclV = Teuchos::null;
      // compute the Householder QR factorization of the current right basis
      // Vhat = W*R
      double *VA; 
      int VLDA, info;
      int lwork = curRank_;
      std::vector<double> tau(curRank_), work(lwork);
      info = lclHV->ExtractView(&VA,&VLDA);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector.");
      lapack.GEQRF(numProc_,curRank_,VA,VLDA,&tau[0],&work[0],lwork,&info);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling GEQRF on current right basis while constructing next pass coefficients.");
      // compute the block representation
      // W = I - Z T Z^T
      lapack.LARFT('F','C',numProc_,curRank_,VA,VLDA,&tau[0],workT_->A(),workT_->LDA());
      // LARFT left upper tri block of Z unchanged
      // note: it should currently contain R factor of V_, which is very close to
      //   diag(\pm 1, ..., \pm 1)
      //
      // we need to set it to:
      //   [1 0 0 ... 0]
      //   [  1 0 ... 0]
      //   [   ....    ]
      //   [          1]
      //
      // see documentation for LARFT
      for (int j=0; j<curRank_; j++) {
        VA[j*VLDA+j] = 1.0;
        for (int i=0; i<j; i++) {
          VA[j*VLDA+i] = 0.0;
        }
      }
      // compute part of A W:  A Z T
      // put this in workAZ_
      // first, A Z (this consumes a pass through A)
      // finish
      curNumPasses_++;
      // second, (A Z) T (in situ, as T is upper triangular)
      // finish

      // set curRank_ = 0
      curRank_ = 0;
    }

    while (numProc_ < numCols) {
      // determine lup
      int lup;
      if (curRank_ == 0) {
        // first step uses startRank_
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
      if (firstPass) {
        // new vectors are just Aplus
        *Unew = *Aplus;
      }
      else {
        // new vectors are Aplus - (A Z T) Z_i^T
        // finish
      }
      Unew = Teuchos::null;
      Aplus = Teuchos::null;

      // perform the incremental step
      incStep(lup);

      // increment the column pointer
      numProc_ += lup;
    }


    // compute W V_ = V_ - Z T Z^T V_
    // finish

    // update pass counter
    curNumPasses_++;

    return 0;
  }

  void ISVDMultiCD::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio
      ) 
  {
    IncSVDPOD::Initialize(params,ss,fileio);
    workAZ_ = Teuchos::rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );
    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    workZ_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );
    workT_  = Teuchos::rcp( new Epetra_SerialDenseMatrix(maxBasisSize_,maxBasisSize_) );
  }

} // end of RBGen namespace
