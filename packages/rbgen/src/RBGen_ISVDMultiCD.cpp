#include "RBGen_ISVDMultiCD.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_LAPACK.h"
#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"

namespace RBGen {

  ISVDMultiCD::ISVDMultiCD() {}

  int ISVDMultiCD::makePass() {
    Epetra_LAPACK lapack;
    Epetra_BLAS   blas;

    bool firstPass = (curRank_ == 0);
    // if maxNumPasses == -1, we get infinite passes
    if (maxNumPasses_ != -1) {
      if (firstPass) {
        // we need only one: passing through A
        if (curNumPasses_ + 1 > maxNumPasses_) return -1;
      }
      else {
        // we need two: one for A V T and one for passing through A - (A V T) V^T
        if (curNumPasses_ + 2 > maxNumPasses_) return -1;
      }
    }
    const int numCols = A_->NumVectors();
    TEST_FOR_EXCEPTION( !firstPass && (numProc_ != numCols), std::logic_error,
        "RBGen::ISVDMultiCD::makePass(): after first pass, numProc should be numCols");

    // compute W = I - Z T Z^T from current V_
    Teuchos::RefCountPtr<Epetra_MultiVector> lclAZT, lclZ;
    double *Z_A, *AZT_A;
    int Z_LDA, AZT_LDA;
    int oldRank = 0;
    double Rerr;
    if (!firstPass) {
      // copy V_ into workZ_
      lclAZT = Teuchos::rcp( new Epetra_MultiVector(::View,*workAZT_,0,curRank_) );
      lclZ   = Teuchos::rcp( new Epetra_MultiVector(::View,*workZ_,0,curRank_) );
      Teuchos::RefCountPtr<Epetra_MultiVector> lclV;
      lclV = Teuchos::rcp( new Epetra_MultiVector(::View,*V_,0,curRank_) );
      *lclZ = *lclV;
      //DBG cout << "V (right hand basis)" << endl;
      //DBG lclV->Print(cout);
      lclV = Teuchos::null;
      // compute the Householder QR factorization of the current right basis
      // Vhat = W*R
      int info, lwork = curRank_;
      std::vector<double> tau(curRank_), work(lwork);
      info = lclZ->ExtractView(&Z_A,&Z_LDA);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector Z.");
      lapack.GEQRF(numCols,curRank_,Z_A,Z_LDA,&tau[0],&work[0],lwork,&info);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling GEQRF on current right basis while constructing next pass coefficients.");
      if (debug_) {
        // we just took the QR factorization of a set of orthonormal vectors
        // they should have an R factor which is diagonal, with unit elements (\pm 1)
        // check it
        Rerr = 0;
        for (int j=0; j<curRank_; j++) {
          for (int i=0; i<j; i++) {
            Rerr += abs(Z_A[j*Z_LDA+i]);
          }
          Rerr += abs(abs(Z_A[j*Z_LDA+j]) - 1.0);
        }
        //DBG cout << "Compressed QR of V (right hand basis)" << endl;
        //DBG lclZ->Print(cout);
      }
      // compute the block representation
      // W = I - Z T Z^T
      lapack.LARFT('F','C',numCols,curRank_,Z_A,Z_LDA,&tau[0],workT_->A(),workT_->LDA());
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
        Z_A[j*Z_LDA+j] = 1.0;
        for (int i=0; i<j; i++) {
          Z_A[j*Z_LDA+i] = 0.0;
        }
      }
      // compute part of A W:  A Z T
      // put this in workAZT_
      // first, A Z (this consumes a pass through A)
      info = lclAZT->Multiply('N','N',1.0,*A_,*lclZ,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for A*Z");
      curNumPasses_++;
      // second, (A Z) T (in situ, as T is upper triangular)
      info = lclAZT->ExtractView(&AZT_A,&AZT_LDA);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector AZ.");
      blas.TRMM('R','U','N','N',numCols,curRank_,1.0,workT_->A(),workT_->LDA(),AZT_A,AZT_LDA);
      // set curRank_ = 0
      oldRank  = curRank_;
      curRank_ = 0;
    }

    numProc_ = 0;
    while (numProc_ < numCols) {
      //
      // determine lup
      //
      // want lup >= lmin
      //      lup <= lmax
      // need lup <= numCols - numProc
      //      lup <= maxBasisSize - curRank
      //
      int lup;
      if (curRank_ == 0) {
        // first step uses startRank_
        // this is not affected by lmin,lmax
        lup = startRank_;
      }
      else {
        // this value minimizes overall complexity, assuming fixed rank
        lup = (int)(curRank_ / Teuchos::ScalarTraits<double>::squareroot(2.0));
        // contrain to [lmin,lmax]
        lup = (lup < lmin_ ? lmin_ : lup);
        lup = (lup > lmax_ ? lmax_ : lup);
      }
      //
      // now cap lup via maxBasisSize and the available data
      // these caps apply to all lup, as a result of memory and data constraints
      //
      // available data
      lup = (lup > numCols - numProc_ ? numCols - numProc_ : lup);
      // available memory
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
        // specifically, Aplus - (A Z T) Z(numProc:numProc+lup-1,1:oldRank)^T
        Epetra_LocalMap lclmap(lup,0,A_->Comm());
        Epetra_MultiVector Zi(::View,lclmap,&Z_A[numProc_],Z_LDA,oldRank);
        *Unew = *Aplus;
        int info = Unew->Multiply('N','T',-1.0,*lclAZT,Zi,1.0);
        TEST_FOR_EXCEPTION(info != 0,std::logic_error,
            "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for A*Wi");
      }
      Unew = Teuchos::null;
      Aplus = Teuchos::null;

      // perform the incremental step
      incStep(lup);

      // increment the column pointer
      numProc_ += lup;
    }

    // compute W V = V - Z T Z^T V
    // Z^T V is oldRank x curRank
    // T Z^T V is oldRank x curRank
    // we need T Z^T V in a local Epetra_MultiVector
    if (!firstPass) {
      Teuchos::RefCountPtr<Epetra_MultiVector> lclV;
      double *TZTV_A;
      int TZTV_LDA;
      int info;
      Epetra_LocalMap lclmap(oldRank,0,A_->Comm());
      // get pointer to current V
      lclV = Teuchos::rcp( new Epetra_MultiVector(::View,*V_,0,curRank_) );
      // create space for T Z^T V
      Epetra_MultiVector TZTV(lclmap,curRank_,false);
      // multiply Z^T V
      info = TZTV.Multiply('T','N',1.0,*lclZ,*lclV,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for Z^T V.");
      // get pointer to data in Z^T V
      info = TZTV.ExtractView(&TZTV_A,&TZTV_LDA);
      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector TZTV.");
      // multiply T (Z^T V)
      blas.TRMM('L','U','N','N',oldRank,curRank_,1.0,workT_->A(),workT_->LDA(),TZTV_A,TZTV_LDA);
      // multiply V - Z (T Z^T V)
      info = lclV->Multiply('N','N',-1.0,*lclZ,TZTV,1.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for W V.");
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
    if (comm->MyPID() == 0) {
      cout 
        << "------------- ISVDMultiCD::makePass() -----------" << endl
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
        if (!firstPass) {
          cout << "|DBG      R-I norm: " << Rerr << endl;
        }
      }
    }

    return 0;
  }

  void ISVDMultiCD::Initialize( 
      const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
      const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
      const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio
      ) 
  {
    workAZT_ = Teuchos::rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );

    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    workZ_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );

    workT_  = Teuchos::rcp( new Epetra_SerialDenseMatrix(maxBasisSize_,maxBasisSize_) );
  }

} // end of RBGen namespace
