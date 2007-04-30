#include "RBGen_ISVDUDV.h"

namespace RBGen {

  ISVDUDV::ISVDUDV() : IncSVDPOD() {}

  void ISVDUDV::expand(int lup) 
  {
    //
    // build B and perform gram-schmidt expansion
    //
    //      k l
    // B = [S C] k
    //     [  Z] l
    //
    Teuchos::RefCountPtr<const Epetra_MultiVector> U1;
    Teuchos::RefCountPtr<Epetra_MultiVector> U2;
    U2 = Teuchos::rcp( new Epetra_MultiVector(View,*U_,curRank_,lup) );
    int info = (*B_).Scale(0.0);
    TEST_FOR_EXCEPTION(info != 0,std::logic_error,
        "RBGen::ISVDUDV::expand(): error calling Epetra_SerialDenseMatrix::scale()");
    for (int i=0; i<curRank_; ++i) {
      (*B_)(i,i) = sigma_[i];
    }
    // get pointer for C,B inside of B, as Teuchos::SerialDenseMatrix objects
    Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> C, Z;
    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,double> > Cteuchos, Zteuchos;
    if (curRank_ > 0) {
      U1 = Teuchos::rcp( new Epetra_MultiVector(View,*U_,0,curRank_) );
      C = Teuchos::rcp( new Epetra_SerialDenseMatrix(View, &((*B_)(0,curRank_)), B_->LDA(), curRank_, lup) );
      Cteuchos = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,double>(Teuchos::View,C->A(),C->LDA(),curRank_,lup) );
    }
    Z = Teuchos::rcp( new Epetra_SerialDenseMatrix(View, &((*B_)(curRank_,curRank_)), B_->LDA(), lup, lup) );
    Zteuchos = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,double>(Teuchos::View, Z->A(), Z->LDA(), lup, lup ) );
    // perform Grams-Schmidt expansion
    int newRank;
    if (curRank_ > 0) {
      newRank = ortho_->projectAndNormalize(*U2,Teuchos::tuple(Cteuchos),Zteuchos,Teuchos::tuple(U1));
    }
    else {
      newRank = ortho_->normalize(*U2,Zteuchos);
    }
    TEST_FOR_EXCEPTION(newRank != lup,logic_error,
                       "RBGen::ISVDUDV::incStep(): Couldn't recover full rank basis.");
    Cteuchos = Teuchos::null;
    Zteuchos = Teuchos::null;
    C = Teuchos::null;
    Z = Teuchos::null;
    U2 = Teuchos::null;
    U1 = Teuchos::null;
    
    curRank_ += lup;
  }


  void ISVDUDV::shrink(int down, std::vector<double> &S, Epetra_SerialDenseMatrix &U, Epetra_SerialDenseMatrix &V) 
  {
    //
    // put RU1 into an Epetra MultiVector
    Epetra_LocalMap LocalMap(curRank_, 0, U_->Map().Comm());
    Epetra_MultiVector Uh1(::View, LocalMap, U.A(), U.LDA(), curRank_-down);
    Epetra_MultiVector Vh1(::View, LocalMap, V.A(), V.LDA(), curRank_-down);

    //
    // update bases
    Teuchos::RefCountPtr<Epetra_MultiVector> newwU, fullU, newU, newwV, fullV, newV;
    fullU = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,0,curRank_) );
    fullV = Teuchos::rcp( new Epetra_MultiVector(::View,*V_,0,curRank_) );
    newwU = Teuchos::rcp( new Epetra_MultiVector(::View,*workU_,0,curRank_-down) );
    newwV = Teuchos::rcp( new Epetra_MultiVector(::View,*workV_,0,curRank_-down) );
    // multiply by U1
    int info = newwU->Multiply('N','N',1.0,*fullU,Uh1,0.0);
    TEST_FOR_EXCEPTION(info != 0,logic_error,"ISVDUDV::shrink(): Error calling EMV::Multiply(U).");
    fullU = Teuchos::null;
    newU = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,0,curRank_-down) );
    *newU = *newwU;
    newU = Teuchos::null;
    newwU = Teuchos::null;

    // multiply by V1
    info = newwV->Multiply('N','N',1.0,*fullV,Vh1,0.0);
    TEST_FOR_EXCEPTION(info != 0,logic_error,"ISVDUDV::shrink(): Error calling EMV::Multiply(V).");
    fullV = Teuchos::null;
    newV = Teuchos::rcp( new Epetra_MultiVector(::View,*V_,0,curRank_-down) );
    *newV = *newwV;
    newV = Teuchos::null;
    newwV = Teuchos::null;

    // save new singular values
    for (int i=0; i<curRank_-down; i++) {
      sigma_[i] = S[i];
    }

    curRank_ = curRank_-down;
  }


  void ISVDUDV::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                            const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
                            const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio) 
  {
    IncSVDPOD::Initialize(params,ss,fileio);
    workU_ = Teuchos::rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );
    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    workV_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );
  }

} // end of RBGen namespace

