#include "RBGen_IncSVDPOD.h"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"

namespace RBGen {

  IncSVDPOD::IncSVDPOD() : 
    isInitialized_(false),
    filter_(Teuchos::null),
    maxBasisSize_(0),
    curRank_(0),
    sigma_(0),
    numProc_(0),
    maxNumPasses_(-1),
    curNumPasses_(0),
    tol_(1e-14),
    lmin_(0),
    lmax_(0),
    startRank_(0),
    timerComp_("Total Elapsed Time")
  {}

  Teuchos::RefCountPtr<const Epetra_MultiVector> IncSVDPOD::getBasis() const {
    if (curRank_ == 0 || isInitialized_ == false) {
      return Teuchos::null;
    }
    return Teuchos::rcp( new Epetra_MultiVector(::View,*U_,0,curRank_) );
  }

  Teuchos::RefCountPtr<const Epetra_MultiVector> IncSVDPOD::getRightBasis() const {
    if (curRank_ == 0 || isInitialized_ == false) {
      return Teuchos::null;
    }
    return Teuchos::rcp( new Epetra_MultiVector(::View,*V_,0,curRank_) );
  }

  std::vector<double> IncSVDPOD::getSingularValues() const { 
    vector<double> ret(sigma_.begin(),sigma_.begin()+curRank_);
    return ret;
  }

  void IncSVDPOD::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                              const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
                              const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio ) {

    using Teuchos::rcp;

    // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList rbmethod_params = params->sublist( "Reduced Basis Method" );

    // Get the maximum basis size
    maxBasisSize_ = rbmethod_params.get<int>("Max Basis Size");
    TEST_FOR_EXCEPTION(maxBasisSize_ < 2,invalid_argument,"""Max Basis Size"" must be at least 2.");

    // Get a filter
    filter_ = rbmethod_params.get<Teuchos::RefCountPtr<Filter<double> > >("Filter",Teuchos::null);
    if (filter_ == Teuchos::null) {
      int k = rbmethod_params.get("Rank",(int)5);
      filter_ = rcp( new RangeFilter<double>(LARGEST,k,k) );
    }

    // Get convergence tolerance
    tol_ = rbmethod_params.get<int>("Converence Tolerance",tol_);

    // Get an Anasazi orthomanager
    if (rbmethod_params.isType<
          Teuchos::RefCountPtr< Anasazi::OrthoManager<double,Epetra_MultiVector> > 
        >("Ortho Manager")
       ) 
    {
      ortho_ = rbmethod_params.get< 
                Teuchos::RefCountPtr<Anasazi::OrthoManager<double,Epetra_MultiVector> >
               >("Ortho Manager");
      TEST_FOR_EXCEPTION(ortho_ == Teuchos::null,invalid_argument,"User specified null ortho manager.");
    }
    else {
      string omstr = rbmethod_params.get("Ortho Manager","DGKS");
      if (omstr == "SVQB") {
        ortho_ = rcp( new Anasazi::SVQBOrthoManager<double,Epetra_MultiVector,Epetra_Operator>() );
      }
      else { // if omstr == "DGKS"
        ortho_ = rcp( new Anasazi::BasicOrthoManager<double,Epetra_MultiVector,Epetra_Operator>() );
      }
    }

    // Lmin,Lmax,Kstart
    lmin_ = rbmethod_params.get("Min Update Size",1);
    TEST_FOR_EXCEPTION(lmin_ < 1 || lmin_ >= maxBasisSize_,invalid_argument,
                       "Method requires 1 <= min update size < max basis size.");
    lmax_ = rbmethod_params.get("Max Update Size",maxBasisSize_);
    TEST_FOR_EXCEPTION(lmin_ > lmax_,invalid_argument,"Max update size must be >= min update size.");

    startRank_ = rbmethod_params.get("Start Rank",lmin_);
    TEST_FOR_EXCEPTION(startRank_ < 1 || startRank_ > maxBasisSize_,invalid_argument,
                       "Starting rank must be in [1,maxBasisSize_)");
    // MaxNumPasses
    maxNumPasses_ = rbmethod_params.get("Maximum Number Passes",maxNumPasses_);
    TEST_FOR_EXCEPTION(maxNumPasses_ != -1 && maxNumPasses_ <= 0, invalid_argument,
                       "Maximum number passes must be -1 or > 0.");

    // Save the pointer to the snapshot matrix
    TEST_FOR_EXCEPTION(ss == Teuchos::null,invalid_argument,"Input snapshot matrix cannot be null.");
    A_ = ss;

    // Allocate space for the factorization
    sigma_.reserve( maxBasisSize_ );
    U_ = Teuchos::null;
    V_ = Teuchos::null;
    U_ = rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );
    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    V_ = rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );
    B_ = rcp( new Epetra_SerialDenseMatrix(maxBasisSize_,maxBasisSize_) );

    // clear counters
    numProc_ = 0;
    curNumPasses_ = 0;

    // we are now initialized, albeit with null rank
    isInitialized_ = true;
  }

  void IncSVDPOD::Reset( const Teuchos::RefCountPtr<Epetra_MultiVector>& new_ss ) {
    // Reset the pointer for the snapshot matrix
    // Note: We will not assume that it is non-null; user could be resetting our
    // pointer in order to delete the original snapshot set
    A_ = new_ss;
  }

  void IncSVDPOD::computeBasis() {

    //
    // perform enough incremental updates to consume the entire snapshot set
    Teuchos::TimeMonitor lcltimer(timerComp_);

    // check that we have a valid snapshot set: user may have cleared 
    // it using Reset()
    TEST_FOR_EXCEPTION(A_ == Teuchos::null,logic_error,
                       "computeBasis() requires non-null snapshot set.");

    // reset state
    curRank_ = 0;    
    numProc_ = 0;
    curNumPasses_ = 0;

    while (makePass() == 0) {
      // i'm making a pass.
      // FINISH: check convergence
    }
  }


  /*
  // old stuff from computeBasis... will go into makePass for ISVDSingle
    while (numProc_ < numCols) {
      // determine lup
      int lup;
      if (curRank_ == 0) {
        // first step
        lup = startRank_;
      }
      else {
        // this value minimizes overall complexity, assuming fixed rank
        lup = (int)(curRank_ / SCT::squareroot(2.0));
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

      // get view up new vectors
      Teuchos::RefCountPtr<const Epetra_MultiVector> Aplus = 
        Teuchos::rcp( new Epetra_MultiVector(::View,*A_,numProc_,lup));

      // perform the incremental step
      incStep(Aplus);
  
      // increment the column pointer
      numProc_ += lup;
    }
    */

  
  void IncSVDPOD::updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss ) {
    // perform enough incremental updates to consume the new snapshots

    /*
    typedef Teuchos::ScalarTraits<double> SCT;

    Teuchos::TimeMonitor lcltimer(timerComp_);

    TEST_FOR_EXCEPTION(update_ss == Teuchos::null,invalid_argument,
                       "updateBasis() requires non-null snapshot set.");

    const int numCols = update_ss->NumVectors();

    int i = 0;
    while (i < numCols) {
      // determine lup
      // this value minimizes overall complexity, assuming fixed rank
      int lup = (int)(curRank_ / SCT::squareroot(2.0));
      // now cap lup via lmin,lmax,maxBasisSize
      // want lup >= lmin
      // need lup <= numCols - i
      //      lup <= lmax
      //      lup <= maxBasisSize - curRank
      lup = (lup > numCols - i ? numCols - i : lup);
      lup = (lup < lmin_ ? lmin_ : lup);
      lup = (lup > lmax_ ? lmax_ : lup);
      lup = (lup > maxBasisSize_ - curRank_ ? maxBasisSize_ - curRank_ : lup);

      // get view up new vectors
      Teuchos::RefCountPtr<const Epetra_MultiVector> Aplus = 
        Teuchos::rcp( new Epetra_MultiVector(::View,*update_ss,i,lup));

      // perform the incremental step
      incStep(Aplus);
  
      // increment the column pointer
      i += lup;
      numProc_ += lup;
    }
    */

  }

  void IncSVDPOD::incStep(int lup) {

    /*
    typedef Teuchos::SerialDenseMatrix<int,double> TSDM;
    int newRank;

    //
    // store new vectors in U
    Teuchos::RefCountPtr< Epetra_MultiVector > U2; 
    Teuchos::RefCountPtr< const Epetra_MultiVector > curU;
    if (curRank_ > 0) {
      curU = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,0       ,curRank_) );
    }
    U2   = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,curRank_,lup     ) );
    *U2 = *Aplus;

    //
    // build R and perform gram-schmidt expansion
    //
    //      k l
    // R = [S C] k
    //     [  B] l
    //
    Epetra_SerialDenseMatrix R(curRank_+lup,curRank_+lup);
    for (int i=0; i<curRank_; ++i) {
      R(i,i) = sigma_[i];
    }
    // get pointer for C,B inside of R, as Teuchos::SerialDenseMatrix objects
    Teuchos::RefCountPtr<TSDM> C, B;
    if (curRank_ > 0) {
      C = Teuchos::rcp( new TSDM(Teuchos::View, &R(0,curRank_), R.LDA(), curRank_, lup) );
    }
    B = Teuchos::rcp( new TSDM(Teuchos::View, &R(curRank_,curRank_), R.LDA(), lup,      lup) );
    // perform Grams-Schmidt expansion
    if (curRank_ > 0) {
      newRank = ortho_->projectAndNormalize(*U2,Teuchos::tuple(C),B,Teuchos::tuple(curU));
    }
    else {
      newRank = ortho_->normalize(*U2,B);
    }
    TEST_FOR_EXCEPTION(newRank != lup,logic_error,
                       "IncSVDPOD::incStep(): Couldn't recover full rank basis.");
    C = Teuchos::null;
    B = Teuchos::null;
    U2 = Teuchos::null;
    curU = Teuchos::null;

    //
    // compute SVD of R
    vector<double> newsig(R.M());
    Epetra_LAPACK lapack;
    vector<double> work(5*R.M());
    int lwork = work.size();
    int info;
    lapack.GESVD('O','N',R.M(),R.N(),R.A(),R.LDA(),
                 &newsig[0],(double*)NULL,1,(double*)NULL,1,
                 &(work[0]),&lwork,&info);
    TEST_FOR_EXCEPTION(info != 0,logic_error,"IncSVDPOD::incStep(): Error calling DGESVD.");

    //
    // filter desired sigmas
    vector<int> f = filter_->filter(newsig);
    newRank = f.size();
    sigma_.resize(newRank);
    Epetra_SerialDenseMatrix RU1(curRank_+lup,newRank);
    for (int i=0; i<newRank; i++) {
      sigma_[i] = newsig[f[i]];
      copy(R[f[i]],R[f[i]]+R.M(),RU1[i]);
    }
    // put RU1 into an Epetra MultiVector
    Epetra_LocalMap LocalMap(RU1.M(), 0, U_->Map().Comm());
    Epetra_MultiVector RU1_Pvec(::Copy, LocalMap, RU1.A(), RU1.LDA(), RU1.N());

    //
    // update bases
    Teuchos::RefCountPtr<Epetra_MultiVector> newwU, fullU, newU;
    fullU = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,0,curRank_+lup) );
    newwU = Teuchos::rcp( new Epetra_MultiVector(::View,*workU_,0,newRank) );
    // multiply by RU1
    info = newwU->Multiply('N','N',1.0,*fullU,RU1_Pvec,0.0);
    TEST_FOR_EXCEPTION(info != 0,logic_error,"IncSVDPOD::incStep(): Error calling EMV::Multiply.");
    fullU = Teuchos::null;
    newU = Teuchos::rcp( new Epetra_MultiVector(::View,*U_,0,newRank) );
    *newU = *newwU;
    newU = Teuchos::null;
    newwU = Teuchos::null;
    curRank_ = newRank;
    */
  }

  
} // end of RBGen namespace
