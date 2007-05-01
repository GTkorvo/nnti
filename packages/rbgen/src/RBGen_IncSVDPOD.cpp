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
    isInitialized_ = false;
  }

  void IncSVDPOD::computeBasis() {

    //
    // perform enough incremental updates to consume the entire snapshot set
    Teuchos::TimeMonitor lcltimer(timerComp_);

    // check that we have a valid snapshot set: user may have cleared 
    // it using Reset()
    TEST_FOR_EXCEPTION(A_ == Teuchos::null,logic_error,
                       "computeBasis() requires non-null snapshot set.");

    // check that we are initialized, i.e., data structures match the data set
    TEST_FOR_EXCEPTION(isInitialized_==false,std::logic_error,
        "RBGen::IncSVDPOD::computeBasis(): Solver must be initialized.");

    // reset state
    curRank_ = 0;    
    numProc_ = 0;
    curNumPasses_ = 0;

    while (makePass() == 0) {
      // i'm making a pass.
      // FINISH: add convergence check
    }
  }

  
  void IncSVDPOD::updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss ) {
    // perform enough incremental updates to consume the new snapshots
    TEST_FOR_EXCEPTION(true,std::logic_error,
        "RBGen::IncSVDPOD::updateBasis(): this routine not yet supported.");
  }

  void IncSVDPOD::incStep(int lup) {

    Epetra_LAPACK lapack;

    // perform gram-schmidt expansion
    this->expand(lup);
    const int lwork = 5*curRank_;
    int info;
    Epetra_SerialDenseMatrix Uhat(curRank_,curRank_), Vhat(curRank_,curRank_);
    std::vector<double> Shat(curRank_), work(lwork);

    // compute the SVD of B
    // Note: this destroys B and actually stores Vhat^T (we remedy this below)
    lapack.GESVD('A','A',curRank_,curRank_,B_->A(),B_->LDA(),&Shat[0],Uhat.A(),Uhat.LDA(),Vhat.A(),Vhat.LDA(),&work[0],&lwork,&info);
    TEST_FOR_EXCEPTION(info!=0,std::logic_error,"RBGen::IncSVDPOD::incStep(): GESVD return info != 0");

    // use filter to determine new rank
    std::vector<int> keepind = filter_->filter(Shat);
    std::vector<int> truncind(curRank_-keepind.size());
    {
      std::vector<int> allind(curRank_);
      for (int i=0; i<curRank_; i++) {
        allind[i] = i;
      }
      set_difference(allind.begin(),allind.end(),keepind.begin(),keepind.end(),truncind.begin());
      
      // Vhat actually contains Vhat^T; correct this here
      Epetra_SerialDenseMatrix Ucopy(Uhat), Vcopy(curRank_,curRank_); 
      std::vector<double> Scopy(Shat);
      for (int j=0; j<curRank_; j++) {
        for (int i=0; i<curRank_; i++) {
          Vcopy(i,j) = Vhat(j,i);
        }
      }
      // put the desired sigmas at the front of Uhat, Vhat
      for (unsigned int j=0; j<keepind.size(); j++) {
        std::copy(&Ucopy(0,keepind[j]),&Ucopy(curRank_,keepind[j]),&Uhat(0,j));
        std::copy(&Vcopy(0,keepind[j]),&Vcopy(curRank_,keepind[j]),&Vhat(0,j));
        Shat[j] = Scopy[keepind[j]];
      }
      for (unsigned int j=0; j<truncind.size(); j++) {
        std::copy(&Ucopy(0,truncind[j]),&Ucopy(curRank_,truncind[j]),&Uhat(0,keepind.size()+j));
        std::copy(&Vcopy(0,truncind[j]),&Vcopy(curRank_,truncind[j]),&Vhat(0,keepind.size()+j));
        Shat[keepind.size()+j] = Scopy[truncind[j]];
      }
    }

    // shrink back down again
    // shrink() will update bases U_ and V_, as well as singular values sigma_ and curReank_
    this->shrink(truncind.size(),Shat,Uhat,Vhat);
  }
    
} // end of RBGen namespace


