#include "RBGen_IncSVDPOD.h"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

namespace RBGen {

  StSVDRTR::StSVDRTR() : 
    isInitialized_(false),
    sigma_(0),
    tol_(1e-12),
    timerComp_("Total Elapsed Time"),
    debug_(false),
    verbLevel_(0),
    resNorms_(0)
  {}

  Teuchos::RefCountPtr<const Epetra_MultiVector> StSVDRTR::getBasis() const {
    if (isInitialized_ == false) {
      return Teuchos::null;
    }
    return U_;
  }

  Teuchos::RefCountPtr<const Epetra_MultiVector> StSVDRTR::getRightBasis() const {
    if (isInitialized_ == false) {
      return Teuchos::null;
    }
    return V_;
  }

  std::vector<double> StSVDRTR::getSingularValues() const { 
    return ret;
  }

  void StSVDRTR::Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                             const Teuchos::RefCountPtr< Epetra_MultiVector >& ss,
                             const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio ) {

    using Teuchos::rcp;

    // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList rbmethod_params = params->sublist( "Reduced Basis Method" );

    // Get rank
    int rank_ = rbmethod_params.get("Rank",(int)5);

    // Get convergence tolerance
    tol_ = rbmethod_params.get<double>("Convergence Tolerance",tol_);

    // Get debugging flag
    debug_ = rbmethod_params.get<bool>("IncSVD Debug",debug_);

    // Get verbosity level
    verbLevel_ = rbmethod_params.get<int>("IncSVD Verbosity Level",verbLevel_);

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

    // Save the pointer to the snapshot matrix
    TEST_FOR_EXCEPTION(ss == Teuchos::null,invalid_argument,"Input snapshot matrix cannot be null.");
    A_ = ss;

    // Allocate space for the factorization
    sigma_.resize( rank_ );
    U_ = Teuchos::null;
    V_ = Teuchos::null;
    U_ = rcp( new Epetra_MultiVector(ss->Map(),rank_,false) );
    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    V_ = rcp( new Epetra_MultiVector(lclmap,rank_,false) );
    resNorms_.resize(rank_);

    // we are now initialized, albeit with null rank
    isInitialized_ = true;
  }

  void StSVDRTR::Reset( const Teuchos::RefCountPtr<Epetra_MultiVector>& new_ss ) {
    // Reset the pointer for the snapshot matrix
    // Note: We will not assume that it is non-null; user could be resetting our
    // pointer in order to delete the original snapshot set
    A_ = new_ss;
    isInitialized_ = false;
  }

  void StSVDRTR::computeBasis() {

    //
    // perform enough incremental updates to consume the entire snapshot set
    Teuchos::TimeMonitor lcltimer(timerComp_);

    // check that we have a valid snapshot set: user may have cleared 
    // it using Reset()
    TEST_FOR_EXCEPTION(A_ == Teuchos::null,logic_error,
                       "computeBasis() requires non-null snapshot set.");

    // check that we are initialized, i.e., data structures match the data set
    TEST_FOR_EXCEPTION(isInitialized_==false,std::logic_error,
        "RBGen::StSVDRTR::computeBasis(): Solver must be initialized.");

    //
    // reset state: finish

    //
    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    while (true) {  // not converged: finish

      // compute gradient: finish

      // check convergence: finish

      // minimize model subproblem: finish

      // compute proposed point: finish

      // evaluate rho: finish

      // accept/reject and adjust trust-region radius: finish

    }
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // update the basis with new snapshots: not supported
  void StSVDRTR::updateBasis( const Teuchos::RefCountPtr< Epetra_MultiVector >& update_ss ) {
    // perform enough incremental updates to consume the new snapshots
    TEST_FOR_EXCEPTION(true,std::logic_error,
        "RBGen::StSVDRTR::updateBasis(): this routine not yet supported.");
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return residual norms
  const std::vector<double> & StSVDRTR::getResNorms() {
    return resNorms_;
  }

    
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TR subproblem solver
  void StSVDRTR::solveTRSubproblem() {

    // return one of:
    // MAXIMUM_ITERATIONS
    // NEGATIVE_CURVATURE
    // EXCEEDED_TR
    // KAPPA_CONVERGENCE
    // THETA_CONVERGENCE

    innerStop_ = MAXIMUM_ITERATIONS;

    // finish
    const int d = (n-k)*k + (m-k)*k + (k-1)*k;
    const double D2 = Delta_*Delta_;

    double d_Hd, alpha, beta, z_r, zold_rold;
    double r0_norm;

    // set eta to zero
    MVT::MvInit(*this->eta_);

    //
    // R_ is A X_ - B X_ diag(theta_)
    //
    // r0 = grad f(X) = 2 P_BX A X = 2 P_BX (A X - B X diag(theta_) = 2 proj(R_)
    // We will do this in place.
    //
    {
      Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
      this->orthman_->project(*this->R_,Teuchos::null,
                              Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                              Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->X_),
                              Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->BX_));
      MVT::MvScale(*this->R_,2.0);
    }
    r0_norm = MAT::squareroot( ginner(*this->R_) );

    // z = Prec^-1 r
    // this must be in the tangent space
    if (this->hasPrec_) {
      Teuchos::TimeMonitor lcltimer( *this->timerPrec_ );
      OPT::Apply(*this->Prec_,*this->R_,*this->Z_);
      // the orthogonalization time counts under Ortho and under Preconditioning
      {
        Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
        this->orthman_->project(*this->Z_,Teuchos::null,
                                Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                                Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->X_),
                                Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->BX_));
      }
      z_r = SCT::real( ginner(*this->R_,*this->Z_) );
    }
    else {
      z_r = SCT::real( ginner(*this->R_) );
    }
    // initial value of <delta,delta>_P
    d_Pd = z_r;

    // delta = -z
    MVT::MvAddMv(-ONE,*this->Z_,ZERO,*this->Z_,*this->delta_);

    // the loop
    for (int i=0; i<d; i++) {

      // 
      // print some debugging info
      if (this->om_->isVerbosity(Debug)) {
        this->om_->stream(Debug) 
          << " Debugging checks: RTR inner iteration " << i << endl
          << " >> (eta  ,eta  )_Prec : " << e_Pe << endl
          << " >> (eta  ,delta)_Prec : " << e_Pd << endl
          << " >> (delta,delta)_Prec : " << d_Pd << endl;
      }

      // [Hdelta,Adelta,Bdelta] = Hess*d = 2 Proj(A*d - B*d*x'*A*x)
      {
        Teuchos::TimeMonitor lcltimer( *this->timerAOp_ );
        OPT::Apply(*this->AOp_,*this->delta_,*this->Adelta_);
        this->counterAOp_ += this->blockSize_;
      }
      if (this->hasBOp_) {
        Teuchos::TimeMonitor lcltimer( *this->timerBOp_ );
        OPT::Apply(*this->BOp_,*this->delta_,*this->Bdelta_);
        this->counterBOp_ += this->blockSize_;
      }
      // put 2*A*d - 2*B*d*theta --> Hd
      MVT::MvAddMv(ONE,*this->Bdelta_,ZERO,*this->Bdelta_,*this->Hdelta_);
      MVT::MvScale(*this->Hdelta_,this->theta_);
      MVT::MvAddMv(2.0,*this->Adelta_,-2.0,*this->Hdelta_,*this->Hdelta_);
      // apply projector
      {
        Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
        this->orthman_->project(*this->Hdelta_,Teuchos::null,
                                Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                                Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->X_),
                                Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->BX_));
      }
      d_Hd = SCT::real( ginner(*this->delta_,*this->Hdelta_) );

      // compute update step
      alpha = z_r/d_Hd;

      if (this->om_->isVerbosity(Debug)) {
        this->om_->stream(Debug) 
          << " >> (z,r)  : " << z_r  << endl
                                << " >> (d,Hd) : " << d_Hd << endl
                                                      << " >> alpha  : " << alpha << endl;
      }

      // <neweta,neweta>_P = <eta,eta>_P + 2*alpha*<eta,delta>_P + alpha*alpha*<delta,delta>_P
      e_Pe_new = e_Pe + 2.0*alpha*e_Pd + alpha*alpha*d_Pd;

      // check truncation criteria: negative curvature or exceeded trust-region
      if (d_Hd <= 0 || e_Pe_new >= D2) {
        double tau = (-e_Pd + MAT::squareroot(e_Pd*e_Pd + d_Pd*(D2-e_Pe))) / d_Pd;
        if (this->om_->isVerbosity(Debug)) {
          this->om_->stream(Debug) 
            << " >> tau  : " << tau << endl;
        }
        MVT::MvAddMv(tau,*this->delta_ ,ONE,*this->eta_ ,*this->eta_);
        MVT::MvAddMv(tau,*this->Adelta_,ONE,*this->Aeta_,*this->Aeta_);
        if (this->hasBOp_) {
          MVT::MvAddMv(tau,*this->Bdelta_,ONE,*this->Beta_,*this->Beta_);
        }
        MVT::MvAddMv(tau,*this->Hdelta_,ONE,*this->Heta_,*this->Heta_);
        if (d_Hd <= 0) {
          innerStop_ = NEGATIVE_CURVATURE;
        }
        else {
          innerStop_ = EXCEEDED_TR;
        }
        if (this->om_->isVerbosity(Debug)) {
          e_Pe_new = e_Pe + 2.0*tau*e_Pd + tau*tau*d_Pd;
          this->om_->stream(Debug) 
            << " >> new (eta,eta)_Prec : " << e_Pe_new << endl
            << " >> (eta,eta)          : " << ginner(*this->eta_) << endl;
        }
        break;
      }

      // compute new eta == eta + alpha*delta
      MVT::MvAddMv(alpha,*this->delta_ ,ONE,*this->eta_ ,*this->eta_);
      MVT::MvAddMv(alpha,*this->Adelta_,ONE,*this->Aeta_,*this->Aeta_);
      if (this->hasBOp_) {
        MVT::MvAddMv(alpha,*this->Bdelta_,ONE,*this->Beta_,*this->Beta_);
      }
      MVT::MvAddMv(alpha,*this->Hdelta_,ONE,*this->Heta_,*this->Heta_);
      // update its P-length
      e_Pe = e_Pe_new;

      // update gradient of m
      MVT::MvAddMv(alpha,*this->Hdelta_,ONE,*this->R_,*this->R_);
      {
        // re-tangentialize r
        Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
        this->orthman_->project(*this->R_,Teuchos::null,
                                Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                                Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->X_),
                                Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->BX_));
      }

      //
      // check convergence
      double r_norm = MAT::squareroot(SCT::real(ginner(*this->R_,*this->R_)));

      //
      // check local convergece 
      //
      // kappa (linear) convergence
      // theta (superlinear) convergence
      //
      double kconv = r0_norm * this->conv_kappa_;
      // insert some scaling here
      // double tconv = r0_norm * MAT::pow(r0_norm/normgradf0_,this->conv_theta_);
      double tconv = MAT::pow(r0_norm,this->conv_theta_+ONE);
      if ( r_norm <= ANASAZI_MIN(tconv,kconv) ) {
        if (tconv <= kconv) {
          innerStop_ = THETA_CONVERGENCE;
        }
        else {
          innerStop_ = KAPPA_CONVERGENCE;
        }
        break;
      }

      // z = Prec^-1 r
      // this must be in the tangent space
      zold_rold = z_r;
      if (this->hasPrec_) {
        Teuchos::TimeMonitor lcltimer( *this->timerPrec_ );
        OPT::Apply(*this->Prec_,*this->R_,*this->Z_);
        // the orthogonalization time counts under Ortho and under Preconditioning
        {
          Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
          this->orthman_->project(*this->Z_,Teuchos::null,
                                  Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                                  Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->X_),
                                  Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(this->BX_));
        }
        z_r = SCT::real( ginner(*this->R_,*this->Z_) );
      }
      else {
        z_r = SCT::real( ginner(*this->R_) );
      }

      // compute new search direction
      beta = z_r/zold_rold;
      MVT::MvAddMv(-ONE,*this->Z_,beta,*this->delta_,*this->delta_);
      // update new P-norms and P-dots
      e_Pd = beta*(e_Pd + alpha*d_Pd);
      d_Pd = z_r + beta*beta*d_Pd;

    } // end of the inner iteration loop

    if (this->om_->isVerbosity(Debug)) {
      this->om_->stream(Debug) 
        << " >> stop reason is " << stopReasons_[innerStop_] << endl
        << endl;
    }

  } // end of solveTRSubproblem

    
} // end of RBGen namespace


