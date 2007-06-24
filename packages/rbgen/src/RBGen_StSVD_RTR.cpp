#include "RBGen_IncSVDPOD.h"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

#define STSVD_DEBUG 1

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

    // Perform initial factorization
    // Make U,V orthonormal
    int ret = ortho_->normalize(*U_,Teuchos::null);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Initial U basis construction failed.");
    ret = ortho_->normalize(*U_,Teuchos::null);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Initial V basis construction failed.");
    // Compute initial singular values
    // finish: compute U'*A*V, compute its singular values
    // finish: retain parts that we need
    // finish: compute residuals: RU = A *V - U*S
    //                            RV = A'*U - V*S

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

    // CG terms
    double d_Hd, alpha, beta, r_r, rold_rold; 
    // terms for monitoring the norm of eta
    double e_e, e_d, d_d, e_e_new;
    // initial residual norm, used for stopping criteria
    double r0_norm;

    // set eta to zero
    etaU_->PutScalar(0.0);
    etaV_->PutScalar(0.0);
    HeU_->PutScalar(0.0);
    HeV_->PutScalar(0.0);
    e_e = 0.0;

    //
    // We have two residuals: RU = A *V - U*S
    //                    and RV = A'*U - V*S
    //
    // Note that grad(U,V) = {proj(A *V*N)} = {proj(RU*N)}
    //                       {proj(A'*U*N)} = {proj(RV*N)}
    // 
    // Project the residuals in RU_ and RV_ to yield the gradient
    //
    // r0 = grad f(X) = Proj R_
    // We will do this in place.
    //
    Proj(U_,V_,RU_,RV_);
    r_r = innerProduct(*RU_,*RV_);
    d_d = r_r;
    r0_norm = SCT::squareroot(r_r);

    // delta = -z
    deltaU_->Update(0.0,*RU_,-1.0);
    deltaV_->Update(0.0,*RV_,-1.0);
    e_d = 0.0;  // because eta = 0

    // the loop
    for (int i=0; i<d; ++i) {

      // Hdelta = Hess*d 
      Hess(*U_,*V_,*deltaU_,*deltaV_,*HdU_,*HdV_);
      d_Hd = ginner(*deltaU_,*deltaV_,*HdU_,*HdV_);

      // compute update step
      alpha = r_r/d_Hd;

#ifdef STSVD_DEBUG
        cout << " >> (r,r)  : " << r_r  << endl
             << " >> (d,Hd) : " << d_Hd << endl
             << " >> alpha  : " << alpha << endl;
#endif

      // <neweta,neweta> = <eta,eta> + 2*alpha*<eta,delta> + alpha*alpha*<delta,delta>
      e_e_new = e_e + 2.0*alpha*e_d + alpha*alpha*d_d;

      // check truncation criteria: negative curvature or exceeded trust-region
      if (d_Hd <= 0 || e_e_new >= D2) {
        double tau = (-e_d + SCT::squareroot(e_d*e_d + d_d*(D2-e_e))) / d_d;
#ifdef STSVD_DEBUG
        cout << " >> tau  : " << tau << endl;
#endif
        // eta = eta + tau*delta
        etaU_->Update(1.0,tau,*deltaU_);
        etaV_->Update(1.0,tau,*deltaV_);
        HeU_->Update(1.0,tau,*HdU_);
        HeV_->Update(1.0,tau,*HdV_);
        if (d_Hd <= 0) {
          innerStop_ = NEGATIVE_CURVATURE;
        }
        else {
          innerStop_ = EXCEEDED_TR;
        }
#ifdef STSVD_DEBUG
        e_e_new = e_e + 2.0*tau*e_d + tau*tau*d_d;
        cout << " >> predicted  (eta,eta) : " << e_e_new << endl
             << " >> actual (eta,eta)     : " << innerProduct(*etaU_,*etaV_) << endl;
#endif
        break;
      }

      // compute new eta == eta + alpha*delta
      etaU_->Update(1.0,alpha,*deltaU_);
      etaV_->Update(1.0,alpha,*deltaV_);
      HeU_->Update(1.0,alpha,*HdU_);
      HeV_->Update(1.0,alpha,*HdV_);
      // update its length
      e_e = e_e_new;

      // update gradient of m
      RU_->Update(1.0,alpha,*HdU_);
      RV_->Update(1.0,alpha,*HdV_);
      Project(*U_,*V_,*RU_,*RV_);

      rold_rold = r_r;
      r_r = innerProduct(*RU_,*RV_);
      double r_norm = SCT::squareroot(r_r);

      //
      // check local convergece 
      //
      // kappa (linear) convergence
      // theta (superlinear) convergence
      //
      double kconv = r0_norm * conv_kappa_;
      // insert some scaling here
      // double tconv = r0_norm * SCT::pow(r0_norm/normgradf0_,this->conv_theta_);
      double tconv = SCT::pow(r0_norm,conv_theta_+1.0);
      double conv = kconv < tconv ? kconv : tconv;
      if ( r_norm <= conv ) {
        if (conv == tconv) {
          innerStop_ = THETA_CONVERGENCE;
        }
        else {
          innerStop_ = KAPPA_CONVERGENCE;
        }
        break;
      }

      // compute new search direction
      beta = r_r/rold_rold;
      deltaU_->Update(beta,-1.0,*RU_);
      deltaV_->Update(beta,-1.0,*RV_);
      // update new norms and dots
      e_d = beta*(e_d + alpha*d_d);
      d_d = r_r + beta*beta*d_d;

    } // end of the inner iteration loop

#ifdef STSVD_DEBUG
    cout << " >> stop reason is " << stopReasons_[innerStop_] << endl
          << endl;
#endif

  } // end of solveTRSubproblem

    
} // end of RBGen namespace
