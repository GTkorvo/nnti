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
    localV_(true)
  {}

  Teuchos::RCP<const Epetra_MultiVector> StSVDRTR::getBasis() const {
    if (isInitialized_ == false) {
      return Teuchos::null;
    }
    return U_;
  }

  Teuchos::RCP<const Epetra_MultiVector> StSVDRTR::getRightBasis() const {
    if (isInitialized_ == false) {
      return Teuchos::null;
    }
    return V_;
  }

  std::vector<double> StSVDRTR::getSingularValues() const { 
    return ret;
  }

  void StSVDRTR::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                             const Teuchos::RCP< Epetra_MultiVector >& ss,
                             const Teuchos::RCP< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio ) {

    using Teuchos::rcp;

    // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList rbmethod_params = params->sublist( "Reduced Basis Method" );

    // Get rank
    int rank_ = rbmethod_params.get("Rank",(int)5);

    // Get convergence tolerance
    tol_ = rbmethod_params.get<double>("Convergence Tolerance",tol_);

    // Get debugging flag
    debug_ = rbmethod_params.get<bool>("StSVD Debug",debug_);

    // Get verbosity level
    verbLevel_ = rbmethod_params.get<int>("StSVD Verbosity Level",verbLevel_);

    // Is V local or distributed
    localV_ = rbmethod_params.get<bool>("StSVD Local V",localV_);

    // Get an Anasazi orthomanager
    if (rbmethod_params.isType<
          Teuchos::RCP< Anasazi::OrthoManager<double,Epetra_MultiVector> > 
        >("Ortho Manager")
       ) 
    {
      ortho_ = rbmethod_params.get< 
                Teuchos::RCP<Anasazi::OrthoManager<double,Epetra_MultiVector> >
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

    // initialize data structures
    initialize();
  }

  // private initialize
  void StSVDRTR::initialize() {
    if (isInitialized_) return;

    // Allocate space for the factorization
    sigma_.resize( rank_ );
    U_ = Teuchos::null;
    V_ = Teuchos::null;
    U_ = rcp( new Epetra_MultiVector(ss->Map(),rank_,false) );
    if (localV_) {
      Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
      V_ = rcp( new Epetra_MultiVector(lclmap,rank_,false) );
    }
    else {
      Epetra_Map vmap(ss->NumVectors(),0,ss->Comm());
      V_ = rcp( new Epetra_MultiVector(vmap,rank_,false) );
    }
    resNorms_.resize(rank_);
    resUNorms_.resize(rank_);
    resVNorms_.resize(rank_);

    // allocate working multivectors
    RU_     = rcp( new Epetra_MultiVector(*U_) );
    etaU    = rcp( new Epetra_MultiVector(*U_) );_
    HeU_    = rcp( new Epetra_MultiVector(*U_) );
    deltaU_ = rcp( new Epetra_MultiVector(*U_) );
    HdU_    = rcp( new Epetra_MultiVector(*U_) );
    RV_     = rcp( new Epetra_MultiVector(*V_) );
    etaV_   = rcp( new Epetra_MultiVector(*V_) );     
    HeV_    = rcp( new Epetra_MultiVector(*V_) );
    deltaV_ = rcp( new Epetra_MultiVector(*V_) );
    HdV_    = rcp( new Epetra_MultiVector(*V_) );
    // allocate working space for DGESVD
    {
      Epetra_LocalMap lclmap(rank_,0,ss->Comm());
      dgesvd_A_ = rcp( new Epetra_MultiVector(lclap,rank_,false) );
      TEST_FOR_EXCEPTION(dgesvd_A_->ConstantStride() == false,std::logic_error,
                         "RBGen::StSVD::initialize(): DGESVD Workspace A was not generated with constant stride!");
    }
    dgesvd_work_.resize(5*rank_);

    // Perform initial factorization
    // Make U,V orthonormal
    int ret = ortho_->normalize(*U_);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Initial U basis construction failed.");
    ret = ortho_->normalize(*V_);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Initial V basis construction failed.");
    // Compute initial singular values
    // compute A*V and A'*U
    {
      int info;
      info = RU_->Multiply('N','N',1.0,*ss,*V_,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }
    {
      int info;
      info = RV_->Multiply('T','N',1.0,*ss,*U_,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }
    // compute (U'*A)*V==RV'*V, compute its singular values
    {
      int info;
      info = dgesvd_A_->Multiply('T','N',1.0,*RV_,*V_,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }
    {
      int info;
      lapack.GESVD('N','N',rank_,rank_,dgesvd_A_->Values(),dgesvd_A_->Stride(),&sigma_[0],
                   NULL,0,NULL,0,&dgesvd_work[0],dgesvd_lwork_.size(),&info);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }
    // compute residuals: RU = A *V - U*S
    //                    RV = A'*U - V*S
    for (int i=0; i<rank_; i++) {
      (*RU_)(i)->Update( *U(i), sigma_[i], 1.0 );
      (*RV_)(i)->Update( *V(i), sigma_[i], 1.0 );
    }
    // compute residual norms
    // for each singular value, we have two residuals:
    //   ru_i = A v_i - u_i s_i 
    // and 
    //   rv_i = A' u_i - v_i s_i 
    // If we have a singular value, these will both be zero
    // We will take the norm of the i-th residual to be
    //   |r_i| = sqrt(|ru_i|^2 + |rv_i|^2)
    // We will scale it by the i-th singular value
    RU_->Norm2(&resUNorms_[0]);
    RV_->Norm2(&resVNorms_[0]);
    for (int i=0; i<rank_; i++) {
      resNorms_[i] = std::sqrt(resUNorms_[i]*resUNorms_[i] + resVNorms_[i]*resVNorms_[i]);
      // sigma_ must be >= 0, because it is a singular value
      // check sigma_ > 0 instead of sigma_ != 0, so the compiler won't complain
      maxScaledNorm_ = 0;
      if (sigma_[i] > 0.0) {
        resNorms_[i] /= sigma_[i];
        maxScaledNorm_ = EPETRA_MAX(resNorms_[i],maxScaledNorm_);
      }
    }

    // we are now initialized, albeit with null rank
    isInitialized_ = true;
  }

  void StSVDRTR::Reset( const Teuchos::RCP<Epetra_MultiVector>& new_ss ) {
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

    //
    // check that we are initialized, i.e., data structures match the data set
    initialize();

    //
    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    while (maxScaledNorm_ > tol_) {

      // status printing: finish

      // compute gradient
      // this is the projected residuals
      Proj(*U_,*V_,*RU_,*RV_);

      // minimize model subproblem
      solveTRSubproblem();

      // compute proposed point
      // R_U(eta) = qf(U+eta)
      // we will store this in eta
      {
        int ret;
        // U
        ret = etaU_->Update(1.0,*U_,1.0);
        TEST_FOR_EXCEPTION(ret != 0,std::logic_error,
            "RBGen::StSVD::computeBasis(): error calling Epetra_MultiVector::Update.");
        ret = ortho_->normalize(*etaU_);
        TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Retraction of etaU failed.");
        // V
        ret = etaV_->Update(1.0,*V_,1.0);
        TEST_FOR_EXCEPTION(ret != 0,std::logic_error,
            "RBGen::StSVD::computeBasis(): error calling Epetra_MultiVector::Update.");
        ret = ortho_->normalize(*etaV_);
        TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Retraction of etaV failed.");
      }

      // evaluate rho: finish

      // accept/reject and adjust trust-region radius: finish

      // compute new residuals and norms: finish

    } // while (not converged)

    // finish
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // update the basis with new snapshots: not supported
  void StSVDRTR::updateBasis( const Teuchos::RCP< Epetra_MultiVector >& update_ss ) {
    // perform enough incremental updates to consume the new snapshots
    TEST_FOR_EXCEPTION(true,std::logic_error,
        "RBGen::StSVDRTR::updateBasis(): this routine not yet supported.");
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return residual norms
  const std::vector<double> & StSVDRTR::getResNorms() {
    if (isInitialized_) {
      return resNorms_;
    }
    else {
      return std::vector(rank_,nan);
    }
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


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Projection onto tangent plane
  void StSVDRTR::Proj( const Epetra_MultiVector &xU, 
                       const Epetra_MultiVector &xV, 
                       Epetra_MultiVector &etaU, 
                       Epetra_MultiVector &etaV ) {
    // perform etaU = Proj_U etaU
    // and     etaV = Proj_V etaV
    // Proj_U E = (I - U U') E + U skew(U' E)
    //          = E - U (U' E) + .5 U (U' E - E' U)
    //          = E - .5 U (S + S')
    // where S = U' E
    static Epetra_LocalMap lclmap(rank_,0,A_->Comm());
    static Epetra_MultiVector S(lclmap,rank_);
    S.Multiply('T','N',1.0,xU,etaU,0.0);
    // set upper tri part of S to S+S'
    for (int j=0; j < rank_; j++) {
      for (int i=0; i <= j; i++) {
        S[j][i] += S[i][j];
      }
    }
    // set lower tri part of S to upper tri part of S
    for (int j=0; j < rank_-1; j++) {
      for (int i=j+1; i < rank_; i++) {
        S[j][i] = S[i][j];
      }
    }
    // E = E - .5 U S
    etaU.Multiply(-0.5,xU,E,1.0);

    S.Multiply('T','N',1.0,xV,etaV,0.0);
    // set upper tri part of S to S+S'
    for (int j=0; j < rank_; j++) {
      for (int i=0; i <= j; i++) {
        S[j][i] += S[i][j];
      }
    }
    // set lower tri part of S to upper tri part of S
    for (int j=0; j < rank_-1; j++) {
      for (int i=j+1; i < rank_; i++) {
        S[j][i] = S[i][j];
      }
    }
    // E = E - .5 V S
    etaV.Multiply(-0.5,xV,E,1.0);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // g(eta,eta)
  double StSVDRTR::innerProduct( 
      const Epetra_MultiVector &etaU, 
      const Epetra_MultiVector &etaV ) 
  {
    // g( (etaU,etaV), (etaU,etaV) ) = trace(etaU'*etaU) + trace(etaV'*etaV)
    std::vector<double> norms(rank_);
    double ip = 0.0;
    etaU.Norm2(&norms);
    for (int i=0; i<rank_; i++) ip += norms[i]*norms[i];
    etaV.Norm2(&norms);
    for (int i=0; i<rank_; i++) ip += norms[i]*norms[i];
    return ip;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // g(eta,zeta)
  double StSVDRTR::innerProduct( 
      const Epetra_MultiVector &etaU, 
      const Epetra_MultiVector &etaV, 
      const Epetra_MultiVector &zetaU, 
      const Epetra_MultiVector &zetaV ) 
  {
    // g( (etaU,etaV), (etaU,etaV) ) = trace(etaU'*etaU) + trace(etaV'*etaV)
    std::vector<double> ips(rank_);
    double ip = 0.0;
    etaU.Dot(zetaU,&ips);
    for (int i=0; i<rank_; i++) ip += ips[i];
    etaV.Dot(zetaV,&ips);
    for (int i=0; i<rank_; i++) ip += ips[i];
    return ip;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Retraction, in situ
  void StSVDRTR::retract( 
      const Epetra_MultiVector &xU, 
      const Epetra_MultiVector &xV, 
      Epetra_MultiVector &etaU, 
      Epetra_MultiVector &etaV ) 
  {
    int ret;
    // R_U(E) = qf(U+E)
    etaU.Update(1.0,xU,1.0);
    ret = ortho_->normalize(etaU);
    TEST_FOR_EXCEPTION(ret != rank_,std::exception,
        "RBGen::StSVD::retract(): Ortho failure in retraction.");
    etaV.Update(1.0,xV,1.0);
    ret = ortho_->normalize(etaV);
    TEST_FOR_EXCEPTION(ret != rank_,std::exception,
        "RBGen::StSVD::retract(): Ortho failure in retraction.");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Apply Hessian
  void StSVDRTR::Hess( 
      const Epetra_MultiVector &xU, 
      const Epetra_MultiVector &xV, 
      const Epetra_MultiVector &etaU, 
      const Epetra_MultiVector &etaV,
      Epetra_MultiVector &HetaU,
      Epetra_MultiVector &HetaV ) 
  {
    // finish
  }


} // end of RBGen namespace
