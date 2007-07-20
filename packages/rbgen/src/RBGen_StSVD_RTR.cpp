#include "RBGen_StSVD_RTR.h"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"

namespace RBGen {

  StSVDRTR::StSVDRTR() : 
    isInitialized_(false),
    sigma_(0),
    tol_(1e-12),
    timerComp_("Total Elapsed Time"),
    debug_(false),
    verbLevel_(0),
    maxOuterIters_(500),
    maxInnerIters_(-1),
    localV_(true)
  {
    // set up array of stop reasons
    stopReasons_.push_back("n/a");
    stopReasons_.push_back("maximum iterations");
    stopReasons_.push_back("negative curvature");
    stopReasons_.push_back("exceeded TR");
    stopReasons_.push_back("kappa convergence");
    stopReasons_.push_back("theta convergence");
    stopReasons_.push_back("");
  }

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
    return sigma_;
  }

  void StSVDRTR::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                             const Teuchos::RCP< const Epetra_MultiVector >& ss,
                             const Teuchos::RCP< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio ) {

    using Teuchos::rcp;

    // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList rbmethod_params = params->sublist( "Reduced Basis Method" );

    // Get rank
    rank_ = rbmethod_params.get("Basis Size",(int)5);

    // Get convergence tolerance
    tol_ = rbmethod_params.get<double>("Convergence Tolerance",tol_);

    // Get debugging flag
    debug_ = rbmethod_params.get<bool>("StSVD Debug",debug_);

    // Get verbosity level
    verbLevel_ = rbmethod_params.get<int>("StSVD Verbosity Level",verbLevel_);

    // Get maximum number of outer iterations
    maxOuterIters_ = rbmethod_params.get<int>("StSVD Max Outer Iters",maxOuterIters_);

    // Get maximum number of inner iterations
    maxInnerIters_ = rbmethod_params.get<int>("StSVD Max Inner Iters",maxInnerIters_);

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
      TEST_FOR_EXCEPTION(ortho_ == Teuchos::null,std::invalid_argument,"User specified null ortho manager.");
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
    TEST_FOR_EXCEPTION(ss == Teuchos::null,std::invalid_argument,"Input snapshot matrix cannot be null.");
    A_ = ss;

    // initialize data structures
    initialize();
  }

  // private initialize
  void StSVDRTR::initialize() {
    if (isInitialized_) return;

    using Teuchos::rcp;

    // Allocate space for the factorization
    sigma_.resize( rank_ );
    U_ = Teuchos::null;
    V_ = Teuchos::null;
    U_ = rcp( new Epetra_MultiVector(A_->Map(),rank_,false) );
    if (localV_) {
      Epetra_LocalMap lclmap(A_->NumVectors(),0,A_->Comm());
      V_ = rcp( new Epetra_MultiVector(lclmap,rank_,false) );
    }
    else {
      Epetra_Map vmap(A_->NumVectors(),0,A_->Comm());
      V_ = rcp( new Epetra_MultiVector(vmap,rank_,false) );
    }
    resNorms_.resize(rank_);
    resUNorms_.resize(rank_);
    resVNorms_.resize(rank_);

    // allocate working multivectors
    RU_     = rcp( new Epetra_MultiVector(*U_) );
    AV_     = rcp( new Epetra_MultiVector(*U_) );
    etaU_   = rcp( new Epetra_MultiVector(*U_) );
    HeU_    = rcp( new Epetra_MultiVector(*U_) );
    deltaU_ = rcp( new Epetra_MultiVector(*U_) );
    HdU_    = rcp( new Epetra_MultiVector(*U_) );
    RV_     = rcp( new Epetra_MultiVector(*V_) );
    AU_     = rcp( new Epetra_MultiVector(*V_) );
    etaV_   = rcp( new Epetra_MultiVector(*V_) );     
    HeV_    = rcp( new Epetra_MultiVector(*V_) );
    deltaV_ = rcp( new Epetra_MultiVector(*V_) );
    HdV_    = rcp( new Epetra_MultiVector(*V_) );
    // allocate working space for DGESVD, UAVNsym, VAUNsym
    {
      Epetra_LocalMap lclmap(rank_,0,A_->Comm());
      dgesvd_A_ = rcp( new Epetra_MultiVector(lclmap,rank_,false) );
      TEST_FOR_EXCEPTION(dgesvd_A_->ConstantStride() == false,std::logic_error,
                         "RBGen::StSVD::initialize(): DGESVD Workspace A was not generated with constant stride!");
      UAVNsym_ = rcp( new Epetra_MultiVector(lclmap,rank_,false) );
      VAUNsym_ = rcp( new Epetra_MultiVector(lclmap,rank_,false) );
    }
    dgesvd_work_.resize(5*rank_);

    // Perform initial factorization
    // Make U,V orthonormal
    int ret = ortho_->normalize(*U_);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Initial U basis construction failed.");
    ret = ortho_->normalize(*V_);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Initial V basis construction failed.");
    // Compute A*V and A'*U
    {
      int info;
      info = AV_->Multiply('N','N',1.0,*A_,*V_,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }
    {
      int info;
      info = AU_->Multiply('T','N',1.0,*A_,*U_,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }

    // compute sigmas, f(x) and UAVNsym,VAUNsym from U,V,AU,AV
    updateF();

    // compute residuals
    updateResiduals();


    //
    // we are now initialized, albeit with null rank
    isInitialized_ = true;
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // reset with a new snapshot matrix
  void StSVDRTR::Reset( const Teuchos::RCP<Epetra_MultiVector>& new_ss ) {
    // Reset the pointer for the snapshot matrix
    // Note: We will not assume that it is non-null; user could be resetting our
    // pointer in order to delete the original snapshot set
    A_ = new_ss;
    isInitialized_ = false;
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // main loop
  void StSVDRTR::computeBasis() {

    using std::cout;
    using std::setprecision;
    using std::vector;
    using std::scientific;
    using std::setw;

    //
    // perform enough incremental updates to consume the entire snapshot set
    Teuchos::TimeMonitor lcltimer(timerComp_);

    // check that we have a valid snapshot set: user may have cleared 
    // it using Reset()
    TEST_FOR_EXCEPTION(A_ == Teuchos::null,std::logic_error,
                       "computeBasis() requires non-null snapshot set.");

    //
    // check that we are initialized, i.e., data structures match the data set
    initialize();

    // some bools
    bool tiny_rhonum, zero_rhoden, neg_rho;
    double rhonum, rhoden, fxnew;

    iter_ = 0;
    numInner_ = -1;
    innerStop_ = NOTHING;
    while (maxScaledNorm_ > tol_ && iter_ < maxOuterIters_) {

      ++iter_;

      // status printing
      if (verbLevel_ == 1) {
        // one line
        // acc TR+   k: %5d     num_inner: %5d     f: %e   |grad|: %e   stop_reason
        if (iter_) {
          cout << (accepted_ ? "accept" : "reject");
        }
        else {
          cout << "<init>";
        }
        cout << " " << tradjust_ 
             << "     k: " << setw(5) << iter_
             << "     num_inner: " << setw(5) << numInner_
             << "     f(x): " << scientific << setprecision(16) << fx_
             << "     |res|: " << scientific << setprecision(16) << maxScaledNorm_
             << stopReasons_[innerStop_] << endl;
      }
      else if (verbLevel_ > 1) {
        // multiline
        // 1:acc TR+   k: %5d     num_inner: %5d     stop_reason
        if (iter_) {
          cout << (accepted_ ? "accept" : "reject");
        }
        else {
          cout << "------";
        }
        cout << " " << tradjust_ 
             << "     k: " << setw(5) << iter_
             << "     num_inner: " << setw(5) << numInner_
             << stopReasons_[innerStop_] << endl;
        // 2:     f(x) : %e     |res| : %e
        cout << "     f(x) : " << scientific << setprecision(16) << fx_
             << "     |res|: " << scientific << setprecision(16) << maxScaledNorm_ << endl;
        // 3:    Delta : %e     |eta| : %e
        cout << "    Delta : " << scientific << setprecision(16) << Delta_
             << "    |eta| : " << scientific << setprecision(16) << etaLen_ << endl;
        if (neg_rho) {
          // 4:  NEGATIVE  rho     : %e
          cout << "  NEGATIVE  rho     : " << scientific << setprecision(16) << rho_ << endl;
        }
        else if (tiny_rhonum) {
          // 4: VERY SMALL rho_num : %e
          cout << " VERY SMALL rho_num : " << scientific << setprecision(16) << rhonum << endl;
        }
        else if (zero_rhoden) {
          // 4:    ZERO    rho_den : %e
          cout << "    ZERO    rho_den : " << scientific << setprecision(16) << rhoden << endl;
        }
        else {
          // 4:      rho : %e
          cout << "      rho : " << scientific << setprecision(16) << rho_ << endl;
        }
      }


      // compute gradient
      // this is the projected residuals
      Proj(*U_,*V_,*RU_,*RV_);

      // minimize model subproblem
      solveTRSubproblem();

      // save length of eta
      etaLen_ = SCT::squareroot(innerProduct(*etaU_,*etaV_));

      // compute proposed point
      // R_U(eta) = qf(U+eta)
      // we will store this in deltaU,deltaV
      {
        int ret;
        // U
        ret = deltaU_->Update(1.0,*U_,1.0,*etaU_,0.0);
        TEST_FOR_EXCEPTION(ret != 0,std::logic_error,
            "RBGen::StSVD::computeBasis(): error calling Epetra_MultiVector::Update.");
        ret = ortho_->normalize(*deltaU_);
        TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Retraction of etaU failed.");
        // V
        ret = deltaV_->Update(1.0,*V_,1.0,*etaV_,0.0);
        TEST_FOR_EXCEPTION(ret != 0,std::logic_error,
            "RBGen::StSVD::computeBasis(): error calling Epetra_MultiVector::Update.");
        ret = ortho_->normalize(*deltaV_);
        TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,"Retraction of etaV failed.");
      }

      //
      // evaluate rho
      //       f(x) - f(R_x(eta))
      // rho = -------------------
      //       m_x(eta) - m_x(eta)
      //
      //               f(x) - f(R_x(eta))
      //     = ------------------------------------
      //       - <eta,Proj(AV,A'U)> - .5 <eta,Heta>
      //
      //               f(x) - f(R_x(eta))
      //     = --------------------------------
      //       - <eta,(AV,A'U)> - .5 <eta,Heta>
      //
      // compute A'*newU into HdV_
      // we don't need A*newV yet
      {
        int info;
        info = HdV_->Multiply('T','N',1.0,*A_,*deltaU_,0.0);
        TEST_FOR_EXCEPTION(info != 0,std::logic_error,
            "RBGen::StSVD::computeBasis(): Error calling Epetra_MultiVector::Muliply.");
      }
      // compute fxnew and rhonum
      // we don't need sigmas, only trace(newU'*A*newV*N)
      // trace(newU'*A*newV*N) = sum_i (newU'*A*newV)_ii N[i]
      //                       = sum_i <newV[i],(A'newU)[i]> N[i]
      tiny_rhonum = neg_rho = zero_rhoden = false;
      {
        std::vector<double> dots(rank_);
        // deltaV_ stores newV, HdV_ stores newU'*A
        deltaV_->Dot(*HdV_,&dots[0]);
        fxnew = 0.0;
        for (int i=0; i<rank_; i++) {
          fxnew += dots[i]*N_[i];
        }
      }
      rhonum = fx_ - fxnew;
      // tiny rhonum means small decrease in f (maybe even negative small)
      // this usually happens near the end of convergence; 
      // this is usually not bad; we will pretend it is very good
      if ( abs(rhonum/fx_) < 1e-12 ) {
        tiny_rhonum = true;
        rho_ = 1.0;
      }
      else {
        // compute rhoden
        // - <eta,(AV,A'U)> - .5 <eta,Heta>
        rhoden = -1.0*innerProduct(*etaU_,*etaV_,*AV_,*AU_) 
                 -0.5*innerProduct(*etaU_,*etaV_,*HeU_,*HeV_);
        if (rhoden == 0) {
          // this is bad
          zero_rhoden = true;
          rho_ = -1.0;
        }
        else {
          rho_ = rhonum / rhoden;
        }
      }
      if (rho_ < 0.0) {
        // this is also bad
        neg_rho = true;
      }

      //
      // accept/reject
      if (rho_ >= rhoPrime_) {
        accepted_ = true;
        // put newx data into x data: U,V,AU,AV,fx,UAVNsym,VAUNsym,sigma

        // etaU,etaV get thrown away
        // deltaU_,deltaV_ go into U_,V_
        *U_ = *deltaU_;
        *V_ = *deltaV_;

        // HdV_ goes into AU_
        // A*V_ gets computed into AV_
        *AU_ = *HdV_;
        {
          int info;
          info = AV_->Multiply('N','N',1.0,*A_,*V_,0.0);
          TEST_FOR_EXCEPTION(info != 0,std::logic_error,
              "RBGen::StSVD::computeBasis(): Error calling Epetra_MultiVector::Muliply.");
        }

        // compute new sigmas, new f(x) and UAVNsym,VAUNsym from U,V,AU,AV
        updateF();

        if (debug_) {
          // check that new fx_ is equal to fxnew
          cout << " DBG: new f(x): " << fx_ << "\t\tnewfx: " << fxnew << endl;
        }
      }
      else {
        accepted_ = false;
      }

      //
      // modify trust-region radius
      tradjust_ = "   ";
      if (this->rho_ < .25) {
        Delta_ = .25 * Delta_;
        tradjust_ = "TR-";
      }
      else if (this->rho_ > .75 && (innerStop_ == NEGATIVE_CURVATURE || innerStop_ == EXCEEDED_TR)) {
        Delta_ = 2.0*Delta_;
        tradjust_ = "TR+";
      }

      //
      // compute new residuals and norms
      // this happens whether we accepted or not, because the trust-region solve
      // destroyed the residual that was in RU,RV
      updateResiduals();
    
    } // while (not converged)

    // nothing else to do.
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // update sigmas, f(x), UAVNsym, VAUNsym
  // this function assumes that U,V,AU,AV
  void StSVDRTR::updateF() {
    //
    // compute (U'*A)*V==RV'*V, compute its singular values
    {
      int info;
      info = dgesvd_A_->Multiply('T','N',1.0,*AU_,*V_,0.0);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }

    //
    // U'*A*V in dgesvd_A_ will be destroyed by GESVD; first, save it
    {
      // set UAVNsym_ = dgesvd_A_  * N
      // set VAUNsym_ = dgesvd_A_' * N
      for (int j=0; j<rank_; j++) {
        for (int i=0; i<rank_; i++) {
          (*UAVNsym_)[j][i] = (*dgesvd_A_)[j][i];
          (*VAUNsym_)[j][i] = (*dgesvd_A_)[i][j];
        }
        (*UAVNsym_)(j)->Scale(N_[j]);
        (*VAUNsym_)(j)->Scale(N_[j]);
      }
      // call Sym on both of these
      Sym(*UAVNsym_);
      Sym(*VAUNsym_);
    }

    //
    // compute f(U,V) = trace(U'*A*V*N) before destroying dgesvd_A_
    fx_ = 0.0;
    for (int j=0; j<rank_; j++) {
      fx_ += (*dgesvd_A_)[j][j]*N_[j];
    }

    //
    // compute the singular values of U'*A*V
    {
      int info;
      lapack.GESVD('N','N',rank_,rank_,dgesvd_A_->Values(),dgesvd_A_->Stride(),&sigma_[0],
                   NULL,0,NULL,0,&dgesvd_work_[0],dgesvd_work_.size(),NULL,&info);
      TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::StSVD::initialize(): Error calling Epetra_MultiVector::Muliply.");
    }
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // update the residuals and norms
  // this function assumes that U,V,AU,AV,sigma are all coherent
  void StSVDRTR::updateResiduals() {
    //
    // compute residuals: RU = A *V - U*S
    //                    RV = A'*U - V*S
    for (int i=0; i<rank_; i++) {
      (*RU_)(i)->Update( 1.0, *(*AV_)(i), sigma_[i], *(*U_)(i), 0.0 );
      (*RV_)(i)->Update( 1.0, *(*AU_)(i), sigma_[i], *(*V_)(i), 0.0 );
    }

    //
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
  std::vector<double> StSVDRTR::getResNorms() {
    if (isInitialized_) {
      return resNorms_;
    }
    else {
      return std::vector<double>(rank_,-1);
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

    int dim = (n_-rank_)*rank_ + (m_-rank_)*rank_ + (rank_-1)*rank_;
    if (maxInnerIters_ != -1 && maxInnerIters_ < dim) {
      dim = maxInnerIters_;
    }
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
    Proj(*U_,*V_,*RU_,*RV_);
    r_r = innerProduct(*RU_,*RV_);
    d_d = r_r;
    r0_norm = SCT::squareroot(r_r);
    // if first outer iteration, save r0_norm as normGrad0_
    if (iter_ == 0) {
      normGrad0_ = r0_norm;
    }

    // delta = -z
    deltaU_->Update(0.0,*RU_,-1.0);
    deltaV_->Update(0.0,*RV_,-1.0);
    e_d = 0.0;  // because eta = 0

    // the loop
    numInner_ = 0;
    for (int i=0; i<dim; ++i) {

      ++numInner_;

      // Hdelta = Hess*d 
      Hess(*U_,*V_,*deltaU_,*deltaV_,*HdU_,*HdV_);
      d_Hd = innerProduct(*deltaU_,*deltaV_,*HdU_,*HdV_);

      // compute update step
      alpha = r_r/d_Hd;

      if (debug_) {
        std::cout << " >> (r,r)  : " << r_r  << std::endl
                  << " >> (d,Hd) : " << d_Hd << std::endl
                  << " >> alpha  : " << alpha << std::endl;
      }

      // <neweta,neweta> = <eta,eta> + 2*alpha*<eta,delta> + alpha*alpha*<delta,delta>
      e_e_new = e_e + 2.0*alpha*e_d + alpha*alpha*d_d;

      // check truncation criteria: negative curvature or exceeded trust-region
      if (d_Hd <= 0 || e_e_new >= D2) {
        double tau = (-e_d + SCT::squareroot(e_d*e_d + d_d*(D2-e_e))) / d_d;
        if (debug_) {
          std::cout << " >> tau  : " << tau << std::endl;
        }
        // eta = eta + tau*delta
        etaU_->Update(tau,*deltaU_,1.0);
        etaV_->Update(tau,*deltaV_,1.0);
        HeU_->Update(tau,*HdU_,1.0);
        HeV_->Update(tau,*HdV_,1.0);
        if (d_Hd <= 0) {
          innerStop_ = NEGATIVE_CURVATURE;
        }
        else {
          innerStop_ = EXCEEDED_TR;
        }
        if (debug_) {
          e_e_new = e_e + 2.0*tau*e_d + tau*tau*d_d;
          std::cout << " >> predicted  (eta,eta) : " << e_e_new << std::endl
                    << " >> actual (eta,eta)     : " << innerProduct(*etaU_,*etaV_) << std::endl;
        }
        break;
      }

      // compute new eta == eta + alpha*delta
      etaU_->Update(alpha,*deltaU_,1.0);
      etaV_->Update(alpha,*deltaV_,1.0);
      HeU_->Update(alpha,*HdU_,1.0);
      HeV_->Update(alpha,*HdV_,1.0);
      // update its length
      e_e = e_e_new;

      // update gradient of model
      RU_->Update(alpha,*HdU_,1.0);
      RV_->Update(alpha,*HdV_,1.0);
      Proj(*U_,*V_,*RU_,*RV_);

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
      // some scaling of r0_norm so that it will be less than zero
      // and theta convergence may actually take over from kappa convergence
      double tconv = r0_norm * SCT::pow(r0_norm/normGrad0_,conv_theta_);
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
      deltaU_->Update(-1.0,*RU_,beta);
      deltaV_->Update(-1.0,*RV_,beta);
      // update new norms and dots
      e_d = beta*(e_d + alpha*d_d);
      d_d = r_r + beta*beta*d_d;

    } // end of the inner iteration loop

  } // end of solveTRSubproblem


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // routine for doing sym(S) in situ
  //
  // sym(S) is .5 (S + S')
  //
  void StSVDRTR::Sym( Epetra_MultiVector &S ) {
    // set strictly upper tri part of S to S+S'
    for (int j=0; j < rank_; j++) {
      for (int i=0; i < j; i++) {
        S[j][i] += S[i][j];
        S[j][i] /= 2.0;
      }
    }
    // set lower tri part of S to upper tri part of S
    for (int j=0; j < rank_-1; j++) {
      for (int i=j+1; i < rank_; i++) {
        S[j][i] = S[i][j];
      }
    }
  }


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
    // set S = .5 (S + S')
    Sym(S);
    // E = E - .5 U (S + S')
    etaU.Multiply(-1.0,xU,S,1.0);

    S.Multiply('T','N',1.0,xV,etaV,0.0);
    // set S = .5(S + S')
    Sym(S);
    // E = E - .5 V (S + S')
    etaV.Multiply(-1.0,xV,S,1.0);
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
    etaU.Norm2(&norms[0]);
    for (int i=0; i<rank_; i++) ip += norms[i]*norms[i];
    etaV.Norm2(&norms[0]);
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
    etaU.Dot(zetaU,&ips[0]);
    for (int i=0; i<rank_; i++) ip += ips[i];
    etaV.Dot(zetaV,&ips[0]);
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
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,
        "RBGen::StSVD::retract(): Ortho failure in retraction.");
    etaV.Update(1.0,xV,1.0);
    ret = ortho_->normalize(etaV);
    TEST_FOR_EXCEPTION(ret != rank_,std::runtime_error,
        "RBGen::StSVD::retract(): Ortho failure in retraction.");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Apply Hessian
  //
  // The Hessian can be applied as follows
  //   Hess f_{U,V}(eU,eV) = Proj [ A  eV N - eU sym(U' A  V N) ]
  //                              [ A' eU N - eV sym(V' A' U N) ]
  //                       = Proj [ A  eV N ] - eU sym(U' A  V N)
  //                              [ A' eU N ] - eV sym(V' A' U N)
  // We will apply the former as follows:
  // 1) Apply A eV into HeU and A' eU into HeV
  // 2) Scale HeU and HeV by N
  // 3) Compute sym(U' A V N) and sym(V' A' U N) from U' A V (
  //
  void StSVDRTR::Hess( 
      const Epetra_MultiVector &xU, 
      const Epetra_MultiVector &xV, 
      const Epetra_MultiVector &etaU, 
      const Epetra_MultiVector &etaV,
      Epetra_MultiVector &HetaU,
      Epetra_MultiVector &HetaV ) 
  {
    // HetaU = A etaV
    HetaU.Multiply('N','N',1.0,*A_,etaV,0.0);
    for (int i=0; i<rank_; i++) {
      HetaU(i)->Scale(N_[i]);
    }
    HetaU.Multiply('N','N',-1.0,etaU,*UAVNsym_,1.0);
    // HetaV = A' etaU
    HetaV.Multiply('T','N',1.0,*A_,etaU,0.0);
    for (int i=0; i<rank_; i++) {
      HetaV(i)->Scale(N_[i]);
    }
    HetaV.Multiply('N','N',-1.0,etaV,*VAUNsym_,1.0);
    // project
    Proj(xU,xV,HetaU,HetaV);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Debugging checks
  //
  // check all code invariants
  // this is gonna hurt.
  // 
  // where == 0: we are outside
  // check that 
  void StSVDRTR::Debug(int where) 
  {
  }

} // end of RBGen namespace
