
#include "Aristos_SQPAlgoDecl.hpp"

namespace Aristos {

typedef Teuchos::RefCountPtr<Aristos::Vector> VectorPtr;


SQPAlgo::SQPAlgo( DataPool &dat, Objective &obj, Constraints &constr, HessVec &hessvec,
                  LagMult  &lagmult, FeasStep &feasstep )
{
  dat_      = Teuchos::rcp(&dat, false);
  obj_      = Teuchos::rcp(&obj, false);
  constr_   = Teuchos::rcp(&constr, false);
  hessvec_  = Teuchos::rcp(&hessvec, false);
  lagmult_  = Teuchos::rcp(&lagmult, false);
  feasstep_ = Teuchos::rcp(&feasstep, false);
}

  
bool SQPAlgo::run(Vector &x, Vector &c, Vector &l, int &iter, int &iflag, Teuchos::ParameterList &parlist){
  // SQP algorithm parameters.
  double Delta = parlist.get("Initial TR Radius", 1.0e2);
  double nu    = parlist.get("Merit Function Penalty", 1.0);
  double zeta  = parlist.get("TR Radius Scaling", 0.8);

  // Stopping tolerances. Termination if:
  // ( norm(fgrad - Jac'*lambda, inf) < gtol AND
  //    norm(c, inf) < ctol) ) OR Delta < stol
  double gtol = parlist.get("Gradient of Lagrangian Tolerance", 1.0e-6);
  double ctol = parlist.get("Constraints Tolerance", 1.0e-6);
  double stol = parlist.get("Min TR Radius", 1.0e-5);

  // Maximum number of outer SQP iterations.
  int maxiter = parlist.get("Max Number of SQP Iterations", 100);

  // Nominal fixed tolerance for inexact solvers. 0 means almost exact computations.
  double tol = 0.0;

  iter = 0;         // Set SQP iteration counter.



  // Vector variables.
  VectorPtr g  = x.createVector();   // gradient of objective
  VectorPtr gl = x.createVector();   // gradient of Lagrangian
  VectorPtr v  = x.createVector();   // quasi-normal step
  VectorPtr w  = x.createVector();   // tangential step


  printf("\n SQP \n");
  printf("    k    f(x_k)     || c(x_k) || ||L_x(x_k,lambda_k)||   Delta_k      ");
  printf("|| n_k ||     || t_k ||    itcg  accept \n");

  // Set first iterate and compute DataPool quantities.
  dat_->computeAll(x);

  // Compute function and constraint values at initial iterate.
  double f = obj_->getValue(x);
  obj_->getGradient(x, *g);
  constr_->getValue(x, c);
  // Compute Lagrange multipliers.
  runMultiplierEst(x, *g, l, tol);    

  constr_->applyJacobian(true, x, l, *gl);
  gl->linComb(1.0, *g, -1.0);
  double normgl = sqrt(gl->innerProd(*gl));
  double normc  = sqrt(c.innerProd(c));

  printf("%5d %12.5e  %12.6e     %12.6e     %12.6e \n", iter, f, normc, normgl, Delta);


  while (   (iter < maxiter) && ((normgl >= gtol) || (normc >= gtol)) && (Delta >= stol) ) {

    int iaccept = 0;  // runAcceptStep return flag
    int itcg = 0;     // number of CG iterations in runTangentialStep

    while ( (iaccept == 0) && (Delta >= stol) ) {
         // Compute the quasi-normal step.
         runQuasiNormalStep(x, c, *v, zeta*Delta, tol);

         // Compute the tangential step.
         int flg;                   // return flag for tangential step computation
         int cgitmax  = 100;        // set maximum number of CG iterations
         double cgtol = 1.e-2;      // set CG relative tolerance
         runTangentialStepInx(x, *gl, *v, l, *w, Delta, cgitmax, cgtol, tol, itcg, flg);

         // Check step, adjust trust-region radius and parameter.
         runAcceptStep(x, l, f, *g, c, *v, *w, Delta, nu, tol, iaccept);

        if ( iaccept == 0 ) {
            // Reset DataPool to old values.
            dat_->computeAll(x);
            printf("%5d %12.5e  %12.6e     %12.6e     %12.6e  %12.6e  %12.6e    %d    %d\n",
                    iter, f, normc, normgl, Delta, sqrt(v->innerProd(*v)), sqrt(w->innerProd(*w)), itcg, iaccept);
        }
    }   // End inner while loop.

    if ( iaccept > 0 ) {
        iter++;

        // Compute function and constraint values at new iterate.
        f = obj_->getValue(x);
        obj_->getGradient(x, *g);
        constr_->getValue(x, c);
        // Compute Lagrange multipliers.
        runMultiplierEst(x, *g, l, tol);    

        constr_->applyJacobian(true, x, l, *gl);
        gl->linComb(1.0, *g, -1.0);
        normgl = sqrt(gl->innerProd(*gl));
        normc  = sqrt(c.innerProd(c));

        printf("%5d %12.5e  %12.6e     %12.6e     %12.6e  %12.6e  %12.6e    %d    %d\n",
                iter, f, normc, normgl, Delta, sqrt(v->innerProd(*v)), sqrt(w->innerProd(*w)), itcg, iaccept);
    }

  } // End outer while loop.

  // Set termination flag.
  if (iter > maxiter)
    iflag = 1;
  else if (Delta < stol)
    iflag = 2;
  else
    iflag = 0;

  return true;
}


bool SQPAlgo::runMultiplierEst(const Vector &x, const Vector &g, Vector &l, double tol)
{
  lagmult_->getValue(x, g, l, tol);
  return true;
}


bool SQPAlgo::runQuasiNormalStep(const Vector &x, const Vector &c, Vector &s, double delta, double tol)
{
  // Create necessary vectors.
  VectorPtr stmp = c.createVector();
  VectorPtr sCP  = x.createVector();
  VectorPtr sN   = x.createVector();

  // Compute Cauchy step
  constr_->applyJacobian(true, x, c, s);                  // s = J(x)'*c
  constr_->applyJacobian(false, x, s, *stmp);             // stmp = J(x)*s
  sCP->Set(-s.innerProd(s)/stmp->innerProd(*stmp), s);    // sCP = - (||s||^2/||stmp||^2) * s

  // If the Cauchy step is outside the trust region, return the
  // scaled Cauchy step.

  double normsCP = sqrt(sCP->innerProd(*sCP));
  if (normsCP >= delta) {
    s.Set(delta/normsCP, *sCP);
    return true;
  }

  // Compute feasibility step, for example, by solving a problem related
  // to finding the the minimum norm solution of   min || Jac*s + c ||.
  feasstep_->getValue(x, c, *sN, tol);
  double normsN = sqrt(sN->innerProd(*sN));

  if (normsN <= delta) {
    // take full feasibility step
    s.Set(1.0, *sN);
  }
  else {
    // take convex combination  s = sCP+tau*(sN-sCP)
    // so that || s || = delta

    // solve scalar quadratic equation: || sCP+tau*(sN-sCP)||_2^2 = delta^2
    s.Set(1.0, *sN);
    s.linComb(-1.0, *sCP, 1.0);
    double a   = s.innerProd(s);
    double b   = s.innerProd(*sCP);
    double c   = normsCP*normsCP - delta*delta;
    double tau = (-b+sqrt(b*b-a*c))/a;

    s.linComb(1.0, *sCP, tau);
  }
  
  return true;

}


bool SQPAlgo::runTangentialStep(const Vector &x, const Vector &g, const Vector &v, const Vector &l, Vector &s,
                                double delta, int cgitmax, double cgtol, double tol, int &cgiter, int &iflag)
{
  bool iprint = true;  // print flag  = true print output;
                        // otherwise no output is generated

  // Create necessary vectors.
  VectorPtr r    = x.createVector();
  VectorPtr z    = x.createVector();
  VectorPtr p    = x.createVector();
  VectorPtr Hp   = x.createVector();
  VectorPtr stmp = x.createVector();


  // Initialization of the CG step.
  double nullsptol[1];
         nullsptol[0] = tol;
  cgiter = 0;
  s.Set(0.0);
  hessvec_->getValue(x, l, v, *r);
  r->linComb(1.0, g, 1.0);

  // Compute projection z=P*r
  constr_->applyNullSp(false, x, *r, *z, nullsptol);

  p->Set(-1.0, *z);
  double rz   = r->innerProd(*z);
  double rz0  = rz;

  if (iprint) {
    printf("\n SQP_horizontal_step \n");
    printf(" iter   || P*r||/|| P*r0||   || D*s ||      delta  \n");
    printf("%5d     %12.5e     %12.5e  %12.5e  \n", 0, 1.0, 0.0, delta);
  }

  // Start CG loop.
  for (cgiter=1; cgiter<=cgitmax; cgiter++) { 

    hessvec_->getValue(x, l, *p, *Hp);

    double pHp = p->innerProd(*Hp);
    if (pHp <= 0.0) {
        // p is a direction of negative/zero curvature.
        // compute theta > 0 such that || s + theta*p || = delta
        
        iflag = 1;
        double a     = p->innerProd(*p);
        double b     = s.innerProd(*p);
        double theta = (-b + sqrt(b*b - a*(s.innerProd(s)-delta*delta))) / a;
        s.linComb(theta, *p, 1.0);
        if (iprint)
           printf(" negative curvature detected \n");
        return true;
    }
  
    double alpha = rz/pHp;

    stmp->Set(1.0, s);
    stmp->linComb(alpha, *p, 1.0);
    double norms = sqrt(stmp->innerProd(*stmp));
    if (norms > delta) {
        // compute theta > 0 such that || s + theta*p || = delta
        
        iflag = 2;
        double a     = p->innerProd(*p);
        double b     = s.innerProd(*p);
        double theta = (-b + sqrt(b*b - a*(s.innerProd(s)-delta*delta))) / a;
        s.linComb(theta, *p, 1.0);
        if (iprint)
           printf(" TR-radius active \n");
        return true;
    }

    // update step and iterate
    s.Set(1.0, *stmp);
    r->linComb(alpha, *Hp, 1.0);

    // Compute projection z =P*r
    constr_->applyNullSp(false, x, *r, *z, nullsptol);

    double rz_new = r->innerProd(*z);

    if (iprint)
       printf("%5d     %12.5e     %12.5e  %12.5e  \n", cgiter, rz_new/rz0, norms, delta);

    if (rz_new <= (cgtol*cgtol)*rz0) {
        iflag = 0;
        if (iprint)
           printf(" || P(g + H*(n+s)) || <= cgtol*|| P(g + H*n)|| \n");
        return true;
    }

    double beta = rz_new/rz;
    p->linComb(-1.0, *z, beta);
    rz     = rz_new;
 
  } // End CG loop.

  iflag = 3;
  if (iprint)
    printf(" maximum number of iterations reached \n");
  
  return true;

}  // End runTangentialStep.


bool SQPAlgo::runTangentialStepInx(const Vector &x, const Vector &g, const Vector &v, const Vector &l, Vector &s,
                                   double delta, int cgitmax, double cgtol, double tol, int &cgiter, int &iflag)
{
  bool iprint = false;  // print flag  = true print output;
                        // otherwise no output is generated
  bool ifixed = false;  // inner solver tolerance flag: if == true use fixed tolerance
                        //                              otherwise use SQP dynamic criteria

  // Create necessary vectors.
  VectorPtr r    = x.createVector();
  VectorPtr z    = x.createVector();
  vector<VectorPtr> pvecs;
  //pvecs.push_back(x.createVector());
  VectorPtr p    = x.createVector();
  VectorPtr d    = x.createVector();
  VectorPtr Hp   = x.createVector();
  VectorPtr stmp = x.createVector();


  // Initialization of the CG step.
  cgiter = 0;
  s.Set(0.0);
  d->Set(0.0);
  hessvec_->getValue(x, l, v, *r);
  r->linComb(1.0, g, 1.0);

  double normg  = sqrt(r->innerProd(*r)); // norm of g = (grad of Lagr) + Hessian*(quasi-normal step)
  double normWr[cgitmax+1];     // will keep track of norms of projected residuals
         normWr[0] = 1e0;
  double normWg = 0.0;          // norm of inexact reduced gradient
  double smallc = 1e-2*cgtol;   // small heuristic constant
  double fixedtol = 1e-12;      // fixed tolerance, used if ifixed==true
  double vartols[4];

  if (iprint) {
    printf("\n SQP_horizontal_step \n");
    printf(" iter   || W*r||/|| W*r0||   || D*s ||      delta  \n");
  }


  // Start CG loop.
  for (cgiter=1; cgiter<=cgitmax; cgiter++) { 

    // Compute (inexact) projection z=W*r
    if (ifixed) {
      vartols[0] = fixedtol;
      vartols[1] = 0.0;     // if 1.0 solver will use relative residuals, if 0.0 absolute
      constr_->applyNullSp(false, x, *r, *z, vartols);
    }
    else {
      if (cgiter == 1) {
        vartols[0] = -1.0; // will be determined later
        vartols[1] = normg; vartols[2] = delta; vartols[3] = 1e-1*smallc;
        constr_->applyNullSp(false, x, *r, *z, vartols);
        normWg = sqrt(z->innerProd(*z));
      }
      else {
        vartols[0] = min(normWg/normg, delta/normg);
        vartols[0] = min(vartols[0], smallc);
        vartols[0] = vartols[0] * min(normWr[cgiter-1], 1.0);
        vartols[1] = 0.0;     // solver will use absolute residuals
        constr_->applyNullSp(false, x, *r, *z, vartols);
      }
    }
    normWr[cgiter] = sqrt(z->innerProd(*z));

    if (iprint)
       printf("%5d     %12.5e     %12.5e  %12.5e  \n",
              cgiter-1, normWr[cgiter]/normWr[1], sqrt(s.innerProd(s)), delta);

    // Check if done.
    if (normWr[cgiter]/normWr[1] < cgtol) {
        iflag = 0;
        cgiter = cgiter-1;
        if (iprint)
           printf(" || W(g + H*(n+s)) || <= cgtol*|| W(g + H*n)|| \n");
        return true;
    }

    // orthogonalize
    d->Set(0.0);
    for (int j=0; j < cgiter-1; j++) {
    //for (int j=max(0,cgiter-2); j < cgiter-1; j++)  {  // use this to simulate true CG
      hessvec_->getValue(x, l, *pvecs[j], *Hp);
      d->linComb( (z->innerProd(*Hp)) / ((pvecs[j])->innerProd(*Hp)), *pvecs[j], 1.0);
    }
    pvecs.push_back(x.createVector());
    (pvecs[cgiter-1])->Set(-1.0, *z);
    (pvecs[cgiter-1])->linComb(1.0, *d, 1.0);
    // set temporary p, simplifies things later
    p = pvecs[cgiter-1];
    hessvec_->getValue(x, l, *p, *Hp);

    ////////////////////////////////// negative curvature stopping condition
    double pHp = p->innerProd(*Hp);
    if (pHp <= 0.0) {
        // p is a direction of negative/zero curvature.
        // compute theta > 0 such that || s + theta*p || = delta
        
        iflag = 1;
        double a     = p->innerProd(*p);
        double b     = s.innerProd(*p);
        double theta = (-b + sqrt(b*b - a*(s.innerProd(s)-delta*delta))) / a;
        s.linComb(theta, *p, 1.0);
        if (iprint)
           printf(" negative curvature detected \n");
        return true;
    }

    /////////////////////////////////////// want to enforce positive alpha's
    if (r->innerProd(*p) >= 0.0) {
        iflag = 4;
        double a     = p->innerProd(*p);
        double b     = s.innerProd(*p);
        double theta = (-b + sqrt(b*b - a*(s.innerProd(s)-delta*delta))) / a;
        s.linComb(theta, *p, 1.0);
        if (iprint)
           printf(" negative search direction coefficient \n");
        return true;
    }
  
    double alpha = - (r->innerProd(*p))/pHp;

    ///////////////////////////////////////////////////////// iterate update
    stmp->Set(1.0, s);
    stmp->linComb(alpha, *p, 1.0);
    double norms = sqrt(stmp->innerProd(*stmp));
    if (norms > delta) {
        // compute theta > 0 such that || s + theta*p || = delta
        
        iflag = 2;
        double a     = p->innerProd(*p);
        double b     = s.innerProd(*p);
        double theta = (-b + sqrt(b*b - a*(s.innerProd(s)-delta*delta))) / a;
        s.linComb(theta, *p, 1.0);
        if (iprint)
           printf(" TR-radius active \n");
        return true;
    }
    s.Set(1.0, *stmp);

    
    ///////////////////////////////////////////////////////// residual update
    r->linComb(alpha, *Hp, 1.0);
  } // End CG loop.

  iflag = 3;
  if (iprint)
    printf(" maximum number of iterations reached \n");
  
  return true;

}  // End runTangentialStep.


bool SQPAlgo::runAcceptStep(Vector &x, const Vector &l, double f, const Vector &g, const Vector &c,
                            const Vector &v, const Vector &w, double &delta, double &rho, double tol,
                            int &iaccept)
{
  bool iprint = true;  // print flag  = true print output
                       // otherwise no output is generated

  double eta  = 1e-8;     // actual/predicted-reduction parameter
  double beta = 1e-8;     // predicted reduction parameter

  iaccept = 0;

  // Create necessary vectors.
  VectorPtr s     = x.createVector();
  VectorPtr Jl    = x.createVector();
  VectorPtr Jsc   = c.createVector();
  VectorPtr tryx  = x.createVector();
  VectorPtr g_new = x.createVector();
  VectorPtr c_new = c.createVector();
  VectorPtr l_new = c.createVector();
  VectorPtr gJl   = x.createVector();
  VectorPtr dl    = c.createVector();
  VectorPtr Hs    = x.createVector();


  // Compute composite step.
  s->Set(1.0, v);
  s->linComb(1.0, w, 1.0);

  // Compute and store some quantities for later use. Necessary
  // because dat_->updateAll updates the constraint derivatives.
  constr_->applyJacobian(true, x, l, *Jl);
  constr_->applyJacobian(false, x, *s, *Jsc);
  Jsc->linComb(1.0, c, 1.0);
  double Jsc_nrm2  = Jsc->innerProd(*Jsc);

  // Compute objective, constraint, etc. values at the trial point.
  tryx->Set(1.0, x);
  tryx->linComb(1.0, *s, 1.0);
  dat_->computeAll(*tryx);
  double f_new = obj_->getValue(*tryx);
  obj_->getGradient(*tryx, *g_new);
  constr_->getValue(*tryx, *c_new);
  runMultiplierEst(*tryx, *g_new, *l_new, tol);

  // Compute actual reduction.
  double ared = (f + l.innerProd(c) + rho*c.innerProd(c)) -
                (f_new + l_new->innerProd(*c_new) + rho*c_new->innerProd(*c_new));

  // Compute predicted reduction.
  gJl->Set(1.0, g);
  gJl->linComb(1.0, *Jl, 1.0);
  dl->Set(1.0, l);
  dl->linComb(1.0, *l_new, -1.0);
  hessvec_->getValue(x, l, *s, *Hs);
  double sHs = s->innerProd(*Hs);
  double pred = (f + l.innerProd(c) + rho*c.innerProd(c)) -
                ( f + l.innerProd(c) + gJl->innerProd(*s) + 0.5*sHs ) +
                dl->innerProd(*Jsc) + rho*Jsc_nrm2;

  // Penalty parameter update.
  if ( pred < 0.5*rho*(c.innerProd(c)-Jsc_nrm2) )
    rho = 2 * ( gJl->innerProd(*s) + sHs + dl->innerProd(*Jsc) ) / ( c.innerProd(c) - Jsc_nrm2 ) + beta;

  // Determine if the step gives sufficient reduction in the merit function,
  // update the trust-region radius.
  double r  = ared/pred;
  double ns = sqrt(s->innerProd(*s));
  if (r >= eta) {
    x.linComb(1.0, *s, 1.0);
    if (r >= 0.9)
        delta = (7*ns <= delta) ? delta : 7*ns;  // delta = max(7*ns, delta);
    else if (r >= 0.3)
        delta = (2*ns <= delta) ? delta : 2*ns;  // delta = max(2*ns, delta);
    iaccept = 1;
  }
  else
    delta = 0.5*ns;
    
  return true;

}


bool SQPAlgo::runDerivativeCheck(const Vector &x, const Vector &l, const Vector &dir) {

  // Vector variables.
  VectorPtr g   = x.createVector();      // gradient of objective
  VectorPtr gl  = x.createVector();      // gradient of Lagrangian
  VectorPtr gln = x.createVector();      // new gradient of Lagrangian
  VectorPtr gldiff = x.createVector();   // difference of old and new gradients of Lagrangian
  VectorPtr h = x.createVector();        // hessvec
  VectorPtr hdiff = x.createVector();    // h - (gln-gl)/delta
  VectorPtr xdd = x.createVector();      // xdd = x + delta*dir
  //VectorPtr w   = x.createVector();    // tangential step
  VectorPtr c   = l.createVector();      // constraint vector
  VectorPtr cn  = l.createVector();      // new constraint vector
  VectorPtr Jd  = l.createVector();      // constraint jacobian * dir
  VectorPtr cdiff  = l.createVector();   // cn-c
  VectorPtr Jdiff  = l.createVector();   // Jd - (cn-c)/delta
  
  // Check the computation of the gradient of the objective function.
  printf(" Gradient check using finite differences (FDs)\n");
  printf(" FD step size      grad'*v      FD approx.   absolute error \n");
  dat_->computeAll(x);
  double f = obj_->getValue(x);
  obj_->getGradient(x, *g);

  double delta = 1.0e-1;
  double f1    = 0.0;
  for (int d = 1; d <= 9; d++) {
    delta = 0.1*delta;
    xdd->Set(1.0, x);
    xdd->linComb(delta, dir, 1.0);
    dat_->computeAll(*xdd);
    f1 = obj_->getValue(*xdd);
    printf(" %12.6e   %12.6e    %12.6e  %12.6e  \n",
            delta, g->innerProd(dir), (f1-f)/delta, abs(g->innerProd(dir) - (f1-f)/delta) );
  }

  // Check the computation of the Jacobian of the constraints.
  printf("\n Jacobian check using finite differences (FDs)\n");
  printf(" FD step size      norm(Jac*v)     norm(FD approx.)   absolute error \n");
  dat_->computeAll(x);
  constr_->getValue(x, *c);
  constr_->applyJacobian(false, x, dir, *Jd); 
  delta       = 1.0e-1;
  for (int d = 1; d <= 9; d++) {
    delta = 0.1*delta;
    xdd->Set(1.0, x);
    xdd->linComb(delta, dir, 1.0);
    dat_->computeAll(*xdd);
    constr_->getValue(*xdd, *cn);
    cdiff->Set(-1.0, *c);
    cdiff->linComb(1.0, *cn, 1.0);
    Jdiff->Set(-1.0/delta, *cdiff);
    Jdiff->linComb(1.0, *Jd, 1.0);
    printf(" %12.6e     %12.6e      %12.6e      %12.6e  \n",
             delta, sqrt(Jd->innerProd(*Jd)), sqrt(cdiff->innerProd(*cdiff))/delta,
             sqrt(Jdiff->innerProd(*Jdiff)));
  }                                  


  // Check the computation of the Hessian of the Lagrangian.
  printf("\n Hessian od Lagrangian check using finite differences (FDs)\n");
  printf(" FD step size      norm(H*v)      norm(FD approx.)    absolute error \n");
  dat_->computeAll(x);
  obj_->getGradient(x, *g);
  constr_->applyJacobian(true, x, l, *gl); 
  gl->linComb(1.0, *g, -1.0);
  delta = 1.0e-1;
  hessvec_->getValue(x, l, dir, *h);
  for (int d = 1; d <= 9; d++) {
    delta = 0.1*delta;
    xdd->Set(1.0, x);
    xdd->linComb(delta, dir, 1.0);
    dat_->computeAll(*xdd);
    obj_->getGradient(*xdd, *g);
    constr_->applyJacobian(true, *xdd, l, *gln); 
    gln->linComb(1.0, *g, -1.0);
    gldiff->Set(-1.0, *gl);
    gldiff->linComb(1.0, *gln, 1.0);
    hdiff->Set(-1.0/delta, *gldiff);
    hdiff->linComb(1.0, *h, 1.0);
    printf(" %12.6e     %12.6e      %12.6e      %12.6e  \n",
             delta, sqrt(h->innerProd(*h)), sqrt(gldiff->innerProd(*gldiff))/delta,
             sqrt(hdiff->innerProd(*hdiff)));
  }

}


}  // End adding to namespace Aristos.
