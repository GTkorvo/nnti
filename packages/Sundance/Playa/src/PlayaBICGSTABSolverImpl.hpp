/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BICGSTABSOLVER_IMPL_HPP
#define PLAYA_BICGSTABSOLVER_IMPL_HPP

#include "PlayaBICGSTABSolverDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverBaseImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaSimpleComposedOpImpl.hpp"
#endif



namespace Playa
{
using namespace Teuchos;


/* */
template <class Scalar> inline
BICGSTABSolver<Scalar>
::BICGSTABSolver(const ParameterList& params)
  : KrylovSolver<Scalar>(params) {;}

/* */
template <class Scalar> inline
BICGSTABSolver<Scalar>::BICGSTABSolver(const ParameterList& params,
  const PreconditionerFactory<Scalar>& precond)
  : KrylovSolver<Scalar>(params, precond) {;}

/* Write to a stream  */
template <class Scalar> inline
void BICGSTABSolver<Scalar>::print(std::ostream& os) const 
{
  os << description() << "[" << std::endl;
  os << this->parameters() << std::endl;
  os << "]" << std::endl;
}

    
template <class Scalar> inline
SolverState<Scalar> BICGSTABSolver<Scalar>
::solveUnprec(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& b,
  Vector<Scalar>& soln) const
{
  int maxiters = this->getMaxiters();
  Scalar tol = this->getTol();
  int verbosity = this->verb();

  Scalar normOfB = sqrt(b.dot(b));

  /* check for trivial case of zero rhs */
  if (normOfB < tol) 
  {
    soln = b.space().createMember();
    soln.zero();
    return SolverState<Scalar>(SolveConverged, "RHS was zero", 0, 0.0);
  }

  /* check for initial zero residual */
  Vector<Scalar> x0 = b.copy();
  Vector<Scalar> r0 = b.space().createMember();
  Vector<Scalar> tmp = b.space().createMember();

  // r0 =  b - op*x0;
  op.apply(x0, tmp);
  r0 = b - tmp;
	
  if (sqrt(r0.dot(r0)) < tol*normOfB) 
  {
    soln = x0;
    return SolverState<Scalar>(SolveConverged, "initial resid was zero", 
      0, 0.0);
  }

  Vector<Scalar> p0 = r0.copy();
  //    p0.randomize();
  Vector<Scalar> r0Hat = r0.copy();
  Vector<Scalar> xMid = b.space().createMember();
  Vector<Scalar> rMid = b.space().createMember();
  Vector<Scalar> ArMid = b.space().createMember();
  Vector<Scalar> x = b.space().createMember();
  Vector<Scalar> r = b.space().createMember();
  Vector<Scalar> s = b.space().createMember();
  Vector<Scalar> ap = b.space().createMember();

  int myRank = MPIComm::world().getRank();

  Scalar resid = -1.0;

  for (int k=1; k<=maxiters; k++)
  {
    // ap = A*p0
    op.apply(p0, ap);

    Scalar den = ap.dot(r0Hat);
    if (Utils::chop(sqrt(fabs(den))/normOfB)==0) 
    {
      SolverState<Scalar> rtn(SolveCrashed, 
        "BICGSTAB failure mode 1", k, resid);
      return rtn;
    }
			
    Scalar a0 = r0.dot(r0Hat)/den;
			
    xMid = x0 + a0*p0;
    //xMid.axpy(a0, p0, x0);

    rMid = r0 - a0*ap;
    //rMid.axpy(-a0, ap, r0);

    // check for convergence
    Scalar resid = rMid.norm2()/normOfB;
    if (resid < tol) 
    {
      soln = xMid; 
      SolverState<Scalar> rtn(SolveConverged, "yippee!!", k, resid);
      return rtn;
    }

    // ArMid = A*rMid
    op.apply(rMid, ArMid);

    den = ArMid.dot(ArMid);
    if (Utils::chop(sqrt(fabs(den))/normOfB)==0)  
    {
      SolverState<Scalar> rtn(SolveCrashed, 
        "BICGSTAB failure mode 2", k, resid);
      return rtn;
    }

    Scalar w = rMid.dot(ArMid)/den;
			
    x = xMid + w*rMid;
    //x.axpy(w, rMid, xMid);
			
    r = rMid - w*ArMid;
    //r.axpy(-w, ArMid, rMid);

    // check for convergence
    resid = sqrt(r.dot(r))/normOfB;
    if (resid < tol) 
    {
      soln = x;
      SolverState<Scalar> rtn(SolveConverged, "yippee!!", k, resid);
      return rtn;
    }

    den = w*(r0.dot(r0Hat));
    if (Utils::chop(sqrt(fabs(den))/normOfB)==0) 
    {
      SolverState<Scalar> rtn(SolveCrashed, 
        "BICGSTAB failure mode 3", k, resid);
      return rtn;
    }
    Scalar beta = a0*(r.dot(r0Hat))/den;

    s = p0 - w*ap;
    p0 = r + beta*s;
//    p0 = r + beta*p0 - beta*w*ap;
    //s.axpy(-w, ap, p0);
    //p0.axpy(beta, s, r);

    r0 = r.copy();
    x0 = x.copy();

    if (myRank==0 && verbosity > 1 ) 
    {
      Out::os() << "BICGSTAB: iteration=";
      Out::os().width(8);
      Out::os() << k;
      Out::os().width(20);
      Out::os() << " resid=" << resid << std::endl;
    }
  }
    
  SolverState<Scalar> rtn(SolveFailedToConverge, 
    "BICGSTAB failed to converge", 
    maxiters, resid);
  return rtn;
}


}

#endif
