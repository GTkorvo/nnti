/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFANASAZIEIGENSOLVER_HPP
#define TSFANASAZIEIGENSOLVER_HPP

#include "SundanceDefs.hpp"
#include "TSFVectorDecl.hpp" 
#include "TSFSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFEigensolverBase.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziThyraAdapter.hpp"


namespace TSFExtended
{
using Teuchos::ParameterList;

/**
 * Object wrapper for Anasazi eigenvalue solver.
 */
template <class Scalar>
class AnasaziEigensolver
  : public EigensolverBase<Scalar>,
    public Sundance::Handleable<EigensolverBase<Scalar> >
{
public:

  /** */
  typedef Thyra::MultiVectorBase<Scalar>         MV;
  /** */
  typedef Thyra::LinearOpBase<Scalar>            OP;
  /** */
  typedef Anasazi::MultiVecTraits<Scalar,MV>     MVT;
  /** */
  typedef Anasazi::OperatorTraits<Scalar,MV,OP>  OPT;

  /** */
  AnasaziEigensolver(const ParameterList& params) 
    : EigensolverBase<Scalar>(params) {;}

  /**
   * Solve a generalized eigenvalue problem \f$ K x = \lambda M x \f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    const LinearOperator<Scalar>& M,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const ;

  /** \name Handleable interface */
  //@{
  /** Return a ref counted pointer to a newly created object */
  virtual RCP<EigensolverBase<Scalar> > getRcp() 
    {return rcp(this);}
    //@}
  
private:

  static Time& solveTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("AnasaziEigensolver solve()"); 
      return *rtn;
    }

  static Time& precondBuildTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("AnasaziEigensolver building preconditioner"); 
      return *rtn;
    }

};


template <class Scalar>  
inline void AnasaziEigensolver<Scalar>::solve(
  const LinearOperator<Scalar>& K,
  const LinearOperator<Scalar>& M,
  Array<Vector<Scalar> >& evecs,
  Array<std::complex<Scalar> >& ew) const 
{
  TimeMonitor timer(solveTimer());
  VectorSpace<Scalar> KDomain = K.domain();
  
  /* Get a Thyra representation of the stiffness matrix */
  RCP<LinearOpBase<Scalar> > KPtr = K.ptr();
  RCP<LinearOpBase<Scalar> > MPtr = M.ptr();
  RCP<const Thyra::VectorSpaceBase<Scalar> > mvSpace = KPtr->domain();
  
  
  // Eigensolver parameters
  string method = this->params().get<string>("Method");
  int numEigs = this->params().get<int>("Number of Eigenvalues");
  int blockSize = this->params().get<int>("Block Size");
  bool usePrec = this->params().get<bool>("Use Preconditioner");
  bool hermitian = this->params().get<bool>("Is Hermitian");

  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
  if (this->params().get<int>("Verbose"))
  {
    verbosity += Anasazi::Debug + Anasazi::IterationDetails;
  }
  
  /* Make a multivector with row space = domain of K, column 
   * space = multiVec Space*/
  RCP<MV> mv = Thyra::createMembers( *mvSpace, blockSize );
  
  /* Fill the multivector with random values */
  MVT::MvRandom( *mv );

  /* Create a preconditioner */
  ParameterList precParams = this->params().sublist("Preconditioner");
  PreconditionerFactory<double> precFactory 
    = new ParameterListPreconditionerFactory(precParams);

  LinearOperator<Scalar> P;
  if (usePrec) 
  {
    TimeMonitor pTimer(precondBuildTimer());
    P = precFactory.createPreconditioner(K).right();
  }

  /* Create eigenproblem */
  RCP<Anasazi::Eigenproblem<Scalar,MV,OP> > problem;

  if (MPtr.get() != 0)
  {
    problem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(KPtr, MPtr, mv) );
  }
  else
  {
    problem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(KPtr, mv) );
  }

  ParameterList eigParams = this->params();
  problem->setHermitian(hermitian);
  problem->setNEV(numEigs);
  if (usePrec) problem->setPrec(P.ptr());
  bool ret = problem->setProblem();
  TEST_FOR_EXCEPTION(!ret, std::runtime_error,
    "Eigenproblem not setup correctly");
  

  // Create the solver manager
  RCP<Anasazi::SolverManager<Scalar,MV,OP> > MySolverMan;
  if (method=="Block Davidson")
  {
    MySolverMan = rcp(new Anasazi::BlockDavidsonSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="Block Krylov Schur")
  {
    MySolverMan = rcp(new Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="Simple LOBPCG")
  {
    MySolverMan = rcp(new Anasazi::SimpleLOBPCGSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="LOBPCG")
  {
    MySolverMan = rcp(new Anasazi::LOBPCGSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else
  {
    TEST_FOR_EXCEPTION(true, std::runtime_error,
      "solver method [" << method << "] not recognized");
  }

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan->solve();
  TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, 
    std::runtime_error, "Anasazi did not converge!");
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<Scalar,MV> sol = problem->getSolution();
  RCP<MV> evecs_mv = sol.Evecs;
  int numev = sol.numVecs;
  
  /* Copy the columns of the eigenvector MV into an array of TSF vectors */
  ew.resize(numev);
  evecs.resize(numev);

  for (int i=0; i<numev; i++)
  {
    Vector<Scalar> tmp = evecs_mv->col(i);
    evecs[i] = KDomain.createMember();
    evecs[i].acceptCopyOf(tmp);
    /* record the associated eigenvalue. The matrix is Hermitian so
     * we know the eigenvalue is real. */
    //evals[i] = sol.Evals[i].realpart;
    // if matrix might not be hermitian
    ew[i].real() = sol.Evals[i].realpart;
    ew[i].imag() = sol.Evals[i].imagpart;
  }
}


}


#endif
