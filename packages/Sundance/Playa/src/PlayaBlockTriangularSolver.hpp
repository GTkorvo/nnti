/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BLOCKTRIANGULARSOLVER_HPP
#define PLAYA_BLOCKTRIANGULARSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverDecl.hpp" 
#include "PlayaLinearCombinationDecl.hpp" 
#include "PlayaCommonOperatorsDecl.hpp" 


namespace Playa
{
  /** */
  template <class Scalar>
  class BlockTriangularSolver : public LinearSolverBase<Scalar>,
                                public Playa::Handleable<LinearSolverBase<Scalar> >
  {
  public:
    /** */
    BlockTriangularSolver(const LinearSolver<Scalar>& solver)
      : LinearSolverBase<Scalar>(ParameterList()), solvers_(tuple(solver)) {;}

    /** */
    BlockTriangularSolver(const Array<LinearSolver<Scalar> >& solvers)
      : LinearSolverBase<Scalar>(ParameterList()), solvers_(solvers) {;}

    /** */
    virtual ~BlockTriangularSolver(){;}

    /** */
    virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                                      const Vector<Scalar>& rhs,
                                      Vector<Scalar>& soln) const ;

    /* */
    GET_RCP(LinearSolverBase<Scalar>);
  private:
    Array<LinearSolver<Scalar> > solvers_;
  };


  template <class Scalar> inline
  SolverState<Scalar> BlockTriangularSolver<Scalar>
  ::solve(const LinearOperator<Scalar>& op,
          const Vector<Scalar>& rhs,
          Vector<Scalar>& soln) const
  {
    int nRows = op.numBlockRows();
    int nCols = op.numBlockCols();

    soln = op.domain().createMember();
    //    bool converged = false;

    TEST_FOR_EXCEPTION(nRows != rhs.space().numBlocks(), std::runtime_error,
                       "number of rows in operator " << op
                       << " not equal to number of blocks on RHS "
                       << rhs);

    TEST_FOR_EXCEPTION(nRows != nCols, std::runtime_error,
                       "nonsquare block structure in block triangular "
                       "solver: nRows=" << nRows << " nCols=" << nCols);

    bool isUpper = false;
    bool isLower = false;

    for (int r=0; r<nRows; r++)
      {
        for (int c=0; c<nCols; c++)
          {
            if (op.getBlock(r,c).ptr().get() == 0 ||
                dynamic_cast<const SimpleZeroOp<Scalar>* >(op.getBlock(r,c).ptr().get()))
              {
                TEST_FOR_EXCEPTION(r==c, std::runtime_error,
                                   "zero diagonal block (" << r << ", " << c 
                                   << " detected in block "
                                   "triangular solver. Operator is " << op);
                continue;
              }
            else
              {
                if (r < c) isUpper = true;
                if (c < r) isLower = true;
              }
          }
      }

    TEST_FOR_EXCEPTION(isUpper && isLower, std::runtime_error, 
                       "block triangular solver detected non-triangular operator "
                       << op);

    bool oneSolverFitsAll = false;
    if ((int) solvers_.size() == 1 && nRows != 1) 
      {
        oneSolverFitsAll = true;
      }

    for (int i=0; i<nRows; i++)
      {
        int r = i;
        if (isUpper) r = nRows - 1 - i;
        Vector<Scalar> rhs_r = rhs.getBlock(r);
        for (int j=0; j<i; j++)
          {
            int c = j;
            if (isUpper) c = nCols - 1 - j;
            if (op.getBlock(r,c).ptr().get() != 0)
              {
                rhs_r = rhs_r - op.getBlock(r,c) * soln.getBlock(c);
              }
          }

        SolverState<Scalar> state;
        Vector<Scalar> soln_r;
        if (oneSolverFitsAll)
          {
            state = solvers_[0].solve(op.getBlock(r,r), rhs_r, soln_r);
          }
        else
          {
            state = solvers_[r].solve(op.getBlock(r,r), rhs_r, soln_r);
          }
        if (nRows > 1) soln.setBlock(r, soln_r);
        else soln = soln_r;
        if (state.finalState() != SolveConverged)
          {
            return state;
          }
      }

    return SolverState<Scalar>(SolveConverged, "block solves converged",
                               0, ScalarTraits<Scalar>::zero());
  }
  
}

#endif
