/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_NONLINEARPROBLEM_H
#define SUNDANCE_NONLINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFNonlinearOperatorBase.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

    /** 
     * 
     */
  class NonlinearProblem 
    : public TSFExtended::ObjectWithVerbosity<NonlinearProblem>,
      public TSFExtended::NonlinearOperatorBase<double>
    {
    public:
      /** Empty ctor */
      NonlinearProblem();

      /** Construct with a mesh, equation set, bcs, test and unknown funcs,
       * and a vector type */
      NonlinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
                       const Expr& test, const Expr& unk, const Expr& u0, 
                       const TSFExtended::VectorType<double>& vecType);

      /** Compute the residual and Jacobian at the current evaluation point */
      LinearOperator<double> computeJacobianAndFunction(Vector<double>& functionValue) const ;

      /** Compute the residual at the current eval point */
      TSFExtended::Vector<double> computeFunctionValue() const ;

      /** Get an initial guess */
      TSFExtended::Vector<double> getInitialGuess() const ;

      /* Handle boilerplate */
      GET_RCP(TSFExtended::NonlinearOperatorBase<double>);

    private:
      
      /** */
      RefCountPtr<Assembler> assembler_;

      /** */
      mutable TSFExtended::LinearOperator<double> J_;

      /** */
      Expr u0_;

      /** */
      mutable DiscreteFunction* discreteU0_;
      
    };
}


#endif
