/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"
#include "Teuchos_HashSet.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    typedef std::set<int> ColSetType;

    /** 
     * 
     */
    class Assembler : public TSFExtended::ObjectWithVerbosity<Assembler>
    {
    public:
      /** */
      Assembler(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn,
                const VectorType<double>& vectorType,
                const VerbositySetting& verb = classVerbosity());
      /** */
      Assembler(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn,
                const VerbositySetting& verb = classVerbosity());
      
      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const 
      {return rowMap_;}

      /** */
      const RefCountPtr<DOFMapBase>& colMap() const 
      {return colMap_;}

      /** */
      const RefCountPtr<DiscreteSpace>& solutionSpace() const 
      {return colSpace_;}

      /** */
      const RefCountPtr<DiscreteSpace>& rowSpace() const 
      {return rowSpace_;}

      /** */
      const RefCountPtr<Set<int> >& bcRows() {return bcRows_;}

      /** */
      void assemble(TSFExtended::LinearOperator<double>& A,
                    TSFExtended::Vector<double>& b) const ;


      /** */
      void assemble(TSFExtended::Vector<double>& b) const ;

      /** */
      void evaluate(double& value,
                    TSFExtended::Vector<double>& gradient) const ;

      /** */
      void evaluate(double& value) const ;

      /** */
      static unsigned int& workSetSize() 
      {static unsigned int rtn = defaultWorkSetSize(); return rtn;}
      
      /** */
      void getGraph(Array<int>& graphData,
                    Array<int>& rowPtrs,
                    Array<int>& nnzPerRow) const ;

      /** */
      void flushConfiguration() 
      {
        vecNeedsConfiguration_ = true;
        matNeedsConfiguration_ = true;
      }

      /** */
      static int& numAssembleCalls() {static int rtn=0; return rtn;}
      
    private:

      /** */
      void init(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn);

      /** */
      void insertLocalMatrixBatch(int cellDim, const Array<int>& workSet,
                                  bool isBCRqc, 
                                  const Array<int>& testIndices,
                                  const Array<int>& unkIndices,
                                  unsigned int nTestNodes, 
                                  unsigned int nUnkNodes,
                                  const Array<int>& testID,
                                  const Array<int>& unkID, 
                                  const Array<double>& localValues,
                                  TSFExtended::LoadableMatrix<double>* mat) const ;

      /** */
      void insertLocalVectorBatch(int cellDim, const Array<int>& workSet,
                                  bool isBCRqc, 
                                  const Array<int>& testIndices,
                                  unsigned int nTestNodes, 
                                  const Array<int>& testID,
                                  const Array<double>& localValues,
                                  TSFExtended::LoadableVector<double>* vec) const ;

      /** */
      void configureMatrix(LinearOperator<double>& A,
                           Vector<double>& b) const ;

      /** */
      void configureVector(Vector<double>& b) const ;

      /** */
      bool isBCRow(int dof) const {return (*isBCRow_)[dof-lowestRow_];}

      /** */
      static int defaultWorkSetSize() {return 100;}
      
      mutable bool matNeedsConfiguration_;

      mutable bool vecNeedsConfiguration_;

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      RefCountPtr<DOFMapBase> rowMap_;

      RefCountPtr<DOFMapBase> colMap_;

      RefCountPtr<DiscreteSpace> rowSpace_;

      RefCountPtr<DiscreteSpace> colSpace_;

      RefCountPtr<Set<int> > bcRows_;

      Array<RegionQuadCombo> rqc_;

      Map<ComputationType, Array<EvalContext> > contexts_;

      Array<int> isBCRqc_;

      Map<ComputationType, Array<Array<IntegralGroup> > > groups_;

      Array<RefCountPtr<StdFwkEvalMediator> > mediators_;

      Map<ComputationType, Array<const EvaluatableExpr*> > evalExprs_;

      RefCountPtr<EvalManager> evalMgr_;

      RefCountPtr<Array<int> > isBCRow_;

      unsigned int lowestRow_;

      VectorType<double> vecType_;

    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
