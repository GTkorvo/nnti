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

#ifndef SUNDANCE_DOFMAPBUILDER_H
#define SUNDANCE_DOFMAPBUILDER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceBasisFamily.hpp"
#include "TSFObjectWithVerbosity.hpp"

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

    /** 
     * 
     */
    class DOFMapBuilder : public TSFExtended::ObjectWithVerbosity<DOFMapBase>
    {
    public:
      /** */
      DOFMapBuilder();
      /** */
      DOFMapBuilder(const Mesh& mesh, const RefCountPtr<EquationSet>& eqn);

      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const {return rowMap_;}

      /** */
      const RefCountPtr<DOFMapBase>& colMap() const {return colMap_;}

      /** */
      const RefCountPtr<Array<int> >& isBCRow() const {return isBCRow_;}

      Array<BasisFamily> testBasisArray() const ;

      Array<BasisFamily> unkBasisArray() const ;

      const Mesh& mesh() const {return mesh_;}

    private:

      bool hasUnks() const ;

      bool unksAreOmnipresent() const ;

      bool testsAreOmnipresent() const ;

      bool regionIsMaximal(int r) const ;

      bool isSymmetric() const ;

      void markBCRows() ;

      const MPIComm& comm() const {return mesh().comm();}

      void init();

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      RefCountPtr<DOFMapBase> rowMap_;

      RefCountPtr<DOFMapBase> colMap_;

      RefCountPtr<Array<int> > isBCRow_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
