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

#ifndef SUNDANCE_CELLPREDICATE_H
#define SUNDANCE_CELLPREDICATE_H



#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellPredicateBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;
  using namespace Internal;
  
    /** 
     * User-level handle for predicates (deriving from CellPredicateBase)
     * used to decide whether
     * a given cell passes a CellFilter.
     */
  class CellPredicate : public TSFExtended::Handle<CellPredicateBase>
  {
    public:
    
    /* handle boilerplate */
    HANDLE_CTORS(CellPredicate, CellPredicateBase);

    /** write to XML */
    XMLObject toXML() const {return ptr()->toXML();}

#ifndef DOXYGEN_DEVELOPER_ONLY

    /** set the mesh on which cells are to be tested */
    void setMesh(const Mesh& mesh, int cellDim) const 
    {ptr()->setMesh(mesh, cellDim);}

    /** compare to another predicate, used for placement in STL containers */
    bool operator<(const CellPredicate& other) const ;


#endif  /* DOXYGEN_DEVELOPER_ONLY */    
    };

}

namespace std
{
  inline ostream& operator<<(ostream& os, const SundanceStdFwk::CellPredicate& pred)
  {
    os << pred.toXML() << endl;
    return os;
  }
}


#endif
