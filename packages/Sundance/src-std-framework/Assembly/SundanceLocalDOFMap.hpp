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



#ifndef SUNDANCE_LOCALDOFMAP_H
#define SUNDANCE_LOCALDOFMAP_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"


namespace SundanceStdFwk
{
namespace Internal
{
/** 
 * LocalDOFMap bundles several compact tables used for fast lookup of
 * local DOFs. 
 */
class LocalDOFMap
{
public:
  /** */
  LocalDOFMap(int numBlocks, int verb);

  /** */
  void setCells(int cellDim, int maxCellDim,
    const RefCountPtr<const Array<int> >& cellLID);

  /** */
  int nCells() const ;

  /** */
  bool isUsed(int b) const {return isUsed_[b];}

  /** */
  bool isUnused(int b) const {return !isUsed(b);}

  /** */
  bool isUnused() const ;

  /** */
  void markAsUnused() ;

  /** */
  bool hasCells() const {return hasCells_;}

  /** */
  const RefCountPtr<const Array<int> >& cellLIDs() const {return cellLID_;}

  /** */
  void markAsUsed(int b) {isUsed_[b] = true;}

  /** */
  int numBlocks() const {return mapStruct_->size();}

  /** */
  const Array<int>& nLocalNodesPerChunk(int b) const 
    {return (*nLocalNodesPerChunk_)[b];}

  /** */
  const RefCountPtr<const MapStructure>& mapStruct(int b) const 
    {return (*mapStruct_)[b];}

  /** */
  const Array<Array<int> >& localDOFs(int b) const 
    {return (*localDOFs_)[b];}

  /** */
  std::ostream& print(std::ostream& os) const ;

  /** */
  void fillBlock(int b, const RefCountPtr<DOFMapBase>& globalMap,
    const Array<Set<int> >& requiredFunc);

private:


  /** */
  Array<int>& nLocalNodesPerChunk(int b)
    {return (*nLocalNodesPerChunk_)[b];}

  /** */
  RefCountPtr<const MapStructure>& mapStruct(int b) 
    {return (*mapStruct_)[b];}

  /** */
  Array<Array<int> >& localDOFs(int b) 
    {return (*localDOFs_)[b];}

  /** */
  void verifyValidBlock(int b) const ;


  int verb_;
  Array<int> isUsed_;
  bool hasCells_;
  RefCountPtr<Array<Array<int> > > nLocalNodesPerChunk_;
  RefCountPtr<Array<RefCountPtr<const MapStructure> > > mapStruct_;  
  RefCountPtr<Array<Array<Array<int> > > > localDOFs_;
  RefCountPtr<const Array<int> > cellLID_;
  int activeCellDim_;
  int maxCellDim_;
};

}
}

namespace SundanceStdFwk{
namespace Internal
{
inline std::ostream& operator<<(std::ostream& os,
  const LocalDOFMap& m)
{
  return m.print(os);
}
}
}



#endif
