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

#ifndef SUNDANCE_FUNCELEMENTBASE_H
#define SUNDANCE_FUNCELEMENTBASE_H


#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
using namespace SundanceUtils;
using namespace Teuchos;

using std::string;
using std::ostream;

namespace Internal
{
/** 
 * FuncElementBase defines the interface for scalar-valued elements
 * of Sundance functions. At the user level, Sundance functions can be
 * list (e.g, vector or tensor) valued; internally, however, compound
 * expressions use only scalar functions deriving from the 
 * FuncElementBase class. 
 *
 * <h4> Function ID </h4>
 *
 * Every function element has a unique integer 
 * identifier, or <b> funcID, </b>
 * which is used internally as a shorthand for the function. The funcID
 * is assigned at construction using the private static method
 * nextID(), and because it is private with no mutators, can never
 * be changed. The nextID() method is responsible for guaranteeing that 
 * funcIDs are unique.
 */
class FuncElementBase : virtual public ScalarExpr
{
public:
  /** */
  FuncElementBase(const string& rootName,
    const string& suffix,
    int sharedID);
  /** */
  FuncElementBase(const string& rootName,
    int sharedID);

  /** virtual destructor */
  virtual ~FuncElementBase() {;}

  /** Return an integer ID which uniquely identifies this
   * vector component of this function */
  int funcComponentID() const {return componentID_;}

  /** Return an integer ID shared by all components of a vector-valued
   * function. */
  int sharedFuncID() const {return sharedID_;}

  /** Append to the set of func IDs present in this expression. */
  virtual void accumulateFuncSet(Set<int>& funcIDs, 
    const Set<int>& activeSet) const ;

  /** Return the name of this function */
  const string& name() const {return name_;}

  /** Return the root name of this function */
  const string& rootName() const {return rootName_;}

  /** Return the root name of this function */
  const string& suffix() const {return suffix_;}

  /** Write self in text form */
  virtual ostream& toText(ostream& os, bool paren) const ;

  /** Write self in Latex form */
  virtual ostream& toLatex(ostream& os, bool paren) const ;
      
  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


protected:
  /** Determine whether this function is in the given active set */
  bool isInActiveSet(const Set<MultiSet<int> >& activeFuncIDs) const ;
          
private:

  string name_;

  string rootName_;

  string suffix_;

  int componentID_;

  int sharedID_;

  static int nextComponentID() {static int rtn = 0; rtn++; return rtn;}

};
}
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
