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

#ifndef SUNDANCE_CELLFILTER_H
#define SUNDANCE_CELLFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellPredicate.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceOrderedHandle.hpp"
#include "SundanceHandle.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
using namespace Teuchos;
  

class CellFilterBase;

  
/** 
 * CellFilter is a user-level object representing a 
 * filter that selects from a mesh all cells that
 * satisfy some condition. CellFilters are used to identify subdomains
 * on which equations or boundary conditions apply. 
 * Examples of cell filters are to identify
 * all cells on the boundary of the mesh, or all cells whose node positions
 * satisfy some mathematical equation or inequality. 
 *
 * <h4> Use of CellFilter </h4>
 *
 * \code
 * // Define a filter that will pick out all maximal cells located 
 * // entirely in the left half-plane x < 0 
 * CellFilter elements = new MaximalCellFilter();
 * CellFilter leftHalf = elements.subset( x <= 0.0 );
 * 
 * // Apply the leftHalf filter to a mesh, thus enumerating
 * // all the cells of that mesh living in the left half-plane
 * CellSet leftCells = leftHalf.getCells(mesh);
 * \endcode
 *
 * <h4> Set operations with CellFilter objects </h4>
 *
 * Operations on cell filters can produce new filters. 
 *
 * The subset() and labeledSubset() operators
 * produce new CellFilters that pick out a subset of the cells
 * satisfying an additional condition given in the argument
 * to the subset methods. 
 *
 * Binary set operations can also produce new filters.
 * Suppose
 * <tt>a</tt> and <tt>b</tt> are CellFilters whose <tt>getCells()</tt>
 * methods produce
 * CellSets \f$\{A\}\f$ and \f$\{B\}\f$, respectively. There exist
 * operators for the following binary operations:
 * <ul>
 * <li> The <b>union</b> operator <tt>a+b.</tt> The result of a union
 * operation is a filter that will produce the union of the two operand's
 * cell sets, 
 * \f[{\tt a+b} \rightarrow \{A\} \cup \{B\}, \f]
 * i.e., all cells that are in either \f$\{A\}\f$ or \f$\{B\}\f$
 * <li> The <b>intersection </b> operator <tt>a.intersection(b)</tt> 
 * The result of an intersection
 * operation is a filter that will produce the intersection
 * of the two operand's
 * cell sets, 
 * \f[{\tt a.intersection(b)} \rightarrow \{A\} \cap \{B\}, \f]
 * i.e., all cells that are in both \f$\{A\}\f$ and \f$\{B\}\f$
 * <li> The <b>exclusion </b> operator <tt>a-b.</tt>  
 * The result of an exclusion
 * operation is a filter that will produce the exclusion
 * of the two operand's
 * cell sets, 
 * \f[{\tt a - b} \rightarrow \{A\} \setminus \{B\}, \f]
 * i.e., all cells that are in \f$\{A\}\f$ but not in \f$\{B\}\f$
 * </ul>
 * \code
 * CellFilter elements = new MaximalCellFilter();
 * CellFilter leftHalf = elements.subset( x <= 0.0 );
 * CellFilter topHalf = elements.subset( x >= 0.0 );
 * CellFilter topLeftQuarter = leftHalf.intersection(topHalf);
 * CellFilter 
 * \endcode
 *

*/
class CellFilter : public OrderedHandle<CellFilterStub>
{
public:
  ORDERED_HANDLE_CTORS(CellFilter, CellFilterStub);

  /** Find the cells passing this filter on the given mesh */
  CellSet getCells(const Mesh& mesh) const ;

  /** Return the dimension of this cell set on the given mesh */
  int dimension(const Mesh& mesh) const ;

  /** Return a filter that returns the set union of the sets
   * produced by the two operand filters */
  CellFilter operator+(const CellFilter& other) const ;
    
  /** Return a filter that returns the set difference between the sets
   * produced by the two operand filters */
  CellFilter operator-(const CellFilter& other) const ;

  /** Return a filter that returns the set intersection of the sets
   * produced by the two operand filters */
  CellFilter intersection(const CellFilter& other) const ;

  // /** Return a filter that returns the
  //      *  subset of cells for which the logical expression 
  //      * is true */
  //     CellFilter subset(const LogicalExpr& expr) const ;

    
  /** Return a filter that will return the subset of cells having
   * the given label */
  CellFilter labeledSubset(int label) const ;

    
  /** Return a filter that will return the subset of cells for which
   * the given predicate is true */
  CellFilter subset(const CellPredicate& test) const ;
    
  /** Return a filter that will return the subset of cells for which
   * the given predicate is true */
  CellFilter subset(const RCP<CellPredicateFunctorBase>& test) const ;

  /** Return a compilation of all registered subsets of this filter */
  const Set<CellFilter>& knownSubsets() const ;

  /** Return a compilation of all filters registered as disjoint with
   * this filter */
  const Set<CellFilter>& knownDisjoints() const ;

  /** Indicate whether this filter is known to be a subset of
   * the specified filter. Note that a negative answer here does NOT
   * mean that it is not a subset, only that it can't be determined
   * to be one through structural properties alone. If this function
   * returns false, then to get a definitive answer one must do a test
   * using an actual realization on a mesh. */
  bool isKnownSubsetOf(const CellFilter& other) const ;

  /** Indicate whether this filter is known to be disjoint with
   * the specified filter. Note that a negative answer here does NOT
   * mean that it is not disjoint, only that it can't be determined
   * to be one through structural properties alone. If this function
   * returns false, then to get a definitive answer one must do a test
   * using an actual realization on a mesh. */
  bool isKnownDisjointWith(const CellFilter& other) const ;

  /** Do a brute-force check of whether this filter is a subset of
   * the specified filter. */
  bool isSubsetOf(const CellFilter& other, const Mesh& mesh) const ;

  /** Do a brute-force check of whether this filter is disjoint with
   * the specified filter. */
  bool isDisjointWith(const CellFilter& other, const Mesh& mesh) const ;

    
  /** Indicate whether this is a null cell filter */
  bool isNullCellFilter() const ;


  /** Determine whether all cells identified by this filter are
   * facets of cells identified by the other filter */
  bool areFacetsOf(const CellFilter& other, const Mesh& mesh) const ;

  /** Indicate whether this cell set is null */
  bool isNull() const ;

  /** */
  bool operator==(const CellFilter& other) const ;

  /** */
  bool operator!=(const CellFilter& other) const ;
    

  /** */
  XMLObject toXML() const ;

  /** */
  string toString() const ;

  /** */
  void setName(const string& name) ;

#ifndef DOXYGEN_DEVELOPER_ONLY

  /** */
  const CellFilterBase* cfbPtr() const ;
  /** */
  CellFilterBase* nonConstCfbPtr();
    
private:
#endif  /* DOXYGEN_DEVELOPER_ONLY */

};

/** \relates CellFilter 
    \brief Create an array with one entry 
*/
inline
Array<CellFilter> List(const CellFilter& a)
{
  Array<CellFilter> rtn(1, a);
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with two entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b)
{
  Array<CellFilter> rtn(2);
  rtn[0] = a;
  rtn[1] = b;
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with three entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c)
{
  Array<CellFilter> rtn(3);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with four entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d)
{
  Array<CellFilter> rtn(4);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with five entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d, const CellFilter& e)
{
  Array<CellFilter> rtn(5);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  return rtn;
}


/** \relates CellFilter 
    \brief Create an array with six entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d, const CellFilter& e,
  const CellFilter& f)
{
  Array<CellFilter> rtn(6);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with seven entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d, const CellFilter& e,
  const CellFilter& f, const CellFilter& g)
{
  Array<CellFilter> rtn(7);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with eight entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d, const CellFilter& e,
  const CellFilter& f, const CellFilter& g, const CellFilter& h)
{
  Array<CellFilter> rtn(8);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  return rtn;
}

/** \relates CellFilter 
    \brief Create an array with nine entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d, const CellFilter& e,
  const CellFilter& f, const CellFilter& g, const CellFilter& h, const CellFilter& i)
{
  Array<CellFilter> rtn(9);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  return rtn;
}


/** \relates CellFilter 
    \brief Create an array with ten entries 
*/
inline
Array<CellFilter> List(const CellFilter& a, const CellFilter& b, const CellFilter& c, const CellFilter& d, const CellFilter& e,
  const CellFilter& f, const CellFilter& g, const CellFilter& h, const CellFilter& i, const CellFilter& j)
{
  Array<CellFilter> rtn(10);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  return rtn;
}

/** \relates CellFilter */
typedef Array<CellFilter> CellFilterArray;

}

namespace std
{
inline std::ostream& operator<<(std::ostream& os, const Sundance::CellFilter& f)
{
  os << f.toString();
  return os;
}
}

#endif
