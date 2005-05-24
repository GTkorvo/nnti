// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
// Lesser General Public License for more details//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_STD_OPS_DECL_HPP
#define THYRA_VECTOR_STD_OPS_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {

/** \defgroup Thyra_Op_Vec_VectorStdOps_grp Collection of all vector operations
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */

/** \defgroup Thyra_Op_Vec_VectorStdOpsAll_grp Collection of vector operations for all scalar types.
 *
 * \ingroup Thyra_Op_Vec_VectorStdOps_grp
 */
//@{

/** \brief Sum of vector elements: <tt>result = sum( v(i), i = 1...v.space()->dim() )</tt>.
 */
template<class Scalar>
Scalar sum( const VectorBase<Scalar>& v );

/** \brief Natural norm: <tt>result = sqrt(<v,v>)</tt>.
 *
 * Returns <tt>Teuchos::ScalarTraits<Scalar>::squareroot(v.space()->scalarProd(v,v))</tt>.
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm( const VectorBase<Scalar>& v );

/** \brief One (1) norm: <tt>result = ||v||1</tt>.
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_1( const VectorBase<Scalar>& v );

/** \brief Euclidean (2) norm: <tt>result = ||v||2</tt>.
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_2( const VectorBase<Scalar>& v );

/** \brief Weighted Euclidean (2) norm: <tt>result = sqrt( sum( w(i)*conj(v(i))*v(i)) )</tt>.
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_2( const VectorBase<Scalar> &w, const VectorBase<Scalar>& v );

/** \brief Infinity norm: <tt>result = ||v||inf</tt>.
 */
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
norm_inf( const VectorBase<Scalar>& v_rhs );

/** \brief Dot product: <tt>result = conj(x)'*y</tt>.
 */
template<class Scalar>
Scalar dot( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y );

/** \brief Get single element: <tt>result = v(i)</tt>.
 */
template<class Scalar>
Scalar get_ele( const VectorBase<Scalar>& v, Index i );

/** \brief Set single element: <tt>v(i) = alpha</tt>.
 */
template<class Scalar>
void set_ele( Index i, Scalar alpha, VectorBase<Scalar>* v );

/** \brief Assign all elements to a scalar: <tt>y(i) = alpha, i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void assign( VectorBase<Scalar>* y, const Scalar& alpha );

/** \brief VectorBase assignment: <tt>y(i) = x(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void assign( VectorBase<Scalar>* y, const VectorBase<Scalar>& x );

/** \brief Add a scalar to all elements: <tt>y(i) += alpha, i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void Vp_S( VectorBase<Scalar>* y, const Scalar& alpha );

/** \brief Scale all elements by a scalar: <tt>y(i) *= alpha, i = 1...y->space()->dim()</tt>.
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>y = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 */
template<class Scalar>
void Vt_S( VectorBase<Scalar>* y, const Scalar& alpha );

/** \brief Assign scaled vector: <tt>y(i) = alpha * x(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void V_StV( VectorBase<Scalar>* y, const Scalar& alpha, const VectorBase<Scalar> &x );

/** \brief AXPY update: <tt>y(i) = alpha * x(i) + y(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void Vp_StV( VectorBase<Scalar>* y, const Scalar& alpha, const VectorBase<Scalar>& x );

/** \brief Scale and update: <tt>y(i) = x(i) + beta*y(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void Vp_V( VectorBase<Scalar>* y, const VectorBase<Scalar>& x, const Scalar& beta );

/** \brief Element-wise absolute value<tt>y(i) = abs(x(i)), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void abs( VectorBase<Scalar>* y, const VectorBase<Scalar>& x );

/** \brief Element-wise reciprocal: <tt>y(i) = 1/x(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void reciprocal( VectorBase<Scalar>* y, const VectorBase<Scalar>& x );

/** \brief Element-wise product update: <tt>y(i) += alpha * x(i) * v(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void ele_wise_prod( const Scalar& alpha, const VectorBase<Scalar>& x, const VectorBase<Scalar>& v, VectorBase<Scalar>* y );

/** \brief Element-wise division update: <tt>y(i) += alpha * x(i) / v(i), i = 1...y->space()->dim()</tt>.
 */
template<class Scalar>
void ele_wise_divide( const Scalar& alpha, const VectorBase<Scalar>& x, const VectorBase<Scalar>& v, VectorBase<Scalar>* y );

/** \brief Linear combination: <tt>y(i) = beta*y(i) + sum( alpha[k]*x[k](i), k=0...m-1 ), i = 1...y->space()->dim()</tt>.
 *
 * @param  m          [in] Number of vectors x[]
 * @param  alpha      [in] Array (length <tt>m</tt>) of input scalars.
 * @param  x          [in] Array (length <tt>m</tt>) of input vectors.
 * @param  beta       [in] Scalar multiplier for y
 * @param  y          [in/out] Target vector that is the result of the linear combination.
 *
 * This function implements a general linear combination:
 \verbatim
 y(i) = beta*y(i) + alpha[0]*x[0](i) + alpha[1]*x[1](i) + ... + alpha[m-1]*x[m-1](i), i = 1...y->space()->dim()
 \endverbatim
 */
template<class Scalar>
void linear_combination(
  const int                    m
  ,const Scalar                alpha[]
  ,const VectorBase<Scalar>*   x[]
  ,const Scalar                &beta
  ,VectorBase<Scalar>          *y
  );

/** \brief Seed the random number generator used in <tt>randomize()</tt>.
 *
 * @param  s  [in] The seed for the random number generator.
 *
 * Note, this just calls <tt>Teuchos::ScalarTraits<Scalar>::seedrandom(s)</tt>.
 */
template<class Scalar>
void seed_randomize( unsigned int s );

/** \brief Random vector generation: <tt>v(i) = rand(l,u), , i = 1...v->space()->dim()</tt>.
 * 
 * The elements <tt>v->getEle(i)</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using the above <tt>seed_randomize()</tt> function.
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, VectorBase<Scalar>* v );

//@}

/** \defgroup Thyra_Op_Vec_VectorStdOpsComparable_grp Subset of vector operations for types supporting relational operators.
 *
 * Warning! do not try to instantiate these functions for complex types
 * where relational operators are not defined (i.e. <, >, <=, >= etc).
 *
 * \ingroup Thyra_Op_Vec_VectorStdOps_grp
 */
//@{

/** \brief Min element: <tt>result = min{ x(i), i = 1...n } </tt>.
 */
template<class Scalar>
Scalar min( const VectorBase<Scalar>& x );

/** \brief Min element and its index: Returns <tt>maxEle = x(k)</tt>
 * and <tt>maxIndex = k</tt> such that <tt>x(k) <= x(i)</tt> for all
 * <tt>i=1...n</tt>.
 *
 * @param  x         [in] Input vector.
 * @param  minEle    [out] The minimum element value.
 * @param  maxindex  [out] The global index of the minimum element.
 *                   If there is more than one element with the
 *                   maximum entry then this returns the lowest index
 *                   in order to make the output independent of the order
 *                   of operations.
 *
 * Preconditions:<ul>
 * <li><tt>minEle!=NULL</tt>
 * <li><tt>minIndex!=NULL</tt>
 * </ul>
 */
template<class Scalar>
void min( const VectorBase<Scalar>& x, Scalar *maxEle, Index *maxIndex );

/** \brief Minimum element greater than some bound and its index:
 * Returns <tt>minEle = x(k)</tt> and <tt>minIndex = k</tt> such that
 * <tt>x(k) <= x(i)</tt> for all <tt>i</tt> where <tt>x(i) >
 * bound</tt>.
 *
 * @param  x         [in] Input vector.
 * @param  bound     [in] The upper bound
 * @param  minEle    [out] The minimum element value as defined above.
 * @param  minIndex  [out] The global index of the maximum element.
 *                   If there is more than one element with the
 *                   minimum value then this returns the lowest index
 *                   in order to make the output independent of the order
 *                   of operations.  If no entries are less than <tt>bound</tt>
 *                   then <tt>minIndex < 0</tt> on return.
 *
 * Preconditions:<ul>
 * <li><tt>minEle!=NULL</tt>
 * <li><tt>minIndex!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li>If <tt>*minIndex > 0</tt> then such an element was found.
 * <li>If <tt>*minIndex < 0</tt> then no such element was found.
 * </ul>
 */
template<class Scalar>
void minGreaterThanBound( const VectorBase<Scalar>& x, const Scalar &bound, Scalar *minEle, Index *minIndex );

/** \brief Max element: <tt>result = max{ x(i), i = 1...n } </tt>.
 */
template<class Scalar>
Scalar max( const VectorBase<Scalar>& x );

/** \brief Max element and its index: Returns <tt>maxEle = x(k)</tt>
 * and <tt>maxIndex = k</tt> such that <tt>x(k) >= x(i)</tt> for
 * <tt>i=1...n</tt>.
 *
 * @param  x         [in] Input vector.
 * @param  maxEle    [out] The maximum element value.
 * @param  maxindex  [out] The global index of the maximum element.
 *                   If there is more than one element with the
 *                   maximum value then this returns the lowest index
 *                   in order to make the output independent of the order
 *                   of operations.
 *
 * Preconditions:<ul>
 * <li><tt>maxEle!=NULL</tt>
 * <li><tt>maxIndex!=NULL</tt>
 * </ul>
 */
template<class Scalar>
void max( const VectorBase<Scalar>& x, Scalar *maxEle, Index *maxIndex );

/** \brief Max element less than bound and its index: Returns <tt>maxEle =
 * x(k)</tt> and <tt>maxIndex = k</tt> such that <tt>x(k) >= x(i)</tt> for all
 * <tt>i</tt> where <tt>x(i) < bound</tt>.
 *
 * @param  x         [in] Input vector.
 * @param  bound     [in] The upper bound
 * @param  maxEle    [out] The maximum element value as defined above.
 * @param  maxindex  [out] The global index of the maximum element.
 *                   If there is more than one element with the
 *                   maximum index then this returns the lowest index
 *                   in order to make the output independent of the order
 *                   of operations.  If no entries are less than <tt>bound</tt>
 *                   then <tt>minIndex < 0</tt> on return.
 *
 * Preconditions:<ul>
 * <li><tt>maxEle!=NULL</tt>
 * <li><tt>maxIndex!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li>If <tt>*maxIndex > 0</tt> then such an element was found.
 * <li>If <tt>*maxIndex < 0</tt> then no such element was found.
 * </ul>
 */
template<class Scalar>
void maxLessThanBound( const VectorBase<Scalar>& x, const Scalar &bound, Scalar *maxEle, Index *maxIndex );

//@}

} // end namespace Thyra

// /////////////////////////
// Inline functions

template<class Scalar>
inline
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::norm( const VectorBase<Scalar>& v )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ST::magnitude(ST::squareroot(v.space()->scalarProd(v,v)));
}


#endif // THYRA_VECTOR_STD_OPS_DECL_HPP
