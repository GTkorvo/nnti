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

#ifndef TSFLINEARCOMBINATIONDECL_HPP
#define TSFLINEARCOMBINATIONDECL_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtendedOps
{
  /** 
   *
   */
  template <class Scalar> 
  class ConvertibleToVector
  {
  public:
    /** */
    virtual ~ConvertibleToVector(){;}

    /** */
    virtual Vector<Scalar> eval() const = 0 ;

    /** */
    VectorSpace<Scalar> space() const {return eval().space();}

    /** Return the dimension of the vector  */
    int dim() const {return eval().dim();} 

    /** 
     * Create a new vector that is a copy of this vector 
     */
    Vector<Scalar> copy() const {return eval().copy();}

    /** 
     * Element-by-element product (Matlab dot-star operator)
     */
    Vector<Scalar> dotStar(const Vector<Scalar>& other) const 
    {return eval().dotStar(other);}

    /** 
     * Element-by-element division (Matlab dot-slash operator)
     */
    Vector<Scalar> dotSlash(const Vector<Scalar>& other) const 
    {return eval().dotSlash(other);}

    /** 
     * Return element-by-element reciprocal as a new vector
     */
    Vector<Scalar> reciprocal() const {return reciprocal();}

    /** 
     * Return element-by-element absolute value as a new vector
     */
    Vector<Scalar> abs() const {return abs();} 

    /** */
    Scalar norm2() const {return eval().norm2();}

    /** */
    Scalar norm1() const {return eval().norm1();}

    /** */
    Scalar normInf() const {return eval().normInf();}

    /** */
    Scalar max() const {return eval().max();}

    /** */
    Scalar min() const {return eval().min();}

    /** Return the min element and the corresponding index */
    Scalar min(int& index)const {return eval().min(index);}

    /** Return the max element and the corresponding index */
    Scalar max(int& index)const {return eval().max(index);}
  };

  /** 
   * Class OpTimesLC holds an operator times something convertible to a vector
   */
  template <class Scalar, class Node>
  class OpTimesLC : public ConvertibleToVector<Scalar>
  {
  public:

    /** */
    virtual ~OpTimesLC(){;}

    /** */
    OpTimesLC(const Scalar& alpha, const Node& x);

    /** */
    OpTimesLC(const Scalar& alpha,
              const TSFExtended::LinearOperator<Scalar>& op, 
              const Node& x);

    /** 
     * Evaluate the term into the argument vector, overwriting 
     * the previous value of the argument. */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** Add the term into the argument vector */
    void addInto(TSFExtended::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** Evaluate the term and return its value */
    virtual TSFExtended::Vector<Scalar> eval() const ;

    /** Determine whether this term contains the given vector */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    /** */
    const LinearOperator<Scalar>& op() const {return op_;}

    /** */
    const Scalar& alpha() const {return alpha_;}

    /** */
    const Node& node() const {return x_;}

    /** */
    Scalar norm2() const {return eval().norm2();}

    /** */
    Scalar normInf() const {return eval().normInf();}
    
  private:
    Scalar alpha_;
    
    TSFExtended::LinearOperator<Scalar> op_;

    Node x_;

    /** */
    static Scalar one() {return Teuchos::ScalarTraits<Scalar>::one();}

    /** */
    static Scalar zero() {return Teuchos::ScalarTraits<Scalar>::zero();}
  };


  /**
   * Class LC2 is a 2-term linear combination
   */
  template <class Scalar, class Node1, class Node2>
  class LC2  : public ConvertibleToVector<Scalar>
  {
  public:
    /** */
    virtual ~LC2(){;}

    /** */
    LC2(const Node1& x1, const Node2& x2, LCSign sign = LCAdd);

    /** */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** */
    void addInto(TSFExtended::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** */
    virtual TSFExtended::Vector<Scalar> eval() const ;

    /** */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    
  private:
    Node1 x1_;

    Node2 x2_;

    LCSign sign_;

    /** */
    static Scalar one() {return Teuchos::ScalarTraits<Scalar>::one();}

    /** */
    static Scalar zero() {return Teuchos::ScalarTraits<Scalar>::zero();}
  };

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */


#endif
