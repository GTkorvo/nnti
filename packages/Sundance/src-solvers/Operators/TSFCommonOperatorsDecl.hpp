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

#ifndef TSF_COMMON_OPERATORS_DECL_HPP
#define TSF_COMMON_OPERATORS_DECL_HPP


#include "SundanceDefs.hpp"
#include "TSFSimplifiedLinearOpBase.hpp"
#include "TSFLinearOperatorDecl.hpp"


namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;
using namespace SundanceUtils;
using std::endl;

/** */
template <class Scalar> 
class SimplifiedLinearOpWithSpaces
  : public virtual SimplifiedLinearOpBase<Scalar>,
    public virtual LinearOpBase<Scalar,Scalar>
{
public:
  /** */
  SimplifiedLinearOpWithSpaces(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range);

  /** 
   * \brief Return a smart pointer for the range space 
   * for <tt>this</tt> operator.
   */
  RCP< const VectorSpaceBase<Scalar> > range() const ;

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   */
  RCP< const VectorSpaceBase<Scalar> > domain() const ;

private:
  VectorSpace<Scalar> range_;
  VectorSpace<Scalar> domain_;
  
};


/** */
template <class Scalar>
class SimpleIdentityOp : public SimplifiedLinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleIdentityOp(const VectorSpace<Scalar>& space);


  /** */
  void applyOp(const Thyra::ETransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const ;

  /** */
  std::string description() const ;
};


/** */
template <class Scalar>
class SimpleZeroOp : public SimplifiedLinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleZeroOp(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range);

  /** */
  void applyOp(const Thyra::ETransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  /** */
  std::string description() const ;
};



/**
 * Represent a composed operator A_0 * A_1 * ... * A_n.
 */
template <class Scalar>
class SimpleComposedOp : public SimplifiedLinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleComposedOp(const Array<LinearOperator<Scalar> >& ops);
  
  /** */
  void applyOp(const Thyra::ETransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  
  /** */
  std::string description() const ;

private:
    Array<LinearOperator<Scalar> > ops_;
};

/**
 * Represent a sum of operators A_0 + A_1 + ... + A_n.
 */
template <class Scalar>
class SimpleAddedOp : public SimplifiedLinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleAddedOp(const Array<LinearOperator<Scalar> >& ops);

  /** */
  void applyOp(const Thyra::ETransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;
  
  /** */
  std::string description() const ;

private:
    Array<LinearOperator<Scalar> > ops_;
};


/**
 * Represent the transpose of an operator
 */
template <class Scalar>
class SimpleTransposedOp : public SimplifiedLinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleTransposedOp(const LinearOperator<Scalar>& A);
  
  /** */
  void applyOp(const Thyra::ETransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  
  /** */
  std::string description() const ;

  /** */
  LinearOperator<Scalar> op() const {return A_;}

private:
    LinearOperator<Scalar> A_;
};




}

#endif
