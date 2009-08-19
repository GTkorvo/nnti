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

#ifndef TSF_COMMON_OPERATORS_IMPL_HPP
#define TSF_COMMON_OPERATORS_IMPL_HPP


#include "TSFCommonOperatorsDecl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "TSFLinearCombinationImpl.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#endif




namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;
using namespace SundanceUtils;
using std::endl;


/* ---- Simplified linear op with spaces ------- */

template <class Scalar> inline
SimplifiedLinearOpWithSpaces<Scalar>
::SimplifiedLinearOpWithSpaces(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : range_(range), domain_(domain) {}


template <class Scalar> inline
RCP< const VectorSpaceBase<Scalar> > 
SimplifiedLinearOpWithSpaces<Scalar>::range() const 
{
  return range_.ptr();
}


template <class Scalar> inline
RCP< const VectorSpaceBase<Scalar> > 
SimplifiedLinearOpWithSpaces<Scalar>::domain() const 
{
  return domain_.ptr();
}




/* ---- Identity op ------- */

template <class Scalar> inline
SimpleIdentityOp<Scalar>::SimpleIdentityOp(const VectorSpace<Scalar>& space)
  : SimplifiedLinearOpWithSpaces<Scalar>(space, space) {}


template <class Scalar> inline
void SimpleIdentityOp<Scalar>::applyOp(const Thyra::ETransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab;
  Out::os() << tab << "IdentityOp().apply()" << endl;
  Tabs tab1;
  Out::os() << tab1 << "in=" << in << endl;
      
  out.acceptCopyOf(in);
      
  Out::os() << tab1 << "out=" << out << endl;
}

template <class Scalar> inline
std::string SimpleIdentityOp<Scalar>::description() const 
{return "I(" + this->domain()->description() + ")";}


/* ---- Zero op ------- */

template <class Scalar> inline
SimpleZeroOp<Scalar>::SimpleZeroOp(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : SimplifiedLinearOpWithSpaces<Scalar>(domain, range) {}



template <class Scalar> inline
void SimpleZeroOp<Scalar>::applyOp(const Thyra::ETransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab;
  Out::os() << tab << "ZeroOp().apply()" << endl;
  Tabs tab1;
  Out::os() << tab1 << "in=" << in << endl;
      
  out.zero();
      
  Out::os() << tab1 << "out=" << out << endl;
}

/* */
template <class Scalar> inline
std::string SimpleZeroOp<Scalar>::description() const 
{return "ZeroOp(domain=" 
    + this->domain()->description() 
    + ", range=" + this->range()->description() + ")";}




/*
 * ------ composed operator  
 */
template <class Scalar> inline
SimpleComposedOp<Scalar>::SimpleComposedOp(const Array<LinearOperator<Scalar> >& ops)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    ops[ops.size()-1].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (unsigned int i=1; i<ops_.size(); i++)
  {
    TEST_FOR_EXCEPT(!(ops[i].range() == ops[i-1].domain()));
  }
}
  


template <class Scalar> inline
void SimpleComposedOp<Scalar>::applyOp(const Thyra::ETransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab;
  if (M_trans == Thyra::NOTRANS)
  {
    Out::os() << tab << this->description() << ".apply(): " << endl;
    Vector<Scalar> tmp = in.copy();
    Out::os() << "ops.size() = " << ops_.size() << endl;
    for (unsigned int i=0; i<ops_.size(); i++)
    {
      Tabs tab1;
      int j = ops_.size()-1-i;
      Out::os() << tab1 << "applying op #" << j ;
      Out::os() << tab1 << "op = " 
                << ops_[j].description() << endl;
      Out::os() << tab1 << "tmp in = " << tmp << endl;
      ops_[j].apply(tmp, tmp);
      Out::os() << tab1 << "tmp out = " << tmp << endl;
    }
    Out::os() << tab << "copying tmp vector into output" << endl;
    out.acceptCopyOf(tmp);
  }
  else if (M_trans == Thyra::TRANS)
  {
    Out::os() << tab << this->description() << ".applyTranspose(): " 
              << endl;
    Vector<Scalar> tmp = in.copy();
    for (unsigned int i=0; i<ops_.size(); i++)
    {
      Tabs tab1;
      Out::os() << tab1 << "applying transpose of " << ops_[i].description() << endl;
      Out::os() << tab1 << "tmp in = " << tmp << endl;
      ops_[i].applyTranspose(tmp, tmp);
      tmp = ops_[i] * tmp;
      Out::os() << tab1 << "tmp out = " << tmp << endl;
    }
    Out::os() << tab << "copying tmp vector into output" << endl;
    out.acceptCopyOf(tmp);
  }
  else
  {
    TEST_FOR_EXCEPT(M_trans != Thyra::TRANS && M_trans != Thyra::NOTRANS);
  }
}
  



template <class Scalar> inline
std::string SimpleComposedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (unsigned int i=0; i<ops_.size(); i++)
  {
    if (i > 0U) rtn += "*";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
}

/*
 * Represent a sum of operators A_0 + A_1 + ... + A_n.
 */
template <class Scalar> inline
SimpleAddedOp<Scalar>::SimpleAddedOp(const Array<LinearOperator<Scalar> >& ops)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    ops[0].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (unsigned int i=1; i<ops_.size(); i++)
  {
    TEST_FOR_EXCEPT(!(ops[i].range() == ops[0].range()));
    TEST_FOR_EXCEPT(!(ops[i].domain() == ops[0].domain()));
  }
}
  
/* */
template <class Scalar> inline
void SimpleAddedOp<Scalar>::applyOp(const Thyra::ETransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab;
  Out::os() << tab << this->description() << ".apply(): " << endl;
  Out::os() << tab << "in = " << in << endl;
  Vector<Scalar> tmp=out.copy();
  tmp.zero();
  for (unsigned int i=0; i<ops_.size(); i++)
  {
    Tabs tab1;
    Out::os() << tab1 << "applying " << ops_[i].description() << endl;
    if (M_trans == Thyra::NOTRANS)
      tmp = tmp + ops_[i] * in;
    else if (M_trans == Thyra::TRANS)
      tmp = tmp + ops_[i].transpose() * in;
    else 
      TEST_FOR_EXCEPT(M_trans != Thyra::TRANS && M_trans != Thyra::NOTRANS);
    Out::os() << tab1 << "tmp out = " << tmp << endl;
  }
  out.acceptCopyOf(tmp);
}
  
/* */
template <class Scalar> inline
std::string SimpleAddedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (unsigned int i=0; i<ops_.size(); i++)
  {
    if (i > 0U) rtn += "+";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
}



/*
 * --- transposed op
 */

template <class Scalar> inline
SimpleTransposedOp<Scalar>::SimpleTransposedOp(const LinearOperator<Scalar>& A)
  : SimplifiedLinearOpWithSpaces<Scalar>(
    A.range(), A.domain()
    ) 
  , A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleTransposedOp<Scalar>::applyOp(const Thyra::ETransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab;
  Out::os() << tab << this->description() << ".apply(): " << endl;
  Tabs tab1;
  Out::os() << tab1 << "in = " << in << endl;
  if (M_trans == Thyra::NOTRANS)
    A_.applyTranspose(in, out);
  else if (M_trans == Thyra::TRANS)
    A_.apply(in, out);
  else 
    TEST_FOR_EXCEPT(M_trans !=Thyra::TRANS && M_trans != Thyra::NOTRANS);
  Out::os() << tab1 << "out = " << out << endl;
}
  
/* */
template <class Scalar> inline
std::string SimpleTransposedOp<Scalar>::description() const 
{
  return "(" + A_.description() + "^T)";
}




}

#endif
