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
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SCALED_ADJOINT_LINEAR_OP_DECL_HPP
#define THYRA_SCALED_ADJOINT_LINEAR_OP_DECL_HPP

#include "Thyra_ScaledAdjointLinearOpBaseDecl.hpp"

namespace Thyra {

/** \brief Concrete decorator <tt>LinearOpBase</tt> subclass that wraps a
<tt>LinearOpBase</tt> object and adds on an extra scaling factor and/or a
transpose enum.

This class represents a scaled, adjointed (transposed) linear operator
<tt>M</tt> of the form:

\verbatim
 
  M = scalar * op(Op)
\endverbatim

where <tt>Op</tt> is another <tt>LinearOpBase</tt> object, <tt>scalar</tt> is
a <tt>Scalar</tt>, and the operation <tt>op(Op)</tt> is specified by a
<tt>ETransp</tt> and is given as <tt>op(Op) = Op</tt> (<tt>NOTRANS</tt>), or
<tt>op(Op) = Op^T</tt> (<tt>TRANS</tt>), or <tt>op(Op) = Op^H</tt>
(<tt>CONJTRANS</tt>).  Of course the operator <tt>M</tt> is not constructed
explicitly but instead just applies the decorated operator <tt>Op</tt> by
modifying the <tt>apply()</tt> function that calls <tt>Op.apply()</tt>.

This subclass is designed to allow the efficient handling of multiple implicit
scalings and/or adjoints (transposes) and allow these implicit
transformations to be reversed.  A sequence of scalings/adjoints from some
original <tt>LinearOpBase</tt> object <tt>origOp</tt> is shown as:

\verbatim

 M = scalar_n * op_n( ... scalar_2 * op_2( scalar_1 * op_1( origOp ) ) ... )

   =>

 M  = overallScalar * overall_op(origOp)
\endverbatim

where <tt>overallScalar</tt> and <tt>overall_op(...)</tt> are the cumulative
transformations given as:

\verbatim

 overallScalar = scalar_n * ... * scalar_2 * ... * scalar_1

 overall_op(origOp) = op_n( ... op_2( op_1( ... op_n( origOp ) ... ) ) )
\endverbatim
 
Each individual transformation pair <tt>(scalar_i,op_i(...))</tt> is specified
with arguments <tt>Scalar scalar</tt> and <tt>ETransp transp</tt>.  The
overall scaling is returned using <tt>this->overallScalar()</tt>, the overall
adjoint enum is returned using <tt>this->overallTransp()</tt>, and the
original linear operator is returned from <tt>this->getOrigOp()</tt>.

The operator most recently wrapped is returned from <tt>this->getOp()</tt>.
The individual scalings and transformations are not exposed from this
interface but can be viewed by calling <tt>this->description()</tt> with a
verbosity level of ???.  The arguments passed into the constructor
<tt>ScaledAdjointLinearOp()</tt> or <tt>initialize()</tt> can always be
extracted using <tt>uninitialize()</tt>.

This subclass keeps track of all of the individual scalings <tt>scalar_i</tt>
and adjoining operations <tt>op_i(...)</tt> so that the <tt>description()</tt>
function can print this out and also so that the operations can be reversed.

The copy constructor and assignment operators are declared private since some
thought needs to be put into what they should mean.

Note: This class does not maintain any specific cached information about the
original operator so it is safe, in some cases, to manipulate the original
operator and still keep <tt>this</tt> intact and automatically updated for any
changes that are made.

\ingroup Thyra_Op_Vec_ANA_Development_grp

*/
template<class Scalar>
class ScaledAdjointLinearOp : virtual public ScaledAdjointLinearOpBase<Scalar> {
public:

  /** \brief . */
  using ScaledAdjointLinearOpBase<Scalar>::apply;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  ScaledAdjointLinearOp();

  /** \brief Calls <tt>initialize()</tt>.
   *
   * Note, instead of calling this constructor directly consider using the
   * non-member functions described \ref Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp "here"
   * which create dynamically allocated <tt>Teuchos::RefCountPtr</tt> objects.
   */
  ScaledAdjointLinearOp(
    const Scalar                                               &scalar
    ,const ETransp                                             &transp
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &Op
    );

  /** \brief Initialize with an operator with by defining adjoint (transpose) and
   * scaling arguments.
   *
   * @param  scalar     [in] Scalar argument defining <tt>scalar_0</tt> (see introduction).
   * @param  transp     [in] Transpose argument defining <tt>op_0(...)</tt> (see introduction).
   * @param  Op         [in] Smart pointer to linear operator (persisting relationship).
   *
   * Preconditions:<ul>
   * <li><tt>Op.get() != NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li>ToDo: Fill these in!!!!
   * </ul>
   */
  void initialize(
    const Scalar                                               &scalar
    ,const ETransp                                             &transp
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &Op
    );

  /** \brief Set to uninitialized and (optionally) extract the objects passed into <tt>initialize()</tt>.
   *
   * @param  scalar     [out] Scalar argument passed into <tt>initialize()</tt>
   * @param  transp     [out] Transpose argument passed into <tt>initialize()</tt>.
   * @param  Op         [out] Smart pointer to linear operator passed into <tt>initialize()</tt>.
   *
   * Postconditions:<ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  void uninitialize(
    Scalar                                              *scalar  = NULL
    ,ETransp                                            *transp  = NULL
    ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  *Op      = NULL
    );

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Outputs
   * <tt>ScaledAdjointLinearOp<Scalar>{this->getOrigOp().description())</tt>
   * along with the dimensions.
   */
  std::string description() const;

  /** \brief Prints out the original operator as well as all of the scalings
   * and transpositions in the order that they occurred.
   *
   * This function outputs different levels of detail based on the value passed in
   * for <tt>verbLevel</tt>:
   *
   * ToDo: Finish documentation!
   */
  std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const;

  //@}

  /** @name Overridden from OpBase */
  //@{

  /** \brief Return the range space of the logical linear operator.
   *
   * Simply returns: \code

   return ( this->overallTransp()==NOTRANS ? this->getOrigOp()->range() : this->getOrigOp()->domain() );
   \endcode
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;

  /** \brief Return the domain space of the logical linear operator.
   *
   * Simply returns: \code

   return ( this->overallTransp()==NOTRANS ? this->getOrigOp()->domain() : this->getOrigOp()->range() );
   \endcode
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;

  /** Return if the operation is supported on the logical linear operator.
   *
   * Simply returns: \code

   return this->getOrigOp()->opSupported(trans_trans(this->overallTransp(),M_trans));
   \endcode
   */
  bool opSupported(ETransp M_trans) const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Apply the linear operator (or its transpose) to a multi-vector :
   * <tt>Y = alpha*op(M)*X + beta*Y</tt>.
   *
   * Simply calls: \code

   this->getOrigOp()->apply(trans_trans(M_trans,this->overallTransp()),X,Y,(this->overallScalar()*alpha),beta)
   \endcode
   */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

  /** \name Overridden from ScaledAdointLinearOpBase */
  //@{

  /** \brief . */
  Scalar overallScalar() const;

  /** \brief . */
  ETransp overallTransp() const;

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > getOrigOp() const;

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > getOp() const;

  //@}

private:

  // ////////////////////////////////
  // Private types

  template <class Scalar2>
  struct ScalarETransp {
    Scalar2   scalar;
    ETransp   transp;
    ScalarETransp()
      {}
    ScalarETransp( const Scalar2 &_scalar, const ETransp &_transp )
      : scalar(_scalar), transp(_transp)
      {}
  };

  typedef std::vector<ScalarETransp<Scalar> >  allScalarETransp_t;
  
  // ////////////////////////////////
  // Private data members

  Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  origOp_;
  Scalar                                             overallScalar_;
  ETransp                                            overallTransp_;
  int                                                my_index_;

  Teuchos::RefCountPtr<allScalarETransp_t>           allScalarETransp_;
  
  // ////////////////////////////////
  // Private member functions

  void assertInitialized() const;

};

/** \defgroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp Helper functions for creating scaled/adjoined linear operators.

These are helper functions for creating and manipulating
<tt>ScaledAdjointLinearOp</tt> objects.  In general, clients should create
dynamically allocated <tt>ScaledAdjointLinearOp</tt> objects by calling these
helper functions and not by directly calling the constructor for
<tt>ScaledAdjointLinearOp<tt>.

\ingroup Thyra_Op_Vec_ANA_Development_grp

*/

/** \brief Build an implicit non-<tt>const</tt> scaled linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
scale( const Scalar &scalar, const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit <tt>const</tt> scaled linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,rcp_const_cast<LinearOpBase<Scalar> >(Op))</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
scale( const Scalar &scalar, const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit non-<tt>const</tt> adjoined linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
adjoint( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit <tt>const</tt> adjoined linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,rcp_const_cast<LinearOpBase<Scalar> >(Op))</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
adjoint( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit non-<tt>const</tt> transposed linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),TRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
transpose( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit <tt>const</tt> transposed linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),TRANS,rcp_const_cast<LinearOpBase<Scalar> >(Op))</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
transpose( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit non-<tt>const</tt> scaled and/or adjoined (transposed) linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scale,transp,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
scaleAndAdjoint( const Scalar &scalar, const ETransp &transp, const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op );

/** \brief Build an implicit <tt>const</tt> scaled and/or adjoined (transposed) linear operator.
 *
 * Returns <tt>Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scale,transp,rcp_const_cast<LinearOpBase<Scalar> >(Op))</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
scaleAndAdjoint( const Scalar &scalar, const ETransp &transp, const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op );

// /////////////////////////////////
// Inline members

template<class Scalar>
inline
ScaledAdjointLinearOp<Scalar>::ScaledAdjointLinearOp()
  :overallScalar_(Teuchos::ScalarTraits<Scalar>::zero())
  ,overallTransp_(NOTRANS)
{}

template<class Scalar>
inline
ScaledAdjointLinearOp<Scalar>::ScaledAdjointLinearOp(
  const Scalar                                               &scalar
  ,const ETransp                                             &transp
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &Op
  )
  :overallScalar_(Teuchos::ScalarTraits<Scalar>::zero())
  ,overallTransp_(NOTRANS)
{
  this->initialize(scalar,transp,Op);
}

template<class Scalar>
inline
void ScaledAdjointLinearOp<Scalar>::assertInitialized() const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( origOp_.get() == NULL );
#endif
}

}	// end namespace Thyra

// /////////////////////////////////
// Inline non-members

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::scale( const Scalar &scalar, const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,Op));
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::scale( const Scalar &scalar, const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,Teuchos::rcp_const_cast<LinearOpBase<Scalar> >(Op)));
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::adjoint( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Op));
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::adjoint( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(
                        Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Teuchos::rcp_const_cast<LinearOpBase<Scalar> >(Op)));
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::transpose( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),TRANS,Op));
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::transpose( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(
                        Teuchos::ScalarTraits<Scalar>::one(),TRANS,Teuchos::rcp_const_cast<LinearOpBase<Scalar> >(Op)));
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::scaleAndAdjoint( const Scalar &scalar, const ETransp &transp, const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scalar,transp,Op));
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::scaleAndAdjoint( const Scalar &scalar, const ETransp &transp, const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &Op )
{
  return Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>(scalar,transp,Teuchos::rcp_const_cast<LinearOpBase<Scalar> >(Op)));
}

#endif	// THYRA_SCALED_ADJOINT_LINEAR_OP_DECL_HPP
