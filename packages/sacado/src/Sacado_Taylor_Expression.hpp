// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_TAYLOR_EXPRESSION_HPP
#define SACADO_TAYLOR_EXPRESSION_HPP

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace Taylor {

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT> 
    class Expr {

    public:

      //! Typename of values
      typedef typename ExprT::value_type value_type;

      //! Constructor with given expression \c expr
      explicit Expr(const ExprT& expr) : expr_(expr) {}

       //! Allocate coefficient cache
      void allocateCache(unsigned int d) const { expr_.allocateCache(d); }

      //! Return degree of polynomial
      unsigned int degree() const {return expr_.degree();}
      
      //! Return if expression has fast access
      bool hasFastAccess(unsigned int d) const { 
	return expr_.hasFastAccess(d);}

      //! Return degree \c i term of expression
      value_type coeff(unsigned int i) const { return expr_.coeff(i);}

      //! Return degree \c i term of expression
      value_type fastAccessCoeff(unsigned int i) const { return 
	  expr_.fastAccessCoeff(i);}
      
    protected:

      //! Disallow default constructor
      Expr() {}

      //! Expression
      ExprT expr_;

    }; // class Expr

    //! Constant expression template
    /*!
     * This template class represents a constant expression.
     */
    template <typename ConstT> 
    class ConstExpr {

    public:

      //! Typename of argument values
      typedef ConstT value_type;

      //! Constructor
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Allocate coefficient cache
      void allocateCache(unsigned int d) const {}

      //! Return degree of polynomial
      unsigned int degree() const { return 0; }
      
      //! Return if operation has fast access
      bool hasFastAccess(unsigned int d) const { return 1; }

      value_type value() const { return constant_; }

      //! Return degree \c i term of expression
      value_type coeff(unsigned int i) const { 
	return i==0 ? constant_ : value_type(0); }
      
      //! Return degree \c i term of expression
      value_type fastAccessCoeff(unsigned int i) const { 
	return i==0 ? constant_ : value_type(0); }

    protected:
      
      //! The constant
      ConstT constant_;

    }; // class ConstExpr

    //! Unary expression template
    /*!
     * This template class represents a unary operation of the form
     * op(a) where a is the argument of type \c ExprT and op is the 
     * operation represented by type \c Op. The operation is evaluated by the 
     * static methods Op::computeValue() and Op::computeDx().
     */
    template <typename ExprT, template<typename> class Op> 
    class UnaryExpr {

    public:

      //! Typename of argument value
      typedef typename ExprT::value_type value_type;

      //! Constructor
      UnaryExpr(const ExprT& expr) : expr_(expr), op_(expr) {}

      //! Allocate coefficient cache
      void allocateCache(unsigned int d) const { 
	expr_.allocateCache(d); 
	op_.allocateCache(d);
      }

      //! Return degree of polynomial
      unsigned int degree() const {return expr_.degree();}
      
      //! Return if operation has fast access
      bool hasFastAccess(unsigned int d) const { 
	return expr_.hasFastAccess(d); }

      //! Return degree \c i term of expression
      value_type coeff(unsigned int i) const { 
	return op_.computeCoeff(i,expr_); }
      
      //! Return derivative component \c i of operation
      value_type fastAccessCoeff(unsigned int i) const { 
	return op_.computeFastAccessCoeff(i,expr_); 
      }

    protected:
      
      //! Left argument
      ExprT expr_;

      //! Operator
      Op<ExprT> op_;

    }; // class UnaryExpr

    //! Binary expression template
    /*!
     * This template class represents a binary operation of the form
     * op(a1,a2) where a1 is the left argument of type \c ExprT1, r is 
     * the right argument of type \c ExprT2, and op is the operation 
     * represented by type \c Op. The operation is evaluated by the static 
     * methods Op::computeValue() and Op::computeDx().
     */
    template <typename ExprT1, typename ExprT2, 
	      template<typename,typename> class Op> 
    class BinaryExpr {
		
    public:

      //! Typename of the first argument value
      typedef typename ExprT1::value_type value_type_1;

      //! Typename of the second argument value
      typedef typename ExprT2::value_type value_type_2;

      //! Typename of the expression values
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      //! Constructor
      BinaryExpr(const ExprT1& expr1, const ExprT2& expr2) : 
	expr1_(expr1), expr2_(expr2), op_(expr1,expr2) {}

      //! Allocate coefficient cache
      void allocateCache(unsigned int d) const { 
	expr1_.allocateCache(d); 
	expr2_.allocateCache(d); 
	op_.allocateCache(d);
      }

      //! Return degree of polynomial
      unsigned int degree() const {
	unsigned int d1 = expr1_.degree(), d2 = expr2_.degree(); 
	return d1 > d2 ? d1 : d2;
      }
      
      //! Return if operation has fast access
      bool hasFastAccess(unsigned int d) const { 
	return expr1_.hasFastAccess(d) && expr2_.hasFastAccess(d);}

      //! Return degree \c i term of expression
      value_type coeff(unsigned int i) const { 
	return op_.computeCoeff(i,expr1_,expr2_); }
      
      //! Return degree \c i term of expression
      value_type fastAccessCoeff(unsigned int i) const { 
	return op_.computeFastAccessCoeff(i,expr1_,expr2_); 
      }

    protected:
      
      //! Left argument
      ExprT1 expr1_;

      //! Right argument
      ExprT2 expr2_;

      //! Operator
      Op<ExprT1,ExprT2> op_;

    }; // class BinaryExpr

    //! Binary expression template with first argument constant
    /*!
     * This template class represents a binary operation of the form
     * op(a1,a2) where a1 is the left argument of type \c ExprT1, r is 
     * the right argument of type \c ExprT2, and op is the operation 
     * represented by type \c Op. The operation is evaluated by the static 
     * methods Op::computeValue() and Op::computeDx().
     */
    template <typename ExprT2, template<typename,typename> class Op> 
    class BinaryExpr<ConstExpr<typename ExprT2::value_type>, ExprT2, Op> {
		
    public:

      //! Typename of constant expression
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      //! Typename of the second argument value
      typedef typename ExprT2::value_type value_type;

      //! Constructor
      BinaryExpr(const ExprT1& expr1, const ExprT2& expr2) : 
	expr1_(expr1), expr2_(expr2), op_(expr1,expr2) {}

      //! Allocate coefficient cache
      void allocateCache(unsigned int d) const { 
	expr2_.allocateCache(d); 
	op_.allocateCache(d);
      }

      //! Return degree of polynomial
      unsigned int degree() const {
	return expr2_.degree(); 
      }
      
      //! Return if operation has fast access
      bool hasFastAccess(unsigned int d) const { 
	return expr2_.hasFastAccess(d);}

      //! Return degree \c i term of expression
      value_type coeff(unsigned int i) const { 
	return op_.computeCoeff(i,expr1_,expr2_); }
      
      //! Return degree \c i term of expression
      value_type fastAccessCoeff(unsigned int i) const { 
	return op_.computeFastAccessCoeff(i,expr1_,expr2_); 
      }

    protected:
      
      //! Left argument
      ExprT1 expr1_;

      //! Right argument
      ExprT2 expr2_;

      //! Operator
      Op<ExprT1,ExprT2> op_;

    }; // class BinaryExpr

    //! Binary expression template with second argument constant
    /*!
     * This template class represents a binary operation of the form
     * op(a1,a2) where a1 is the left argument of type \c ExprT1, r is 
     * the right argument of type \c ExprT2, and op is the operation 
     * represented by type \c Op. The operation is evaluated by the static 
     * methods Op::computeValue() and Op::computeDx().
     */
    template <typename ExprT1, template<typename,typename> class Op> 
    class BinaryExpr<ExprT1,ConstExpr<typename ExprT1::value_type>, Op> {
		
    public:

      //! Typename of constant expression
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      //! Typename of the second argument value
      typedef typename ExprT1::value_type value_type;

      //! Constructor
      BinaryExpr(const ExprT1& expr1, const ExprT2& expr2) : 
	expr1_(expr1), expr2_(expr2), op_(expr1,expr2) {}

      //! Allocate coefficient cache
      void allocateCache(unsigned int d) const { 
	expr1_.allocateCache(d); 
	op_.allocateCache(d);
      }

      //! Return degree of polynomial
      unsigned int degree() const {
	return expr1_.degree(); 
      }
      
      //! Return if operation has fast access
      bool hasFastAccess(unsigned int d) const { 
	return expr1_.hasFastAccess(d);}

      //! Return degree \c i term of expression
      value_type coeff(unsigned int i) const { 
	return op_.computeCoeff(i,expr1_,expr2_); }
      
      //! eturn degree \c i term of expression
      value_type fastAccessCoeff(unsigned int i) const { 
	return op_.computeFastAccessCoeff(i,expr1_,expr2_); 
      }

    protected:
      
      //! Left argument
      ExprT1 expr1_;

      //! Right argument
      ExprT2 expr2_;

      //! Operator
      Op<ExprT1,ExprT2> op_;

    }; // class BinaryExpr

  } // namespace Taylor

} // namespace Sacado

#endif // SACADO_TAYLOR_EXPRESSION_HPP
