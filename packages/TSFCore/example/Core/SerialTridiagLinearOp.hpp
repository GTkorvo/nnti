// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

#ifndef TSFCORE_SERIAL_TRIDIAG_LINEAR_OP_HPP
#define TSFCORE_SERIAL_TRIDIAG_LINEAR_OP_HPP

#include "TSFCoreSerialLinearOpBase.hpp"

///
/** \brief Simple example subclass for serial tridiagonal matrices.
 *
 * This subclass form a linear operator for tridiagonal matrices
 * of the form:
 *
 \f[

 M=
 \left[\begin{array}{ccccc}
 d_{(1)} & u_{(1)}  \\
 l_{(1)} & d_{(2)} & u_{(2)} \\
         & \ddots  & \ddots    & \ddots \\
         &         & l_{(n-2)} & d_{(n-1)}  & u_{(n-1)} \\
         &         &           & l_{(n-1)}  & d_{(n)}
 \end{array}\right]
 \f]
 *
 * given the arrays <tt>lower[]</tt>, <tt>diag[]</tt>, and
 * <tt>upper[]</tt> of dimension <tt>dim-1</tt>, <tt>dim</tt> and <tt>dim-1</tt>
 * respectively (see <tt>initialize()</tt>).
 *
 * Note, this is just an exmaple class and for the sake if simplifed
 * presentation the private members are listed first and inclass
 * declarations are used which are not a good idea in code.  However,
 * in this case, they make the example code easier to read and
 * maintaining encapsulation and a well defined interface are
 * unecessary here.
 *
 * See the source code for this simple example by clicking on the
 * link to the defintion below.
 *
 * \ingroup TSFCore_examples_power_method_serial_grp
 * \ingroup TSFCore_examples_cg_serial_grp
 */
template<class Scalar>
class SerialTridiagLinearOp : virtual public TSFCore::SerialLinearOpBase<Scalar> {
private:

	TSFCore::Index       dim_;
	std::vector<Scalar>  lower_;   // size dim - 1    
	std::vector<Scalar>  diag_;    // size dim
	std::vector<Scalar>  upper_;   // size dim - 1

public:

	/// Calls <tt>initialize()</tt>.
	SerialTridiagLinearOp( const TSFCore::Index dim, const Scalar lower[], const Scalar diag[], const Scalar upper[] )
		{ this->initialize(dim,lower,diag,upper);	}
	
	///
	/** Initialize given lower, diagonal and upper arrays of data.
	 *
	 * @param  dim    [in] Dimension of this matrix (must be >= 2).
	 * @param  lower  [in] Array (length <tt>dim-1</tt>) of the lower diagonal elements
	 * @param  diag   [in] Array (length <tt>dim</tt>) of the central diagonal elements
	 * @param  upper  [in] Array (length <tt>dim-1</tt>) of the upper diagonal elements
	 *
	 * Preconditions:<ul>
	 * <li><tt>dim >= 2</tt>
	 * </ul>
	 *
	 * Postonditions:<ul>
	 * <li>Should be obvious!
	 * <li>See <tt>setDimensions()</tt>
	 * </ul>
	 */
	void initialize(
		const TSFCore::Index            dim      // >= 2
		,const Scalar                   lower[]  // size == dim - 1
		,const Scalar                   diag[]   // size == dim
		,const Scalar                   upper[]  // size == dim - 1
		)
		{
			TEST_FOR_EXCEPT( dim < 2 );
			setDimensions(dim,dim); // We must tell the base class our dimension to setup range() and domain()
			dim_ = dim;
			lower_.resize(dim-1);  for( int k = 0; k < dim-1; ++k ) lower_[k] = lower[k];
			diag_.resize(dim);     for( int k = 0; k < dim;   ++k ) diag_[k]  = diag[k];
			upper_.resize(dim-1);  for( int k = 0; k < dim-1; ++k ) upper_[k] = upper[k];
		}

protected:

	// Overridden from SeiralLinearOpBase

	void euclideanApply(
		const TSFCore::ETransp                       M_trans
		,const RTOpPack::SubVectorT<Scalar>          &x_in
		,const RTOpPack::MutableSubVectorT<Scalar>   *y_out
		,const Scalar                                alpha
		,const Scalar                                beta
		) const
		{
			typedef Teuchos::ScalarTraits<Scalar> ST;
			// Handle transpose
			const std::vector<Scalar>
				&lower = ( M_trans == TSFCore::NOTRANS ? lower_ : upper_ ),
				&upper = ( M_trans == TSFCore::NOTRANS ? upper_ : lower_ );
			// Get raw pointers to vector data to make me feel better!
			const Scalar *x = x_in.values();
			Scalar       *y = y_out->values();
			// Perform operation (if beta==0 then we must be careful since y could be uninitialized on input!)
			TSFCore::Index k = 0;
			if( beta == ST::zero() ) {
				y[k] = alpha * ( diag_[k]*x[k] + upper[k]*x[k+1] );                                    // First row
				for( k = 1; k < dim_ - 1; ++k )
					y[k] = alpha * ( lower[k-1]*x[k-1] + diag_[k]*x[k] + upper[k]*x[k+1] );              // Middle rows
				y[k] = alpha * ( lower[k-1]*x[k-1] + diag_[k]*x[k] );                                  // Last row
			}
			else {
				y[k] = alpha * ( diag_[k]*x[k] + upper[k]*x[k+1] ) + beta*y[k];                        // First row
				for( k = 1; k < dim_ - 1; ++k )
					y[k] = alpha * ( lower[k-1]*x[k-1] + diag_[k]*x[k] + upper[k]*x[k+1] ) + beta*y[k];  // Middle rows
				y[k] = alpha * ( lower[k-1]*x[k-1] + diag_[k]*x[k] ) + beta*y[k];                      // Last row
			}
		}

};	// end class SerialTridiagLinearOp

#endif	// TSFCORE_SERIAL_TRIDIAG_LINEAR_OP_HPP
