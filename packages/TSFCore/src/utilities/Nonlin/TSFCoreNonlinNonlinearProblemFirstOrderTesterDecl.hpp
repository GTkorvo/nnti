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

// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemFirstOrderTesterDecl.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_DECL_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_DECL_HPP

#include "TSFCoreNonlinTypes.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Unit testing class for <tt>NonlinearProblemFirstOrder</tt>.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NonlinearProblemFirstOrderTester {
public:

	///
	enum EOutputLevel { BASIC_OUTPUT, VERBOSE_OUTPUT };

	///
	/** Test a <tt>NonlinearProblemFirstOrder</tt> object.
	 *
	 * ToDo: Finish documentation!
	 */
	bool doTest( NonlinearProblemFirstOrder<Scalar>* prob, std::ostream* out = NULL, EOutputLevel outputLevel = BASIC_OUTPUT ) const;

}; // NonlinearProblemFirstOrderTester

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_DECL_HPP
