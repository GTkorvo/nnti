// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
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

#include "Range1D.hpp"
#include "Teuchos_TestForException.hpp"

namespace RangePack {

const Range1D Range1D::Invalid(Range1D::INVALID);

//#ifdef _DEBUG

void Range1D::assert_valid_range(int lbound, int ubound) const {
	TEST_FOR_EXCEPTION(
		lbound < 1, std::range_error
		,"Range1D::assert_valid_range(): Error, lbound ="<<lbound<<" must be greater than 0." );
	TEST_FOR_EXCEPTION(
		lbound > ubound, std::range_error
		,"Range1D::assert_valid_range(): Error, lbound = "<<lbound<<" > ubound = "<<ubound );
}

//#endif // _DEBUG

} // end namespace RangePack
