// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_MATRIX_HPP
#define ANASAZI_MATRIX_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziConfigDefs.hpp"

namespace Anasazi {

/*!	\class Anasazi::Matrix

	\brief Anasazi's templated pure virtual class for constructing the matrix/operator that is
	used by the eigensolver.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

template <class TYPE>
class Matrix {
public:

	//@{ \name Constructor/Destructor.
	//! %Anasazi::Matrix constructor.
	Matrix() {
//		cout << "ctor:Anasazi::Matrix " << this << endl; 
	}
	//! %Anasazi::Matrix destructor.
	virtual ~Matrix() {
//		cout << "dtor:Anasazi::Matrix " << this << endl; 
	};
	//@}
	
	//@{ \name Matrix/Operator application method.

	/*! \brief This routine takes the %Anasazi::MultiVec \c x and applies the matrix/operator
	to it resulting in the %Anasazi::MultiVec \c y, which is returned.
	*/
	virtual ReturnType ApplyMatrix (const MultiVec<TYPE>& x, 
						      MultiVec<TYPE>& y ) const = 0;
	//@}
};

} // end Anasazi namespace
#endif
// end of file AnasaziMatrix.hpp
