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

#ifndef ANASAZI_RETURN_TYPE_HPP
#define ANASAZI_RETURN_TYPE_HPP

/*! \file AnasaziReturnType.hpp
  \brief Enumerated type used to pass back information from the eigensolver or linear algebra interfaces
*/

  /*!  \enum Anasazi::ReturnType
    
  \brief Enumerated type used to pass back information from the eigensolver or linear algebra interfaces.
  
  */
  
namespace Anasazi {

	enum ReturnType		 
	  {	
	    Ok,           /*!< Computation completed sucessfully */	
	    Undefined,    /*!< This operation is not defined */
	    Unconverged,  /*!< This operation returned unconverged */ 
	    Failed	  /*!< Any other numerical failure in the computation */
	  };
  
}

#endif
// end of file AnasaziReturnType.hpp
