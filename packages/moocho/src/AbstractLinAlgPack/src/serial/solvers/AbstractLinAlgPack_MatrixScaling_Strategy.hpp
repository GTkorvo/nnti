// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef MATRIX_SCALING_STRATEGY_H
#define MATRIX_SCALING_STRATEGY_H

namespace AbstractLinAlgPack {

/** \brief Abstract interface for sparse matrix scaling strategies
 *
 * ToDo: Finish documentation!
 */
class MatrixScaling_Strategy {
public:

  /** \brief . */
  virtual ~MatrixScaling_Strategy() {}

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /// Scale the matrix and save the scalings for later use for rhs and lhs.
  virtual void scale_matrix(
    index_type m, index_type n, index_type nz
    ,const index_type row_i[], const index_type col_j[]
    ,bool new_matrix, value_type A[]
    ) = 0;
  
  /// Scale the rhs vector
  virtual void scale_rhs( BLAS_Cpp::Transp trans, value_type b[] ) const = 0;
  
  /// Scale the lhs vector
  virtual void scale_lhs( BLAS_Cpp::Transp trans, value_type x[] ) const = 0;
  
  //@}

}; // end class MatrixScaling_Strategy

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SCALING_STRATEGY_H
