// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_FULLY_ASSEMBLED_OPERATOR_HPP
#define STOKHOS_FULLY_ASSEMBLED_OPERATOR_HPP

#include "Stokhos_SGOperator.hpp"
#include "EpetraExt_BlockCrsMatrix.h"
#include <vector>
#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the block stochastic Galerkin
   * operator generated by fully assembling the matrix.
   */
  class FullyAssembledOperator : 
    public Stokhos::SGOperator, 
    public EpetraExt::BlockCrsMatrix {
      
  public:

    //! Constructor 
    FullyAssembledOperator(
      const Teuchos::RCP<const Epetra_CrsMatrix>& base_matrix,
      const Teuchos::RCP<const std::vector< std::vector<int> > >& rowStencil,
      const Teuchos::RCP<const std::vector<int> >& rowIndex,
      const Teuchos::RCP<const Epetra_Comm>& sg_comm,
      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    
    //! Destructor
    virtual ~FullyAssembledOperator();

    /** \name Stokhos::SGOperator methods */
    //@{

    //! Setup operator
    virtual void setupOperator(
      const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& poly,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk);

    //! Get SG polynomial
    virtual Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > 
    getSGPolynomial();

    //! Get SG polynomial
    virtual Teuchos::RCP<const Stokhos::VectorOrthogPoly<Epetra_Operator> > 
    getSGPolynomial() const;

    //! Get triple product tensor
    virtual Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
    getTripleProduct() const;

    //@}

  private:
    
    //! Private to prohibit copying
    FullyAssembledOperator(const FullyAssembledOperator&);
    
    //! Private to prohibit copying
    FullyAssembledOperator& operator=(const FullyAssembledOperator&);
    
  protected:

    //! Stores triple product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;

    //! Stores operators
    Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> > block_ops;

    //! Flag indicating whether operator be scaled with <\psi_i^2>
    bool scale_op;

    //! Flag indicating whether to include mean term
    bool include_mean;
    
    //! Flag indicating whether to only use linear terms
    bool only_use_linear;

  }; // class FullyAssembledOperator
  
} // namespace Stokhos

#endif // STOKHOS_FULLY_ASSEMBLED_OPERATOR_HPP
