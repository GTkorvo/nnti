// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

#ifndef STOKHOS_MATRIX_FREE_EPETRA_OP_HPP
#define STOKHOS_MATRIX_FREE_EPETRA_OP_HPP

#include "Teuchos_RCP.hpp"

#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the block stochastic Galerkin
   * operator.
   */
  class MatrixFreeEpetraOp : public Epetra_Operator {
      
  public:

    //! Constructor 
    MatrixFreeEpetraOp(
     const Teuchos::RCP<const Epetra_Map>& base_map,
     const Teuchos::RCP<const Epetra_Map>& sg_map,
     const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
     const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
     const Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Operator> > >& ops);
    
    //! Destructor
    virtual ~MatrixFreeEpetraOp();

    //! Get operator blocks
    virtual std::vector< Teuchos::RCP<Epetra_Operator> >&
    getOperatorBlocks();
    
    //! Set to true if the transpose of the operator is requested
    virtual int SetUseTranspose(bool UseTranspose);
    
    /*! 
     * \brief Returns the result of a Epetra_Operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int Apply(const Epetra_MultiVector& Input, 
                      Epetra_MultiVector& Result) const;

    /*! 
     * \brief Returns the result of the inverse of the operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;
    
    //! Returns an approximate infinity norm of the operator matrix.
    virtual double NormInf() const;
    
    //! Returns a character string describing the operator
    virtual const char* Label () const;
  
    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;
    
    /*! 
     * \brief Returns true if the \e this object can provide an 
     * approximate Inf-norm, false otherwise.
     */
    virtual bool HasNormInf() const;

    /*! 
     * \brief Returns a reference to the Epetra_Comm communicator 
     * associated with this operator.
     */
    virtual const Epetra_Comm & Comm() const;

    /*!
     * \brief Returns the Epetra_Map object associated with the 
     * domain of this matrix operator.
     */
    virtual const Epetra_Map& OperatorDomainMap () const;

    /*! 
     * \brief Returns the Epetra_Map object associated with the 
     * range of this matrix operator.
     */
    virtual const Epetra_Map& OperatorRangeMap () const;

  private:
    
    //! Private to prohibit copying
    MatrixFreeEpetraOp(const MatrixFreeEpetraOp&);
    
    //! Private to prohibit copying
    MatrixFreeEpetraOp& operator=(const MatrixFreeEpetraOp&);
    
  protected:
    
    //! Label for operator
    string label;
    
    //! Stores base map
    Teuchos::RCP<const Epetra_Map> base_map;

    //! Stores SG map
    Teuchos::RCP<const Epetra_Map> sg_map;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Stores triple product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;

    //! Stores operators
    Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Operator> > > block_ops;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Number of blocks
    unsigned int num_blocks;

    //! BlockMultiVector for Apply() input
    mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_input;

    //! BlockMultiVector for Apply() result
    mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_result;

    //! MultiVectors for each block for Apply() input
    mutable std::vector< Teuchos::RCP<Epetra_MultiVector> > input_block;

    //! MultiVectors for each block for Apply() result
    mutable std::vector< Teuchos::RCP<Epetra_MultiVector> > result_block;

    //! Temporary multivector
    mutable Teuchos::RCP<Epetra_MultiVector> tmp;

  }; // class MatrixFreeEpetraOp
  
} // namespace Stokhos

#endif // STOKHOS_MATRIX_FREE_EPETRA_OP_HPP
