/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
 /* @HEADER@ */

#ifndef TSFEPETRAMATRIX_HPP
#define TSFEPETRAMATRIX_HPP

#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraMatrixFactory.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFEuclideanOpWithBackwardsCompatibleApply.hpp"
#include "TSFRowAccessibleOp.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFILUFactorizableOp.hpp"
#include "Epetra_CrsMatrix.h"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"

namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

class EpetraMatrix : virtual public EuclideanOpWithBackwardsCompatibleApply<double>,
                     public LoadableMatrix<double>,
                     public RowAccessibleOp<double>,
                     public Printable,
                     public ILUFactorizableOp<double>,
                     virtual public EpetraLinearOpBase
{
public:

  /** Construct an empty EpetraMatrix structured according to the graph 
   * argument */
  EpetraMatrix(const Epetra_CrsGraph& graph,
    const RefCountPtr<const EpetraVectorSpace>& domain,
    const RefCountPtr<const EpetraVectorSpace>& range);

  /** Wrap an existing Epetra CRS Matrix */
  EpetraMatrix(const RefCountPtr<Epetra_CrsMatrix>& mat,
    const RefCountPtr<const EpetraVectorSpace>& domain,
    const RefCountPtr<const EpetraVectorSpace>& range);

  /** */
  RefCountPtr< const VectorSpaceBase<double> > domain() const {return domain_;}

  /** */
  RefCountPtr< const VectorSpaceBase<double> > range() const {return range_;}


  /**
   * generalApply() applies either the operator or the transpose
   * according to the value of the transpose flag. This method is
   * backwards compatible with TSFCore-based code.
   */
  virtual void generalApply(const Thyra::ETransp M_trans,
    const Thyra::VectorBase<double>    &x,
    Thyra::VectorBase<double>          *y,
    const double            alpha,
    const double            beta) const ;


  virtual bool applyTransposeSupports(const Thyra::EConj conj) const
    {
      if (conj==Thyra::NONCONJ_ELE) return true;
      return false;
    }


  /** \name Epetra-Thyra adapter interface */
  //@{

  /** */
  void getEpetraOpView(RefCountPtr<Epetra_Operator> *epetraOp,
    Thyra::ETransp *epetraOpTransp,
    Thyra::EApplyEpetraOpAs *epetraOpApplyAs,
    Thyra::EAdjointEpetraOp *epetraOpAdjointSupport) ;

  /** */
  void getEpetraOpView(RefCountPtr<const Epetra_Operator> *epetraOp,
    Thyra::ETransp *epetraOpTransp,
    Thyra::EApplyEpetraOpAs *epetraOpApplyAs,
    Thyra::EAdjointEpetraOp *epetraOpAdjointSupport) const ;


  /// Returns <tt>this->mpiRange()</tt>
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<double> > rangeScalarProdVecSpc() const;
  /// Returns <tt>this->mpiDomain()</tt>
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<double> > domainScalarProdVecSpc() const;


  //@}

  /** \brief . */
	bool opSupported(Thyra::ETransp M_trans) const;

  /** Insert a set of elements in a row, adding to any previously
   * existing values. 
   * @param globalRowIndex the global index of the row to which these
   * elements belong.
   * @param nElemsToInsert the number of elements being inserted in this
   * step
   * @param globalColumnIndices array of column indices. Must 
   * be nElemsToInsert in length. 
   * @param elements array of element values. Must be nElemsToInsert in
   * length
   */
  virtual void addToRow(int globalRowIndex,
    int nElemsToInsert,
    const int* globalColumnIndices,
    const double* elementValues) ;


  /** 
   * Add to a batch of elements
   */
  virtual void addToElementBatch(int numRows, 
    int rowBlockSize,
    const int* globalRowIndices,
    int numColumnsPerRow,
    const int* globalColumnIndices,
    const double* values,
    const int* skipRow);

  /** Set all elements to zero, preserving the existing structure */
  virtual void zero() ;


  /** \name incomplete factorization preconditioning interface */
  //@{
  /** create an incomplete factorization. 
   * @param fillLevels number of levels of fill on the local processor
   * @param overlapFill number of levels of fill on remote processors
   * @param relaxationValue fraction of dropped values to be added to the
   * diagonal
   * @param relativeThreshold relative diagonal perutrbation
   * @param absoluteThreshold absolute diagonal perturbation
   * @param leftOrRight whether this preconditioner is to be applied
   * from the left or right 
   * @param rtn newly created preconditioner, returned 
   * by reference argument.
   */
  virtual void getILUKPreconditioner(int fillLevels,
    int overlapFill,
    double relaxationValue,
    double relativeThreshold,
    double absoluteThreshold,
    LeftOrRight leftOrRight,
    Preconditioner<double>& rtn) const ;

  /** Printable interface */
  virtual void print(ostream& os) const ;


  /** */
  std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    , const std::string                   indentSpacer
    ) const 
    {
      out << leadingIndent << indentSpacer << this->description() << std::endl;
      return out;
    }
    

  /** */
  string description() const ;

  /** */
  static Epetra_CrsMatrix& getConcrete(const LinearOperator<double>& A);

  /** */
  static RefCountPtr<const Epetra_CrsMatrix> getConcretePtr(const LinearOperator<double>& A);

  /** 
   * Read-only access to the underlying crs matrix. Needed for Ifpack.
   */
  const Epetra_CrsMatrix* crsMatrix() const ;

protected:

  /** Get the specified row as defined by RowAccessible  */
  void getRow(const int& row, 
		Teuchos::Array<int>& indices, 
		Teuchos::Array<double>& values) const;

private:

  Epetra_CrsMatrix* crsMatrix();

  RefCountPtr<Epetra_CrsMatrix> matrix_;

  RefCountPtr<const VectorSpaceBase<double> > range_;

  RefCountPtr<const VectorSpaceBase<double> > domain_;

  const Epetra_Map& getRangeMap() const;
  const Epetra_Map& getDomainMap() const;
};
}

#endif
