// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef __sun

#include "Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp"
#include "Thyra_DiagonalLinearOp.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace Thyra {

bool DiagonalEpetraLinearOpWithSolveFactory::isCompatible(
  const LinearOpBase<double> &fwdOp
  ) const
{
  const EpetraLinearOp *eFwdOp = NULL;
  if( ! (eFwdOp = dynamic_cast<const EpetraLinearOp*>(&fwdOp)) )
    return false;
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*eFwdOp->epetra_op()) )
    return false;
  return true;
}

Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
DiagonalEpetraLinearOpWithSolveFactory::createOp() const
{
  return Teuchos::rcp(new DiagonalLinearOp<double>());
}

void DiagonalEpetraLinearOpWithSolveFactory::initializeOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,LinearOpWithSolveBase<double>                             *Op
  ) const
{
  TEST_FOR_EXCEPT(Op==NULL);
  const EpetraLinearOp   &eFwdOp = Teuchos::dyn_cast<const EpetraLinearOp>(*fwdOp);
  const Epetra_RowMatrix &eRMOp  = Teuchos::dyn_cast<const Epetra_RowMatrix>(*(eFwdOp.epetra_op()));
  const Epetra_Map &map = eRMOp.OperatorDomainMap();
  Teuchos::RefCountPtr<Epetra_Vector>
    e_diag = Teuchos::rcp(new Epetra_Vector(map));
  eRMOp.ExtractDiagonalCopy(*e_diag);
  Teuchos::RefCountPtr< const MPIVectorSpaceBase<double> >
    space = create_MPIVectorSpaceBase(Teuchos::rcp(new Epetra_Map(map)));
  Teuchos::RefCountPtr< const MPIVectorBase<double> >
    diag = create_MPIVectorBase(e_diag,space);
  Teuchos::set_extra_data<Teuchos::RefCountPtr<const LinearOpBase<double> > >(
    fwdOp, "Thyra::DiagonalEpetraLinearOpWithSolveFactory::fwdOp", &diag
    );
  Teuchos::dyn_cast< DiagonalLinearOp<double> >(*Op).initialize(diag);
}

void DiagonalEpetraLinearOpWithSolveFactory::uninitializeOp(
  LinearOpWithSolveBase<double>                          *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >     *fwdOp
  ) const
{
  using Teuchos::get_extra_data;
  TEST_FOR_EXCEPT(Op==NULL);
  DiagonalLinearOp<double> &diagOp = Teuchos::dyn_cast<DiagonalLinearOp<double> >(*Op);
  Teuchos::RefCountPtr< const VectorBase<double> > diag = diagOp.diag();
  if( fwdOp ) {
    if(diag.get()) {
      *fwdOp =
        get_extra_data<Teuchos::RefCountPtr<const LinearOpBase<double> > >(
          diag,"Thyra::DiagonalEpetraLinearOpWithSolveFactory::fwdOp"
          );
    }
  }
  else {
    *fwdOp = Teuchos::null;
  }
}

} // namespace Thyra

#endif // __sun










