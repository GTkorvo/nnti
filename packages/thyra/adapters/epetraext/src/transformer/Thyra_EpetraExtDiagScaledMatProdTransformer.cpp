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


#include "Thyra_EpetraExtDiagScaledMatProdTransformer.hpp"
#include "Thyra_MultipliedLinearOpBase.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Overridden from LinearOpTransformerBase


bool EpetraExtDiagScaledMatProdTransformer::isCompatible(
      const LinearOpBase<double> &op_in) const
{
   TEST_FOR_EXCEPT(true);
   return false;
}


RCP<LinearOpBase<double> >
EpetraExtDiagScaledMatProdTransformer::createOutputOp() const
{
   return nonconstEpetraLinearOp();
}


void EpetraExtDiagScaledMatProdTransformer::transform(
   const LinearOpBase<double> &op_in,
   const Ptr<LinearOpBase<double> > &op_inout) const
{
   using Thyra::unwrap;
   using EpetraExt::MatrixMatrix;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::dyn_cast;
 
   //
   // A) Get the component Thyra objects for M = op(B) * D * G
   //
 
   const MultipliedLinearOpBase<double> &multi_op =
         dyn_cast<const MultipliedLinearOpBase<double> >(op_in);

   bool haveDiagScaling = (multi_op.numOps()==3);
 
   // get properties of first operator: Transpose, scaler multiply...etc
   const RCP<const LinearOpBase<double> > op_B = multi_op.getOp(0);
   double B_scalar = 0.0;
   EOpTransp B_transp = NOTRANS;
   RCP<const LinearOpBase<double> > B;
   unwrap( op_B, &B_scalar, &B_transp, &B );
   TEUCHOS_ASSERT(B_transp==NOTRANS || B_transp==CONJTRANS); // sanity check
 
   // get diagonal scaling
   RCP<const VectorBase<double> > d;
   if(haveDiagScaling) {
      const RCP<const LinearOpBase<double> > D = multi_op.getOp(1);
      d = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(D, true)->getDiag();
   }
 
   // get properties of third operator: Transpose, scaler multiply...etc
   const RCP<const LinearOpBase<double> > op_G = multi_op.getOp(haveDiagScaling ? 2 : 1);
   double G_scalar = 0.0;
   EOpTransp G_transp = NOTRANS;
   RCP<const LinearOpBase<double> > G;
   unwrap( op_G, &G_scalar, &G_transp, &G );
   TEUCHOS_ASSERT(G_transp==NOTRANS || G_transp==CONJTRANS); // sanity check

   //
   // B) Extract out the Epetra_CrsMatrix objects and the vector
   //
   
   // convert second operator to an Epetra_CrsMatrix
   const RCP<const Epetra_CrsMatrix> epetra_B =
         rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*B), true);
   // TEUCHOS_ASSERT( B_transp == NOTRANS ); // ToDo: Handle the transpose
   
   // extract dagonal
   RCP<const Epetra_Vector> epetra_d;
   if(haveDiagScaling) {
      epetra_d = (B_transp==CONJTRANS ? get_Epetra_Vector(epetra_B->OperatorRangeMap(), d)
                                      : get_Epetra_Vector(epetra_B->OperatorDomainMap(), d));
   }
 
   // convert third operator to an Epetra_CrsMatrix
   const RCP<const Epetra_CrsMatrix> epetra_G =
     rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*G), true);
   
   // determine row map for final operator
   const Epetra_Map op_inout_row_map 
         = (B_transp==CONJTRANS ? epetra_B->ColMap() : epetra_B->RowMap());
   const Epetra_Map op_inout_col_map 
         = (G_transp==CONJTRANS ? epetra_B->RowMap() : epetra_B->ColMap());
 
   //
   // C) Do the explicit multiplication
   //
  
   // allocate space for final product: 3 steps
   //   1. Get destination EpetraLinearOp
   //   2. Extract RCP to destination Epetra_CrsMatrix
   //   3. If neccessary, allocate new Epetra_CrsMatrix
   EpetraLinearOp &thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
   RCP<Epetra_CrsMatrix>  epetra_op =
         rcp_dynamic_cast<Epetra_CrsMatrix>(thyra_epetra_op_inout.epetra_op());
   if(is_null(epetra_op)) {
      epetra_op = Teuchos::rcp(
            new Epetra_CrsMatrix(::Copy, op_inout_row_map, 0));
   }
    
   // if necessary scale B by diagonal vector
   RCP<const Epetra_CrsMatrix> epetra_BD;
   if(haveDiagScaling) {
      // create a temporary to get around const issue
      RCP<Epetra_CrsMatrix> epetra_BD_temp = rcp(new Epetra_CrsMatrix(*epetra_B));

      // scale matrix depending on properties of B
      if(B_transp==CONJTRANS) 
         epetra_BD_temp->LeftScale(*epetra_d);
      else
         epetra_BD_temp->RightScale(*epetra_d);
 
      epetra_BD = epetra_BD_temp;
   }
   else
      epetra_BD = epetra_B;
 
   // perform multiply
   TEUCHOS_ASSERT(
       MatrixMatrix::Multiply( *epetra_BD,  B_transp==CONJTRANS,
                               *epetra_G,   G_transp==CONJTRANS, *epetra_op)==0);

   // scale the whole thing if neccessary
   if(B_scalar*G_scalar!=1.0)
      epetra_op->Scale(B_scalar*G_scalar);
 
   // set output operator to use newly create epetra_op
   thyra_epetra_op_inout.initialize(epetra_op);
}


} // namespace Thyra
