/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_WrapMultiLevelPreconditioner.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 02/06|
 *----------------------------------------------------------------------*/
MOERTEL::ConstrainedPreconditioner::ConstrainedPreconditioner(
                          RefCountPtr<ML_Epetra::MultiLevelPreconditioner> mlprec,
                          RefCountPtr<Epetra_CrsMatrix> I,
                          RefCountPtr<Epetra_CrsMatrix> WT,
                          RefCountPtr<Epetra_CrsMatrix> B) :
mlprec_(mlprec),
I_(I),
WT_(WT),
B_(B)                          
{
  label_  = "MOERTEL::ConstrainedPreconditioner";
  return;
}

/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 02/06|
 *----------------------------------------------------------------------*/
int MOERTEL::ConstrainedPreconditioner::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  if (mlprec_==null)
  {
    cout << "MOERTEL: ***WRN*** MOERTEL::ConstrainedPreconditioner::ApplyInverse:\n"
         << "MOERTEL: ***WRN*** ML preconditioner is NULL\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  
  // make X (residual) satisfy constraints 
  Epetra_Vector x(View,X,0);
  Epetra_Vector xtmp(x.Map(),false);
  WT_->Multiply(false,x,xtmp);
  B_->Multiply(false,xtmp,xtmp);
  x.Update(-1.0,xtmp,1.0);
  
  // Apply ML
  int err = mlprec_->ApplyInverse(X,Y);
  if (err)
  {
    cout << "MOERTEL: ***WRN*** MOERTEL::ConstrainedPreconditioner::ApplyInverse:\n"
         << "MOERTEL: ***WRN*** ML preconditioner returned " << err << "\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  
  // make Y (search direction) satisfy constraints
  Epetra_Vector y(View,Y,0);
  Epetra_Vector& ytmp = xtmp;
  B_->Multiply(true,y,ytmp);
  WT_->Multiply(true,ytmp,ytmp);
  y.Update(-1.0,ytmp,1.0);
  
  return err;
}
