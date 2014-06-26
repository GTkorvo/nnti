// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_InitPostOps.hpp"
#include "Teuchos_Assert.hpp"
//#include "Teuchos_Exceptions.hpp"
#include "Epetra_Map.h"
#include "EpetraExt_MatrixMatrix.h"

#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"

FEApp::ResidualOp::ResidualOp(
        const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
	      const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
	      const Teuchos::RCP<Epetra_Vector>& overlapped_f) :
  xdot(overlapped_xdot),
  x(overlapped_x),
  f(overlapped_f)
{
}

FEApp::ResidualOp::~ResidualOp()
{
}

void
FEApp::ResidualOp::elementInit(const FEApp::AbstractElement& e,
                               unsigned int neqn,
                               std::vector<double>* elem_xdot,
                               std::vector<double>& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = x->Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = x->Map().LID(static_cast<int>(node_GID*neqn));
#endif
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = (*x)[firstDOF+j];
      if (elem_xdot != NULL)
        (*elem_xdot)[neqn*i+j] = (*xdot)[firstDOF+j];
    }
  }
}

void
FEApp::ResidualOp::elementPost(const FEApp::AbstractElement& e,
                               unsigned int neqn,
                               std::vector<double>& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual into global residual
  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
      f->SumIntoGlobalValue(row, 0, elem_f[lrow]);
    }
  }
}

void
FEApp::ResidualOp::nodeInit(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector<double>* node_xdot,
                            std::vector<double>& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = x->Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = x->Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy node solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = (*x)[firstDOF+j];
    if (node_xdot != NULL)
      (*node_xdot)[j] = (*xdot)[firstDOF+j];
  }
}

void
FEApp::ResidualOp::nodePost(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector<double>& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  // Replace global residual
  for (unsigned int j=0; j<offsets.size(); j++) {
    int row = static_cast<int>(firstDOF + offsets[j]);
    if (bc.isOwned())
      f->ReplaceGlobalValue(row, 0, node_f[offsets[j]]);
    else if (bc.isShared())
      f->ReplaceGlobalValue(row, 0, 0.0);
  }
}

FEApp::JacobianOp::JacobianOp(
             double alpha, double beta,
             const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
             const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
             const Teuchos::RCP<Epetra_Vector>& overlapped_f,
             const Teuchos::RCP<Epetra_CrsMatrix>& overlapped_jac) :
  m_coeff(alpha),
  j_coeff(beta),
  xdot(overlapped_xdot),
  x(overlapped_x),
  f(overlapped_f),
  jac(overlapped_jac)
{
}

FEApp::JacobianOp::~JacobianOp()
{
}

void
FEApp::JacobianOp::elementInit(const FEApp::AbstractElement& e,
                               unsigned int neqn,
                               std::vector< FadType >* elem_xdot,
                               std::vector< FadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = x->Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = x->Map().LID(static_cast<int>(node_GID*neqn));
#endif
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = 
        FadType(ndof, (*x)[firstDOF+j]);
      elem_x[neqn*i+j].fastAccessDx(neqn*i+j) = j_coeff;
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j] = 
          FadType(ndof, (*xdot)[firstDOF+j]);
        (*elem_xdot)[neqn*i+j].fastAccessDx(neqn*i+j) = m_coeff;
      }
    }
  }
}

void
FEApp::JacobianOp::elementPost(
			    const FEApp::AbstractElement& e,
			    unsigned int neqn,
			    std::vector< FadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  // Loop over nodes in element
  int row, col;
  unsigned int lrow, lcol;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
        f->SumIntoGlobalValue(row, 0, elem_f[lrow].val());
	
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {

        // Loop over nodes in element
        for (unsigned int node_col=0; node_col<nnode; node_col++){
          
          // Loop over equations per node
          for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
            lcol = neqn*node_col+eq_col;
            
            // Global column
            col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);
            
            // Sum Jacobian
            jac->SumIntoGlobalValues(row, 1, 
                                     &(elem_f[lrow].fastAccessDx(lcol)),
                                     &col);
            
          } // column equations
          
        } // column nodes
        
      } // has fast access
      
    } // row equations
    
  } // row node
  
}

void
FEApp::JacobianOp::nodeInit(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector< FadType >* node_xdot,
                            std::vector< FadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = x->Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = x->Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = FadType(neqn, (*x)[firstDOF+j]);
    node_x[j].fastAccessDx(j) = j_coeff;
    if (node_xdot != NULL) {
      (*node_xdot)[j] = FadType(neqn, (*xdot)[firstDOF+j]);
      (*node_xdot)[j].fastAccessDx(j) = m_coeff;
    }
  }
}

void
FEApp::JacobianOp::nodePost(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector< FadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int num_entries;
  double* row_view = 0;
  int row, col;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
        f->ReplaceGlobalValue(row, 0, node_f[offsets[eq_row]].val());
      else if (bc.isShared())
        f->ReplaceGlobalValue(row, 0, 0.0);
    }

    // Always zero out row (This takes care of the not-owned case)
    if (bc.isOwned() || bc.isShared()) {
      jac->ExtractGlobalRowView(row, num_entries, row_view);
      for (int k=0; k<num_entries; k++)
        row_view[k] = 0.0;
    }
	
    // Check derivative array is nonzero
    if (node_f[offsets[eq_row]].hasFastAccess()) {
	    
      // Loop over equations per node
      for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	      
        // Global column
        col = static_cast<int>(firstDOF + eq_col);

        // Replace Jacobian
        if (bc.isOwned())
          jac->ReplaceGlobalValues(row, 1, 
                                   &(node_f[eq_row].fastAccessDx(eq_col)),
                                   &col);

      } // column equations
	
    } // has fast access
    
  } // row equations

}

FEApp::TangentOp::TangentOp(
        double alpha, double beta, bool sum_derivs_,
        const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
        const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
        const Teuchos::RCP<ParamVec>& p,
        const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vx,
        const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vxdot,
        const Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> >& Vp_,
        const Teuchos::RCP<Epetra_Vector>& overlapped_f,
        const Teuchos::RCP<Epetra_MultiVector>& overlapped_JV,
        const Teuchos::RCP<Epetra_MultiVector>& overlapped_fp) :
  m_coeff(alpha),
  j_coeff(beta),
  sum_derivs(sum_derivs_),
  xdot(overlapped_xdot),
  x(overlapped_x),
  params(p),
  Vx(overlapped_Vx),
  Vxdot(overlapped_Vxdot),
  Vp(Vp_),
  f(overlapped_f),
  JV(overlapped_JV),
  fp(overlapped_fp),
  num_cols_x(0),
  num_cols_p(0),
  num_cols_tot(0),
  param_offset(0)
{
  if (Vx != Teuchos::null)
    num_cols_x = Vx->NumVectors();
  else if (Vxdot != Teuchos::null)
    num_cols_x = Vxdot->NumVectors();
  if (params != Teuchos::null) {
    if (Vp != Teuchos::null)
      num_cols_p = Vp->numCols();
    else
      num_cols_p = params->size();
  }
  if (sum_derivs) {
    num_cols_tot = num_cols_x;
    param_offset = 0;
  }
  else {
    num_cols_tot = num_cols_x + num_cols_p;
    param_offset = num_cols_x;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(sum_derivs && (num_cols_x != 0) && (num_cols_p != 0) && 
                     (num_cols_x != num_cols_p),
                     std::logic_error,
                     "Seed matrices Vx and Vp must have the same number " << 
                     " of columns when sum_derivs is true and both are "
                     << "non-null!" << std::endl);
}

FEApp::TangentOp::~TangentOp()
{
}

void
FEApp::TangentOp::elementInit(
			 const FEApp::AbstractElement& e,
			 unsigned int neqn,
			 std::vector< FadType >* elem_xdot,
			 std::vector< FadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = x->Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = x->Map().LID(static_cast<int>(node_GID*neqn));
#endif

    for (unsigned int j=0; j<neqn; j++) {
      if (Vx != Teuchos::null && j_coeff != 0.0) {
        elem_x[neqn*i+j] = 
          FadType(num_cols_tot, (*x)[firstDOF+j]);
        for (int k=0; k<num_cols_x; k++)
          elem_x[neqn*i+j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
      }
      else
        elem_x[neqn*i+j] = 
          FadType((*x)[firstDOF+j]);
      
      if (elem_xdot != NULL) {
        if (Vxdot != Teuchos::null && m_coeff != 0.0) {
          (*elem_xdot)[neqn*i+j] = 
            FadType(num_cols_tot, (*xdot)[firstDOF+j]);
          for (int k=0; k<num_cols_x; k++)
            (*elem_xdot)[neqn*i+j].fastAccessDx(k) = 
              m_coeff*(*Vxdot)[k][firstDOF+j];
        }
        else
          (*elem_xdot)[neqn*i+j] = 
            FadType((*xdot)[firstDOF+j]);
      }   
    }
  }

  if (params != Teuchos::null) {
    FadType p;
    for (unsigned int i=0; i<params->size(); i++) {
      p = FadType(num_cols_tot, (*params)[i].baseValue);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::TangentType>(p);
    }
  }
}

void
FEApp::TangentOp::elementPost(const FEApp::AbstractElement& e,
                              unsigned int neqn,
                              std::vector< FadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
        f->SumIntoGlobalValue(row, 0, elem_f[lrow].val());
	
      // Extract derivative components for each column
      if (JV != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
          JV->SumIntoGlobalValue(row, col, elem_f[lrow].dx(col));
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
          fp->SumIntoGlobalValue(row, col, elem_f[lrow].dx(col+param_offset));
      
    } // row equations
    
  } // row node

}

void
FEApp::TangentOp::nodeInit(const FEApp::NodeBC& bc,
                           unsigned int neqn,
                           std::vector< FadType >* node_xdot,
                           std::vector< FadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = x->Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = x->Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    if (Vx != Teuchos::null && j_coeff != 0.0) {
      node_x[j] = FadType(num_cols_tot, (*x)[firstDOF+j]);
      for (int k=0; k < num_cols_x; k++)
        node_x[j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
    }
    else
      node_x[j] = FadType((*x)[firstDOF+j]);

    if (node_xdot != NULL) {
      if (Vxdot != Teuchos::null && m_coeff != 0.0) {
        (*node_xdot)[j] = 
          FadType(num_cols_tot, (*xdot)[firstDOF+j]);
        for (int k=0; k < num_cols_x; k++)
          (*node_xdot)[j].fastAccessDx(k) = m_coeff*(*Vxdot)[k][firstDOF+j];
      }
      else
        (*node_xdot)[j] = 
          FadType((*xdot)[firstDOF+j]);
    }
  }
  
  if (params != Teuchos::null) {
    FadType p;
    for (unsigned int i=0; i<params->size(); i++) {
      p = FadType(num_cols_tot, (*params)[i].baseValue);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::TangentType>(p);
    }
  }
}

void
FEApp::TangentOp::nodePost(const FEApp::NodeBC& bc,
			    unsigned int neqn,
			    std::vector< FadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int row;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
        f->ReplaceGlobalValue(row, 0, node_f[offsets[eq_row]].val());
      else if (bc.isShared())
        f->ReplaceGlobalValue(row, 0, 0.0);
    }
    
    if (bc.isOwned()) {
      
      // Extract derivative components for each column
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
          JV->ReplaceGlobalValue(row, col, node_f[offsets[eq_row]].dx(col));
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
          fp->ReplaceGlobalValue(row, col, 
                                 node_f[offsets[eq_row]].dx(col+param_offset));
    }
    else if (bc.isShared()) {
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
          JV->ReplaceGlobalValue(row, col, 0.0);
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
          fp->ReplaceGlobalValue(row, col, 0.0);
    }
    
  } // row equations

}

FEApp::SGResidualOp::SGResidualOp(
   const Teuchos::RCP< Stokhos::OrthogPolyExpansion<int,double> >& expansion_,
   const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly >& xdot_,
   const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly >& x_,
   const Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly >& f_) :
  expansion(expansion_),
  nblock(x_->size()),
  xdot(xdot_),
  x(x_),
  f(f_)
{
}

FEApp::SGResidualOp::~SGResidualOp()
{
}

void
FEApp::SGResidualOp::elementInit(const FEApp::AbstractElement& e,
				 unsigned int neqn,
				 std::vector<SGType>* elem_xdot,
				 std::vector<SGType>& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Allocate appropriate sizes for coefficients
  for (unsigned int i=0; i<nnode; i++) {
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j].reset(expansion);
      elem_x[neqn*i+j].copyForWrite();
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j].reset(expansion);
        (*elem_xdot)[neqn*i+j].copyForWrite();
      }
    }
  }

  for (int block=0; block<nblock; block++) {

    // Copy element solution
    for (unsigned int i=0; i<nnode; i++) {
      node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
      firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif
      for (unsigned int j=0; j<neqn; j++) {
        elem_x[neqn*i+j].fastAccessCoeff(block) = (*x)[block][firstDOF+j];
        if (elem_xdot != NULL)
          (*elem_xdot)[neqn*i+j].fastAccessCoeff(block) = 
            (*xdot)[block][firstDOF+j];
      }
    }
    
  }
}

void
FEApp::SGResidualOp::elementPost(const FEApp::AbstractElement& e,
                                 unsigned int neqn,
                                 std::vector<SGType>& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  int row;
  unsigned int lrow;

  for (int block=0; block<nblock; block++) {

    // Sum element residual into global residual
    for (unsigned int node_row=0; node_row<nnode; node_row++) {
      for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
        lrow = neqn*node_row+eq_row;
        row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
        (*f)[block].SumIntoGlobalValue(row, 0, elem_f[lrow].coeff(block));
      }
    }
    
  }
}

void
FEApp::SGResidualOp::nodeInit(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector<SGType>* node_xdot,
                              std::vector<SGType>& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Allocate appropriate sizes for coefficients
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j].reset(expansion);
    node_x[j].copyForWrite();
    if (node_xdot != NULL) {
      (*node_xdot)[j].reset(expansion);
      (*node_xdot)[j].copyForWrite();
    }
  }

  for (int block=0; block<nblock; block++) {

    // Copy node solution
    for (unsigned int j=0; j<neqn; j++) {
      node_x[j].fastAccessCoeff(block) = (*x)[block][firstDOF+j];
      if (node_xdot != NULL)
        (*node_xdot)[j].fastAccessCoeff(block) = (*xdot)[block][firstDOF+j];
    }

  }
}

void
FEApp::SGResidualOp::nodePost(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector<SGType>& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  for (int block=0; block<nblock; block++) {

    // Replace global residual
    for (unsigned int j=0; j<offsets.size(); j++) {
      int row = static_cast<int>(firstDOF + offsets[j]);
      if (bc.isOwned())
        (*f)[block].ReplaceGlobalValue(row, 0, 
				       node_f[offsets[j]].coeff(block));
      else if (bc.isShared())
        (*f)[block].ReplaceGlobalValue(row, 0, 0.0);
    }

  }
}

void
FEApp::SGResidualOp::finalizeFill()
{
}

FEApp::SGJacobianOp::SGJacobianOp(
   const Teuchos::RCP< Stokhos::OrthogPolyExpansion<int,double> >& expansion_,
   double alpha, double beta,
   const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly >& xdot_,
   const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly >& x_,
   const Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly >& f_,
   const Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_CrsMatrix> >& jac_) :
  expansion(expansion_),
  nblock(x_->size()),
  m_coeff(alpha),
  j_coeff(beta),
  xdot(xdot_),
  x(x_),
  f(f_),
  jac(jac_)
{
}

FEApp::SGJacobianOp::~SGJacobianOp()
{
}

void
FEApp::SGJacobianOp::elementInit(const FEApp::AbstractElement& e,
				 unsigned int neqn,
				 std::vector< SGFadType >* elem_xdot,
				 std::vector< SGFadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = SGFadType(ndof, 0.0);
      elem_x[neqn*i+j].fastAccessDx(neqn*i+j) = j_coeff;
      elem_x[neqn*i+j].val().reset(expansion);
      elem_x[neqn*i+j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
        elem_x[neqn*i+j].val().fastAccessCoeff(block) = 
          (*x)[block][firstDOF+j];
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j] = SGFadType(ndof, 0.0);
        (*elem_xdot)[neqn*i+j].fastAccessDx(neqn*i+j) = m_coeff;
        (*elem_xdot)[neqn*i+j].val().reset(expansion);
        (*elem_xdot)[neqn*i+j].val().copyForWrite();
        for (int block=0; block<nblock; block++)
          (*elem_xdot)[neqn*i+j].val().fastAccessCoeff(block) = 
            (*xdot)[block][firstDOF+j];
      }
      
    }
  }

}

void
FEApp::SGJacobianOp::elementPost(const FEApp::AbstractElement& e,
                                 unsigned int neqn,
                                 std::vector< SGFadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  int row, col;
  unsigned int lrow, lcol;
  double c;
  int nblock_jac = jac->size();
  
  // Loop over nodes in element
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
	for (int block=0; block<nblock; block++)
          (*f)[block].SumIntoGlobalValue(row, 0, 
					 elem_f[lrow].val().coeff(block));
      
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {
        
	// Loop over nodes in element
	for (unsigned int node_col=0; node_col<nnode; node_col++){
          
	  // Loop over equations per node
	  for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	    lcol = neqn*node_col+eq_col;
            
	    // Global column
	    col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);

	    TEUCHOS_TEST_FOR_EXCEPTION(elem_f[lrow].fastAccessDx(lcol).size() >
			       jac->size(), 
			       std::logic_error,
			       "Jacobian entry polynomial has size " <<
			       elem_f[lrow].fastAccessDx(lcol).size() <<
			       ", but Jacobian polynomial only has size " <<
			       jac->size() << "!")
            
	    // Sum Jacobian
	    for (int block=0; block<nblock_jac; block++) {
	      c = elem_f[lrow].fastAccessDx(lcol).coeff(block);
	      (*jac)[block].SumIntoGlobalValues(row, 1, &c, &col);
	    }
            
	  } // column equations
          
	} // column nodes
        
      } // has fast access
      
    } // row equations
    
  } // row node

}

void
FEApp::SGJacobianOp::nodeInit(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector< SGFadType >* node_xdot,
                              std::vector< SGFadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = SGFadType(neqn, 0.0);
    node_x[j].fastAccessDx(j) = j_coeff;
    node_x[j].val().reset(expansion);
    node_x[j].val().copyForWrite();
    for (int block=0; block<nblock; block++)
      node_x[j].val().fastAccessCoeff(block) = (*x)[block][firstDOF+j];
    if (node_xdot != NULL) {
      (*node_xdot)[j] = SGFadType(neqn, 0.0);
      (*node_xdot)[j].fastAccessDx(j) = m_coeff;
      (*node_xdot)[j].val().reset(expansion);
      (*node_xdot)[j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
        (*node_xdot)[j].val().fastAccessCoeff(block) = 
          (*xdot)[block][firstDOF+j];
    }
  }
}

void
FEApp::SGJacobianOp::nodePost(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector< SGFadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int num_entries;
  double* row_view = 0;
  int row, col;
  double c;
  int nblock_jac = jac->size();

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
	for (int block=0; block<nblock; block++)
          (*f)[block].ReplaceGlobalValue(row, 0, 
					 node_f[offsets[eq_row]].val().coeff(block));
      else if (bc.isShared())
	for (int block=0; block<nblock; block++)
          (*f)[block].ReplaceGlobalValue(row, 0, 0.0);
    }
    
    // Always zero out row (This takes care of the not-owned case)
    if (bc.isOwned() || bc.isShared()) {
      for (int block=0; block<nblock_jac; block++) {
	(*jac)[block].ExtractGlobalRowView(row, num_entries, row_view);
	for (int k=0; k<num_entries; k++)
	  row_view[k] = 0.0;
      }
    }
    
    // Check derivative array is nonzero
    if (node_f[offsets[eq_row]].hasFastAccess()) {
	    
      // Loop over equations per node
      for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	      
	// Global column
	col = static_cast<int>(firstDOF + eq_col);
        
	// Replace Jacobian
	if (bc.isOwned()) {
	  for (int block=0; block<nblock_jac; block++) {
            c = node_f[eq_row].fastAccessDx(eq_col).coeff(block);
	    (*jac)[block].ReplaceGlobalValues(row, 1, &c, &col);
	  }
	}
        
      } // column equations
      
    } // has fast access
    
  } // row equations

}

void
FEApp::SGJacobianOp::finalizeFill()
{
}

FEApp::SGTangentOp::SGTangentOp(
  const Teuchos::RCP< Stokhos::OrthogPolyExpansion<int,double> >& expansion_,
  double alpha, double beta, bool sum_derivs_,
  const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly >& xdot_,
  const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly >& x_,
  const Teuchos::RCP<ParamVec>& p,
  const Teuchos::RCP<const Epetra_MultiVector>& Vx_,
  const Teuchos::RCP<const Epetra_MultiVector>& Vxdot_,
  const Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> >& Vp_,
  const Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly >& f_,
  const Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly >& JV_,
  const Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly >& fp_) :
  expansion(expansion_),
  nblock(x_->size()),
  m_coeff(alpha),
  j_coeff(beta),
  sum_derivs(sum_derivs_),
  xdot(xdot_),
  x(x_),
  params(p),
  Vx(Vx_),
  Vxdot(Vxdot_),
  Vp(Vp_),
  f(f_),
  JV(JV_),
  fp(fp_),
  num_cols_x(0),
  num_cols_p(0),
  num_cols_tot(0),
  param_offset(0)
{
  if (Vx != Teuchos::null)
    num_cols_x = Vx->NumVectors();
  else if (Vxdot != Teuchos::null)
    num_cols_x = Vxdot->NumVectors();
  if (params != Teuchos::null) {
    if (Vp != Teuchos::null)
      num_cols_p = Vp->numCols();
    else
      num_cols_p = params->size();
  }
  if (sum_derivs) {
    num_cols_tot = num_cols_x;
    param_offset = 0;
  }
  else {
    num_cols_tot = num_cols_x + num_cols_p;
    param_offset = num_cols_x;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(sum_derivs && (num_cols_x != 0) && (num_cols_p != 0) && 
                     (num_cols_x != num_cols_p),
                     std::logic_error,
                     "Seed matrices Vx and Vp must have the same number " << 
                     " of columns when sum_derivs is true and both are "
                     << "non-null!" << std::endl);
}

FEApp::SGTangentOp::~SGTangentOp()
{
}

void
FEApp::SGTangentOp::elementInit(
			 const FEApp::AbstractElement& e,
			 unsigned int neqn,
			 std::vector< SGFadType >* elem_xdot,
			 std::vector< SGFadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

    for (unsigned int j=0; j<neqn; j++) {
      if (Vx != Teuchos::null && j_coeff != 0.0) {
        elem_x[neqn*i+j] = SGFadType(num_cols_tot, 0.0);
        for (int k=0; k<num_cols_x; k++)
          elem_x[neqn*i+j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
      }
      else
        elem_x[neqn*i+j] = SGFadType(0.0);
      elem_x[neqn*i+j].val().reset(expansion);
      elem_x[neqn*i+j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
	elem_x[neqn*i+j].val().fastAccessCoeff(block) = 
	  (*x)[block][firstDOF+j];
      
      if (elem_xdot != NULL) {
        if (Vxdot != Teuchos::null && m_coeff != 0.0) {
          (*elem_xdot)[neqn*i+j] = SGFadType(num_cols_tot, 0.0);
          for (int k=0; k<num_cols_x; k++)
            (*elem_xdot)[neqn*i+j].fastAccessDx(k) = 
              m_coeff*(*Vxdot)[k][firstDOF+j];
        }
        else
          (*elem_xdot)[neqn*i+j] = SGFadType(0.0);
	(*elem_xdot)[neqn*i+j].val().reset(expansion);
	(*elem_xdot)[neqn*i+j].val().copyForWrite();
	for (int block=0; block<nblock; block++)
	  (*elem_xdot)[neqn*i+j].val().fastAccessCoeff(block) = 
	    (*xdot)[block][firstDOF+j];
      }   
    }
  }

  if (params != Teuchos::null) {
    SGFadType p;
    SGType pval;
    for (unsigned int i=0; i<params->size(); i++) {
      pval = (*params)[i].family->getValue<FEApp::SGTangentType>().val();
      p = SGFadType(num_cols_tot, pval);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::SGTangentType>(p);
    }
  }
}

void
FEApp::SGTangentOp::elementPost(const FEApp::AbstractElement& e,
                              unsigned int neqn,
                              std::vector< SGFadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
	for (int block=0; block<nblock; block++)
	  (*f)[block].SumIntoGlobalValue(row, 0, 
					 elem_f[lrow].val().coeff(block));
	
      // Extract derivative components for each column
      if (JV != Teuchos::null)
	for (int col=0; col<num_cols_x; col++)
	  for (int block=0; block<nblock; block++)
	    (*JV)[block].SumIntoGlobalValue(row, col, 
					    elem_f[lrow].dx(col).coeff(block));
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
	  for (int block=0; block<nblock; block++)
	    (*fp)[block].SumIntoGlobalValue(row, col, 
					    elem_f[lrow].dx(col+param_offset).coeff(block));
      
    } // row equations
    
  } // row node

}

void
FEApp::SGTangentOp::nodeInit(const FEApp::NodeBC& bc,
                           unsigned int neqn,
                           std::vector< SGFadType >* node_xdot,
                           std::vector< SGFadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    if (Vx != Teuchos::null && j_coeff != 0.0) {
      node_x[j] = SGFadType(num_cols_tot, 0.0);
      for (int k=0; k < num_cols_x; k++)
        node_x[j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
    }
    else
      node_x[j] = SGFadType(0.0);
    node_x[j].val().reset(expansion);
    node_x[j].val().copyForWrite();
    for (int block=0; block<nblock; block++)
      node_x[j].val().fastAccessCoeff(block) = (*x)[block][firstDOF+j];

    if (node_xdot != NULL) {
      if (Vxdot != Teuchos::null && m_coeff != 0.0) {
        (*node_xdot)[j] = SGFadType(num_cols_tot, 0.0);
        for (int k=0; k < num_cols_x; k++)
          (*node_xdot)[j].fastAccessDx(k) = m_coeff*(*Vxdot)[k][firstDOF+j];
      }
      else
        (*node_xdot)[j] = SGFadType(0.0);
      (*node_xdot)[j].val().reset(expansion);
      (*node_xdot)[j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
        (*node_xdot)[j].val().fastAccessCoeff(block) = 
          (*xdot)[block][firstDOF+j];
    }
  }
  
  if (params != Teuchos::null) {
    SGFadType p;
    SGType pval;
    for (unsigned int i=0; i<params->size(); i++) {
      pval = (*params)[i].family->getValue<FEApp::SGTangentType>().val();
      p = SGFadType(num_cols_tot, pval);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::SGTangentType>(p);
    }
  }
}

void
FEApp::SGTangentOp::nodePost(const FEApp::NodeBC& bc,
			    unsigned int neqn,
			    std::vector< SGFadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int row;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
	for (int block=0; block<nblock; block++)
	  (*f)[block].ReplaceGlobalValue(row, 0, 
					 node_f[offsets[eq_row]].val().coeff(block));
      else if (bc.isShared())
	for (int block=0; block<nblock; block++)
	  (*f)[block].ReplaceGlobalValue(row, 0, 0.0);
    }
    
    if (bc.isOwned()) {
      
      // Extract derivative components for each column
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
	  for (int block=0; block<nblock; block++)
	    (*JV)[block].ReplaceGlobalValue(row, col, 
					    node_f[offsets[eq_row]].dx(col).coeff(block));
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
	  for (int block=0; block<nblock; block++)
	    (*fp)[block].ReplaceGlobalValue(row, col, 
					    node_f[offsets[eq_row]].dx(col+param_offset).coeff(block));
    }
    else if (bc.isShared()) {
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
	  for (int block=0; block<nblock; block++)
	    (*JV)[block].ReplaceGlobalValue(row, col, 0.0);
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
	  for (int block=0; block<nblock; block++)
	    (*fp)[block].ReplaceGlobalValue(row, col, 0.0);
    }
    
  } // row equations

}

FEApp::MPResidualOp::MPResidualOp(
   const Teuchos::RCP<const Stokhos::ProductEpetraVector >& xdot_,
   const Teuchos::RCP<const Stokhos::ProductEpetraVector >& x_,
   const Teuchos::RCP< Stokhos::ProductEpetraVector >& f_) :
  nblock(x_->size()),
  xdot(xdot_),
  x(x_),
  f(f_)
{
}

FEApp::MPResidualOp::~MPResidualOp()
{
}

void
FEApp::MPResidualOp::elementInit(const FEApp::AbstractElement& e,
				 unsigned int neqn,
				 std::vector<MPType>* elem_xdot,
				 std::vector<MPType>& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Allocate appropriate sizes for coefficients
  for (unsigned int i=0; i<nnode; i++) {
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j].reset(nblock);
      elem_x[neqn*i+j].copyForWrite();
      if (elem_xdot != NULL) {
	(*elem_xdot)[neqn*i+j].reset(nblock);
        (*elem_xdot)[neqn*i+j].copyForWrite();
      }
    }
  }

  for (int block=0; block<nblock; block++) {

    // Copy element solution
    for (unsigned int i=0; i<nnode; i++) {
      node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
      firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif
      for (unsigned int j=0; j<neqn; j++) {
        elem_x[neqn*i+j].fastAccessCoeff(block) = (*x)[block][firstDOF+j];
        if (elem_xdot != NULL)
          (*elem_xdot)[neqn*i+j].fastAccessCoeff(block) = 
            (*xdot)[block][firstDOF+j];
      }
    }
    
  }
}

void
FEApp::MPResidualOp::elementPost(const FEApp::AbstractElement& e,
                                 unsigned int neqn,
                                 std::vector<MPType>& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  int row;
  unsigned int lrow;

  for (int block=0; block<nblock; block++) {

    // Sum element residual into global residual
    for (unsigned int node_row=0; node_row<nnode; node_row++) {
      for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
        lrow = neqn*node_row+eq_row;
        row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
        (*f)[block].SumIntoGlobalValue(row, 0, elem_f[lrow].coeff(block));
      }
    }
    
  }
}

void
FEApp::MPResidualOp::nodeInit(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector<MPType>* node_xdot,
                              std::vector<MPType>& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Allocate appropriate sizes for coefficients
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j].reset(nblock);
    node_x[j].copyForWrite();
    if (node_xdot != NULL) {
      (*node_xdot)[j].reset(nblock);
      (*node_xdot)[j].copyForWrite();
    }
  }

  for (int block=0; block<nblock; block++) {

    // Copy node solution
    for (unsigned int j=0; j<neqn; j++) {
      node_x[j].fastAccessCoeff(block) = (*x)[block][firstDOF+j];
      if (node_xdot != NULL)
        (*node_xdot)[j].fastAccessCoeff(block) = (*xdot)[block][firstDOF+j];
    }

  }
}

void
FEApp::MPResidualOp::nodePost(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector<MPType>& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  for (int block=0; block<nblock; block++) {

    // Replace global residual
    for (unsigned int j=0; j<offsets.size(); j++) {
      int row = static_cast<int>(firstDOF + offsets[j]);
      if (bc.isOwned())
        (*f)[block].ReplaceGlobalValue(row, 0, 
				       node_f[offsets[j]].coeff(block));
      else if (bc.isShared())
        (*f)[block].ReplaceGlobalValue(row, 0, 0.0);
    }

  }
}

void
FEApp::MPResidualOp::finalizeFill()
{
}

FEApp::MPJacobianOp::MPJacobianOp(
   double alpha, double beta,
   const Teuchos::RCP<const Stokhos::ProductEpetraVector >& xdot_,
   const Teuchos::RCP<const Stokhos::ProductEpetraVector >& x_,
   const Teuchos::RCP< Stokhos::ProductEpetraVector >& f_,
   const Teuchos::RCP< Stokhos::ProductContainer<Epetra_CrsMatrix> >& jac_) :
  nblock(x_->size()),
  m_coeff(alpha),
  j_coeff(beta),
  xdot(xdot_),
  x(x_),
  f(f_),
  jac(jac_)
{
}

FEApp::MPJacobianOp::~MPJacobianOp()
{
}

void
FEApp::MPJacobianOp::elementInit(const FEApp::AbstractElement& e,
				 unsigned int neqn,
				 std::vector< MPFadType >* elem_xdot,
				 std::vector< MPFadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = MPFadType(ndof, 0.0);
      elem_x[neqn*i+j].fastAccessDx(neqn*i+j) = j_coeff;
      elem_x[neqn*i+j].val().reset(nblock);
      elem_x[neqn*i+j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
        elem_x[neqn*i+j].val().fastAccessCoeff(block) = 
          (*x)[block][firstDOF+j];
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j] = MPFadType(ndof, 0.0);
        (*elem_xdot)[neqn*i+j].fastAccessDx(neqn*i+j) = m_coeff;
	(*elem_xdot)[neqn*i+j].val().reset(nblock);
        (*elem_xdot)[neqn*i+j].val().copyForWrite();
        for (int block=0; block<nblock; block++)
          (*elem_xdot)[neqn*i+j].val().fastAccessCoeff(block) = 
            (*xdot)[block][firstDOF+j];
      }
      
    }
  }

}

void
FEApp::MPJacobianOp::elementPost(const FEApp::AbstractElement& e,
                                 unsigned int neqn,
                                 std::vector< MPFadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  int row, col;
  unsigned int lrow, lcol;
  double c;
  int nblock_jac = jac->size();
  
  // Loop over nodes in element
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
	for (int block=0; block<nblock; block++)
          (*f)[block].SumIntoGlobalValue(row, 0, 
					 elem_f[lrow].val().coeff(block));
      
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {
        
	// Loop over nodes in element
	for (unsigned int node_col=0; node_col<nnode; node_col++){
          
	  // Loop over equations per node
	  for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	    lcol = neqn*node_col+eq_col;
            
	    // Global column
	    col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);

	    TEUCHOS_TEST_FOR_EXCEPTION(elem_f[lrow].fastAccessDx(lcol).size() >
			       jac->size(), 
			       std::logic_error,
			       "Jacobian entry vector has size " <<
			       elem_f[lrow].fastAccessDx(lcol).size() <<
			       ", but Jacobian vector only has size " <<
			       jac->size() << "!")
            
	    // Sum Jacobian
	    for (int block=0; block<nblock_jac; block++) {
	      c = elem_f[lrow].fastAccessDx(lcol).coeff(block);
	      (*jac)[block].SumIntoGlobalValues(row, 1, &c, &col);
	    }
            
	  } // column equations
          
	} // column nodes
        
      } // has fast access
      
    } // row equations
    
  } // row node

}

void
FEApp::MPJacobianOp::nodeInit(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector< MPFadType >* node_xdot,
                              std::vector< MPFadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = MPFadType(neqn, 0.0);
    node_x[j].fastAccessDx(j) = j_coeff;
    node_x[j].val().reset(nblock);
    node_x[j].val().copyForWrite();
    for (int block=0; block<nblock; block++)
      node_x[j].val().fastAccessCoeff(block) = (*x)[block][firstDOF+j];
    if (node_xdot != NULL) {
      (*node_xdot)[j] = MPFadType(neqn, 0.0);
      (*node_xdot)[j].fastAccessDx(j) = m_coeff;
      (*node_xdot)[j].val().reset(nblock);
      (*node_xdot)[j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
        (*node_xdot)[j].val().fastAccessCoeff(block) = 
          (*xdot)[block][firstDOF+j];
    }
  }
}

void
FEApp::MPJacobianOp::nodePost(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector< MPFadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int num_entries;
  double* row_view = 0;
  int row, col;
  double c;
  int nblock_jac = jac->size();

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
	for (int block=0; block<nblock; block++)
          (*f)[block].ReplaceGlobalValue(row, 0, 
					 node_f[offsets[eq_row]].val().coeff(block));
      else if (bc.isShared())
	for (int block=0; block<nblock; block++)
          (*f)[block].ReplaceGlobalValue(row, 0, 0.0);
    }
    
    // Always zero out row (This takes care of the not-owned case)
    if (bc.isOwned() || bc.isShared()) {
      for (int block=0; block<nblock_jac; block++) {
	(*jac)[block].ExtractGlobalRowView(row, num_entries, row_view);
	for (int k=0; k<num_entries; k++)
	  row_view[k] = 0.0;
      }
    }
    
    // Check derivative array is nonzero
    if (node_f[offsets[eq_row]].hasFastAccess()) {
	    
      // Loop over equations per node
      for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	      
	// Global column
	col = static_cast<int>(firstDOF + eq_col);
        
	// Replace Jacobian
	if (bc.isOwned()) {
	  for (int block=0; block<nblock_jac; block++) {
            c = node_f[eq_row].fastAccessDx(eq_col).coeff(block);
	    (*jac)[block].ReplaceGlobalValues(row, 1, &c, &col);
	  }
	}
        
      } // column equations
      
    } // has fast access
    
  } // row equations

}

void
FEApp::MPJacobianOp::finalizeFill()
{
}

FEApp::MPTangentOp::MPTangentOp(
  double alpha, double beta, bool sum_derivs_,
  const Teuchos::RCP<const Stokhos::ProductEpetraVector >& xdot_,
  const Teuchos::RCP<const Stokhos::ProductEpetraVector >& x_,
  const Teuchos::RCP<ParamVec>& p,
  const Teuchos::RCP<const Epetra_MultiVector>& Vx_,
  const Teuchos::RCP<const Epetra_MultiVector>& Vxdot_,
  const Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> >& Vp_,
  const Teuchos::RCP< Stokhos::ProductEpetraVector >& f_,
  const Teuchos::RCP< Stokhos::ProductEpetraMultiVector >& JV_,
  const Teuchos::RCP< Stokhos::ProductEpetraMultiVector >& fp_) :
  nblock(x_->size()),
  m_coeff(alpha),
  j_coeff(beta),
  sum_derivs(sum_derivs_),
  xdot(xdot_),
  x(x_),
  params(p),
  Vx(Vx_),
  Vxdot(Vxdot_),
  Vp(Vp_),
  f(f_),
  JV(JV_),
  fp(fp_),
  num_cols_x(0),
  num_cols_p(0),
  num_cols_tot(0),
  param_offset(0)
{
  if (Vx != Teuchos::null)
    num_cols_x = Vx->NumVectors();
  else if (Vxdot != Teuchos::null)
    num_cols_x = Vxdot->NumVectors();
  if (params != Teuchos::null) {
    if (Vp != Teuchos::null)
      num_cols_p = Vp->numCols();
    else
      num_cols_p = params->size();
  }
  if (sum_derivs) {
    num_cols_tot = num_cols_x;
    param_offset = 0;
  }
  else {
    num_cols_tot = num_cols_x + num_cols_p;
    param_offset = num_cols_x;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(sum_derivs && (num_cols_x != 0) && (num_cols_p != 0) && 
                     (num_cols_x != num_cols_p),
                     std::logic_error,
                     "Seed matrices Vx and Vp must have the same number " << 
                     " of columns when sum_derivs is true and both are "
                     << "non-null!" << std::endl);
}

FEApp::MPTangentOp::~MPTangentOp()
{
}

void
FEApp::MPTangentOp::elementInit(
			 const FEApp::AbstractElement& e,
			 unsigned int neqn,
			 std::vector< MPFadType >* elem_xdot,
			 std::vector< MPFadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
    firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

    for (unsigned int j=0; j<neqn; j++) {
      if (Vx != Teuchos::null && j_coeff != 0.0) {
        elem_x[neqn*i+j] = MPFadType(num_cols_tot, 0.0);
        for (int k=0; k<num_cols_x; k++)
          elem_x[neqn*i+j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
      }
      else
        elem_x[neqn*i+j] = MPFadType(0.0);
      elem_x[neqn*i+j].val().reset(nblock);
      elem_x[neqn*i+j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
	elem_x[neqn*i+j].val().fastAccessCoeff(block) = 
	  (*x)[block][firstDOF+j];
      
      if (elem_xdot != NULL) {
        if (Vxdot != Teuchos::null && m_coeff != 0.0) {
          (*elem_xdot)[neqn*i+j] = MPFadType(num_cols_tot, 0.0);
          for (int k=0; k<num_cols_x; k++)
            (*elem_xdot)[neqn*i+j].fastAccessDx(k) = 
              m_coeff*(*Vxdot)[k][firstDOF+j];
        }
        else
          (*elem_xdot)[neqn*i+j] = MPFadType(0.0);
	(*elem_xdot)[neqn*i+j].val().reset(nblock);
	(*elem_xdot)[neqn*i+j].val().copyForWrite();
	for (int block=0; block<nblock; block++)
	  (*elem_xdot)[neqn*i+j].val().fastAccessCoeff(block) = 
	    (*xdot)[block][firstDOF+j];
      }   
    }
  }

  if (params != Teuchos::null) {
    MPFadType p;
    MPType pval;
    for (unsigned int i=0; i<params->size(); i++) {
      pval = (*params)[i].family->getValue<FEApp::MPTangentType>().val();
      p = MPFadType(num_cols_tot, pval);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::MPTangentType>(p);
    }
  }
}

void
FEApp::MPTangentOp::elementPost(const FEApp::AbstractElement& e,
                              unsigned int neqn,
                              std::vector< MPFadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
	for (int block=0; block<nblock; block++)
	  (*f)[block].SumIntoGlobalValue(row, 0, 
					 elem_f[lrow].val().coeff(block));
	
      // Extract derivative components for each column
      if (JV != Teuchos::null)
	for (int col=0; col<num_cols_x; col++)
	  for (int block=0; block<nblock; block++)
	    (*JV)[block].SumIntoGlobalValue(row, col, 
					    elem_f[lrow].dx(col).coeff(block));
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
	  for (int block=0; block<nblock; block++)
	    (*fp)[block].SumIntoGlobalValue(row, col, 
					    elem_f[lrow].dx(col+param_offset).coeff(block));
      
    } // row equations
    
  } // row node

}

void
FEApp::MPTangentOp::nodeInit(const FEApp::NodeBC& bc,
                           unsigned int neqn,
                           std::vector< MPFadType >* node_xdot,
                           std::vector< MPFadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<long long>(node_GID*neqn));
#else
  unsigned int firstDOF = (*x)[0].Map().LID(static_cast<int>(node_GID*neqn));
#endif

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    if (Vx != Teuchos::null && j_coeff != 0.0) {
      node_x[j] = MPFadType(num_cols_tot, 0.0);
      for (int k=0; k < num_cols_x; k++)
        node_x[j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
    }
    else
      node_x[j] = MPFadType(0.0);
    node_x[j].val().reset(nblock);
    node_x[j].val().copyForWrite();
    for (int block=0; block<nblock; block++)
      node_x[j].val().fastAccessCoeff(block) = (*x)[block][firstDOF+j];

    if (node_xdot != NULL) {
      if (Vxdot != Teuchos::null && m_coeff != 0.0) {
        (*node_xdot)[j] = MPFadType(num_cols_tot, 0.0);
        for (int k=0; k < num_cols_x; k++)
          (*node_xdot)[j].fastAccessDx(k) = m_coeff*(*Vxdot)[k][firstDOF+j];
      }
      else
        (*node_xdot)[j] = MPFadType(0.0);
      (*node_xdot)[j].val().reset(nblock);
      (*node_xdot)[j].val().copyForWrite();
      for (int block=0; block<nblock; block++)
        (*node_xdot)[j].val().fastAccessCoeff(block) = 
          (*xdot)[block][firstDOF+j];
    }
  }
  
  if (params != Teuchos::null) {
    MPFadType p;
    MPType pval;
    for (unsigned int i=0; i<params->size(); i++) {
      pval = (*params)[i].family->getValue<FEApp::MPTangentType>().val();
      p = MPFadType(num_cols_tot, pval);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::MPTangentType>(p);
    }
  }
}

void
FEApp::MPTangentOp::nodePost(const FEApp::NodeBC& bc,
			    unsigned int neqn,
			    std::vector< MPFadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int row;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
	for (int block=0; block<nblock; block++)
	  (*f)[block].ReplaceGlobalValue(row, 0, 
					 node_f[offsets[eq_row]].val().coeff(block));
      else if (bc.isShared())
	for (int block=0; block<nblock; block++)
	  (*f)[block].ReplaceGlobalValue(row, 0, 0.0);
    }
    
    if (bc.isOwned()) {
      
      // Extract derivative components for each column
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
	  for (int block=0; block<nblock; block++)
	    (*JV)[block].ReplaceGlobalValue(row, col, 
					    node_f[offsets[eq_row]].dx(col).coeff(block));
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
	  for (int block=0; block<nblock; block++)
	    (*fp)[block].ReplaceGlobalValue(row, col, 
					    node_f[offsets[eq_row]].dx(col+param_offset).coeff(block));
    }
    else if (bc.isShared()) {
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
	  for (int block=0; block<nblock; block++)
	    (*JV)[block].ReplaceGlobalValue(row, col, 0.0);
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
	  for (int block=0; block<nblock; block++)
	    (*fp)[block].ReplaceGlobalValue(row, col, 0.0);
    }
    
  } // row equations

}
