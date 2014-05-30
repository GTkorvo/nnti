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

#ifndef FEAPP_APPLICATION_HPP
#define FEAPP_APPLICATION_HPP

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "FEApp_NodeBC.hpp"
#include "FEApp_AbstractPDE.hpp"
#include "FEApp_AbstractQuadrature.hpp"
#include "FEApp_AbstractDiscretization.hpp"
#include "FEApp_AbstractProblem.hpp"
#include "FEApp_GlobalFill.hpp"
#include "FEApp_TemplateTypes.hpp"
#include "FEApp_InitPostOps.hpp"

#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterVector.hpp"

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "Stokhos_ParallelData.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "EpetraExt_MultiComm.h"

namespace FEApp {

  class Application {
  public:

    //! Constructor 
    Application(
      const Teuchos::RCP<const Epetra_Comm>& comm,
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      const Epetra_Vector* initial_soln = NULL);

    //! Destructor
    ~Application();

    //! Get DOF map
    Teuchos::RCP<const Epetra_Map> getMap() const;

    //! Get response map
    Teuchos::RCP<const Epetra_Map> getResponseMap() const;

    //! Get Jacobian graph
    Teuchos::RCP<const Epetra_CrsGraph> getJacobianGraph() const;

    //! Get initial solution
    Teuchos::RCP<const Epetra_Vector> getInitialSolution() const;

    //! Get parameter library
    Teuchos::RCP<ParamLib> getParamLib();

    //! Initialize SG expansion data
    void init_sg(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sg_quad,
      const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& sg_exp,
      const Teuchos::RCP<const EpetraExt::MultiComm>& multiComm);

    //! Create new W operator
    Teuchos::RCP<Epetra_Operator> createW() const;

    //! Create new preconditioner operator
    Teuchos::RCP<Epetra_Operator> createPrec() const;

    //! Compute global residual
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalResidual(const Epetra_Vector* xdot,
			       const Epetra_Vector& x,
			       const Teuchos::Array<ParamVec>& p,
			       Epetra_Vector& f);

    //! Compute global Jacobian
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalJacobian(double alpha, double beta,
			       const Epetra_Vector* xdot,
			       const Epetra_Vector& x,
			       const Teuchos::Array<ParamVec>& p,
			       Epetra_Vector* f,
			       Epetra_Operator& jac);
    
    //! Compute global Preconditioner
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalPreconditioner(double alpha, double beta,
				     const Epetra_Vector* xdot,
				     const Epetra_Vector& x,
				     const Teuchos::Array<ParamVec>& p,
				     Epetra_Vector* f,
				     Epetra_Operator& jac);
    
    //! Compute global Tangent
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalTangent(double alpha, double beta,
                              bool sum_derivs,
                              const Epetra_Vector* xdot,
                              const Epetra_Vector& x,
			      const Teuchos::Array<ParamVec>& p,
                              const Teuchos::RCP<ParamVec>& deriv_p,
                              const Epetra_MultiVector* Vx,
                              const Teuchos::SerialDenseMatrix<int,double>* Vp,
                              Epetra_Vector* f,
                              Epetra_MultiVector* JVx,
                              Epetra_MultiVector* fVp);

    //! Evaluate response functions
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void 
    evaluateResponses(const Epetra_Vector* xdot,
                      const Epetra_Vector& x,
                      const Teuchos::Array<ParamVec>& p,
                      Epetra_Vector& g);

    //! Evaluate tangent = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
    /*!
     * Set xdot, dxdot_dp to NULL for steady-state problems
     */
    void 
    evaluateResponseTangents(
      const Epetra_Vector* xdot,
      const Epetra_Vector& x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
      Epetra_Vector* g,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& gt);
    
    //! Evaluate gradient = dg/dx, dg/dxdot, dg/dp
    /*!
     * Set xdot, dg_dxdot to NULL for steady-state problems
     */
    void 
    evaluateResponseGradients(
      const Epetra_Vector* xdot,
      const Epetra_Vector& x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      Epetra_Vector* g,
      Epetra_MultiVector* dg_dx,
      Epetra_MultiVector* dg_dxdot,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dg_dp);

    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >
    getStochasticExpansion() { return sg_expansion; }

    //! Compute global residual for stochastic Galerkin problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalSGResidual(
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      Stokhos::EpetraVectorOrthogPoly& sg_f);

    //! Compute global Jacobian for stochastic Galerkin problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalSGJacobian(
      double alpha, double beta,
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      Stokhos::EpetraVectorOrthogPoly* sg_f,
      Stokhos::EpetraOperatorOrthogPoly& sg_jac);

    //! Compute global Tangent for stochastic Galerkin problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalSGTangent(
      double alpha, double beta, bool sum_derivs,
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      const Teuchos::RCP<ParamVec>& deriv_p,
      const Epetra_MultiVector* Vx,
      const Teuchos::SerialDenseMatrix<int,double>* Vp,
      Stokhos::EpetraVectorOrthogPoly* sg_f,
      Stokhos::EpetraMultiVectorOrthogPoly* sg_JVx,
      Stokhos::EpetraMultiVectorOrthogPoly* sg_fVp);

    //! Evaluate stochastic Galerkin response functions
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void 
    evaluateSGResponses(
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      Stokhos::EpetraVectorOrthogPoly& sg_g);

    //! Evaluate tangent = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
    /*!
     * Set xdot, dxdot_dp to NULL for steady-state problems
     */
    void 
    evaluateSGResponseTangents(
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
      Stokhos::EpetraVectorOrthogPoly* sg_g,
      const Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > >& sg_gt);

    //! Evaluate gradient = dg/dx, dg/dxdot, dg/dp
    /*!
     * Set xdot, dg_dxdot to NULL for steady-state problems
     */
    void 
    evaluateSGResponseGradients(
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      Stokhos::EpetraVectorOrthogPoly* sg_g,
      Stokhos::EpetraMultiVectorOrthogPoly* sg_dg_dx,
      Stokhos::EpetraMultiVectorOrthogPoly* sg_dg_dxdot,
      const Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > >& sg_dg_dp);

    //! Compute global residual for multi-point problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalMPResidual(
      const Stokhos::ProductEpetraVector* mp_xdot,
      const Stokhos::ProductEpetraVector& mp_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& mp_p_index,
      const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
      Stokhos::ProductEpetraVector& mp_f);

    //! Compute global Jacobian for multi-point problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalMPJacobian(
      double alpha, double beta,
      const Stokhos::ProductEpetraVector* mp_xdot,
      const Stokhos::ProductEpetraVector& mp_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& mp_p_index,
      const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
      Stokhos::ProductEpetraVector* mp_f,
      Stokhos::ProductEpetraOperator& mp_jac);

    //! Compute global Tangent for multi-point problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalMPTangent(
      double alpha, double beta, bool sum_derivs,
      const Stokhos::ProductEpetraVector* mp_xdot,
      const Stokhos::ProductEpetraVector& mp_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& mp_p_index,
      const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
      const Teuchos::RCP<ParamVec>& deriv_p,
      const Epetra_MultiVector* Vx,
      const Teuchos::SerialDenseMatrix<int,double>* Vp,
      Stokhos::ProductEpetraVector* mp_f,
      Stokhos::ProductEpetraMultiVector* mp_JVx,
      Stokhos::ProductEpetraMultiVector* mp_fVp);

    //! Evaluate multi-point response functions
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void 
    evaluateMPResponses(
      const Stokhos::ProductEpetraVector* mp_xdot,
      const Stokhos::ProductEpetraVector& mp_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& mp_p_index,
      const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
      Stokhos::ProductEpetraVector& mp_g);

    //! Evaluate tangent = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
    /*!
     * Set xdot, dxdot_dp to NULL for steady-state problems
     */
    void 
    evaluateMPResponseTangents(
      const Stokhos::ProductEpetraVector* mp_xdot,
      const Stokhos::ProductEpetraVector& mp_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& mp_p_index,
      const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
      Stokhos::ProductEpetraVector* mp_g,
      const Teuchos::Array< Teuchos::RCP< Stokhos::ProductEpetraMultiVector > >& mp_gt);

    //! Evaluate gradient = dg/dx, dg/dxdot, dg/dp
    /*!
     * Set xdot, dg_dxdot to NULL for steady-state problems
     */
    void 
    evaluateMPResponseGradients(
      const Stokhos::ProductEpetraVector* mp_xdot,
      const Stokhos::ProductEpetraVector& mp_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& mp_p_index,
      const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      Stokhos::ProductEpetraVector* mp_g,
      Stokhos::ProductEpetraMultiVector* mp_dg_dx,
      Stokhos::ProductEpetraMultiVector* mp_dg_dxdot,
      const Teuchos::Array< Teuchos::RCP< Stokhos::ProductEpetraMultiVector > >& mp_dg_dp);

  private:
    
    //! Private to prohibit copying
    Application(const Application&);

    //! Private to prohibit copying
    Application& operator=(const Application&);

  protected:

    //! Parameter list
    Teuchos::RCP<Teuchos::ParameterList> params;
    
    //! Element discretization
    Teuchos::RCP<FEApp::AbstractDiscretization> disc;

    //! Boundary conditions
    std::vector< Teuchos::RCP<FEApp::NodeBC> > bc;

    //! Quadrature rule
    Teuchos::RCP<const FEApp::AbstractQuadrature> quad;

    //! PDE equations
    FEApp::AbstractPDE_TemplateManager<EvalTypes> pdeTM;

    //! Initial solution vector
    Teuchos::RCP<Epetra_Vector> initial_x;

    //! Importer for overlapped data
    Teuchos::RCP<Epetra_Import> importer;

    //! Exporter for overlapped data
    Teuchos::RCP<Epetra_Export> exporter;

    //! Overlapped solution vector
    Teuchos::RCP<Epetra_Vector> overlapped_x;

    //! Overlapped time derivative vector
    Teuchos::RCP<Epetra_Vector> overlapped_xdot;

    //! Overlapped residual vector
    Teuchos::RCP<Epetra_Vector> overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> overlapped_jac;

    //! Parameter library
    Teuchos::RCP<ParamLib> paramLib;

    //! Response functions
    std::vector< Teuchos::RCP<FEApp::AbstractResponseFunction> > responses;

    //! Map for combined response functions
    Teuchos::RCP<Epetra_Map> response_map;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Stochastic Galerking quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > sg_quad;

    //! Stochastic Galerking expansion
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > sg_expansion;

    //! Parallel SG data
    Teuchos::RCP<const Stokhos::ParallelData> sg_parallel_data;

    //! Overlapped SG map
    Teuchos::RCP<const Epetra_BlockMap> sg_overlap_map;

    //! Overlapped SG map for Jacobian
    Teuchos::RCP<const Epetra_BlockMap> sg_overlap_jac_map;

    //! SG overlapped solution vector
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly >  sg_overlapped_x;

    //! SG overlapped time derivative vector
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > sg_overlapped_xdot;

    //! SG overlapped residual vector
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > sg_overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_CrsMatrix> > sg_overlapped_jac;

    //! SG Residual global fill object
    Teuchos::RCP< FEApp::GlobalFill<FEApp::SGResidualType> > sg_res_global_fill;

    //! SG Jacobian global fill object
    Teuchos::RCP< FEApp::GlobalFill<FEApp::SGJacobianType> > sg_jac_global_fill;

    //! MP overlapped solution vector
    Teuchos::RCP< Stokhos::ProductEpetraVector >  mp_overlapped_x;

    //! MP overlapped time derivative vector
    Teuchos::RCP< Stokhos::ProductEpetraVector > mp_overlapped_xdot;

    //! MP overlapped residual vector
    Teuchos::RCP< Stokhos::ProductEpetraVector > mp_overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RCP< Stokhos::ProductContainer<Epetra_CrsMatrix> > mp_overlapped_jac;

    //! MP Residual global fill object
    Teuchos::RCP< FEApp::GlobalFill<FEApp::MPResidualType> > mp_res_global_fill;

    //! MP Jacobian global fill object
    Teuchos::RCP< FEApp::GlobalFill<FEApp::MPJacobianType> > mp_jac_global_fill;

  };

}

#endif // FEAPP_APPLICATION_HPP
