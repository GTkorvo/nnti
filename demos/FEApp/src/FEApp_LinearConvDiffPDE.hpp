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

#ifndef FEAPP_LINEARCONVDIFFPDE_HPP
#define FEAPP_LINEARCONVDIFFPDE_HPP

#include "Teuchos_RCP.hpp"

#include "FEApp_AbstractPDE.hpp"
#include "FEApp_FunctionFactory.hpp"
#include "FEApp_SourceFunctionFactory.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterEntry.hpp"

namespace FEApp {

  template <typename EvalT>
  class LinearConvDiffPDE : public FEApp::AbstractPDE<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractPDE<EvalT>::ScalarT ScalarT;
  
    //! Constructor
    LinearConvDiffPDE(
       const Teuchos::RCP< const FEApp::AbstractSourceFunction<EvalT> >& src_func,
       const Teuchos::RCP< const FEApp::AbstractFunction<EvalT> >& cnv_func,
       double viscosity,double scalingParameter);

    //! Destructor
    virtual ~LinearConvDiffPDE();

    //! Number of discretized equations
    virtual unsigned int numEquations() const;

    //! Initialize PDE
    virtual void init(unsigned int numQuadPoints, unsigned int numNodes);

    //! Evaluate discretized PDE element-level residual
    virtual void
    evaluateElementResidual(const FEApp::AbstractQuadrature& quadRule,
                            const FEApp::AbstractElement& element,
                            const std::vector<ScalarT>* dot,
                            const std::vector<ScalarT>& solution,
                            std::vector<ScalarT>& residual);

    Teuchos::RCP<Sacado::ScalarParameterEntry<EvalT,EvaluationTraits> > 
    getViscosityPertParameter();
 
  private:

    //! Private to prohibit copying
    LinearConvDiffPDE(const LinearConvDiffPDE&);

    //! Private to prohibit copying
    LinearConvDiffPDE& operator=(const LinearConvDiffPDE&);

  protected:

    //! Pointer to source function
    Teuchos::RCP< const FEApp::AbstractSourceFunction<EvalT> > source;

    //! pointer to convection function
    Teuchos::RCP< const FEApp::AbstractFunction<EvalT> > convection;

    //! Number of quad points
    unsigned int num_qp;
    
    //! Number of nodes
    unsigned int num_nodes;

    double baseViscosity;
    double scalingParameter;
    ScalarT viscPert;

    //! Shape function values
    std::vector< std::vector<double> > phi;

    //! Shape function derivatives
    std::vector< std::vector<double> > dphi;

    //! Element transformation Jacobian
    std::vector<double> jac;

    //! Coordinates of quadrature points
    std::vector<double> x;

    //! Discretized solution
    std::vector<ScalarT> u;

    //! Discretized solution
    std::vector<ScalarT> du;

    //! Discretized time derivative
    std::vector<ScalarT> udot;

    //! Source function values
    std::vector<ScalarT> f;

    //! Convection function values
    std::vector<ScalarT> conv;
  };

  class LinearConvDiffPDE_TemplateBuilder {
  public:
    LinearConvDiffPDE_TemplateBuilder(
		const Teuchos::RCP<Teuchos::ParameterList>& params_,
	        const Teuchos::RCP<ParamLib>& paramLib) :
       src_params(Teuchos::rcp(&(params_->sublist("Source Function")),false)),
       cnv_params(Teuchos::rcp(&(params_->sublist("Convection Function")),false)),
       baseViscosity(params_->get<double>("Base Viscosity",1.0)),
       scalingParameter(params_->get<double>("Perturbation Scaling",1.0)),
       pl(paramLib) {}

    template <typename T>
    Teuchos::RCP<FEApp::AbstractPDE_NTBase> build() const 
    {
       FEApp::SourceFunctionFactory<T> srcFactory(src_params, pl);
       FEApp::FunctionFactory<T> cnvFactory(cnv_params, pl);

       Teuchos::RCP< FEApp::AbstractSourceFunction<T> > source = srcFactory.create();
       Teuchos::RCP< FEApp::AbstractFunction<T> > convection = cnvFactory.create();

       Teuchos::RCP<FEApp::LinearConvDiffPDE<T> > pde 
          = Teuchos::rcp( new FEApp::LinearConvDiffPDE<T>(source,convection,baseViscosity,scalingParameter));

       // add viscosity perturbation parameter to the library
       std::string name = "Viscosity Perturbation";
       if(!pl->isParameter(name))
          pl->addParameterFamily(name, true, false);
       if(!pl->template isParameterForType<T>(name))
          pl->template addEntry<T>(name, pde->getViscosityPertParameter());

       return pde;
    }

  protected:
    Teuchos::RCP<Teuchos::ParameterList> src_params;
    Teuchos::RCP<Teuchos::ParameterList> cnv_params;
    double baseViscosity;
    double scalingParameter;
    Teuchos::RCP<ParamLib> pl;
  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_LinearConvDiffPDEImpl.hpp"
#endif 

#endif // FEAPP_HEATNONLINERASOURCEPDE_HPP
