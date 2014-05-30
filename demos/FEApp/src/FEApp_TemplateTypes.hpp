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

#ifndef FEAPP_TEMPLATETYPES_HPP
#define FEAPP_TEMPLATETYPES_HPP

#include "FEApp_ConfigDefs.h"

// Include the MPL vector
#include "Sacado_mpl_vector.hpp"

// Include all of our AD types
#include "Sacado_MathFunctions.hpp"
#include "Stokhos_Sacado_MathFunctions.hpp"
#include "Sacado.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Sacado_ETV_Vector.hpp"

// Include ScalarParameterLibrary to specialize traits
#include "Sacado_ScalarParameterLibrary.hpp"

// Typedef AD types to standard names
typedef double RealType;
typedef Sacado::ELRFad::DFad<double> FadType;
typedef Stokhos::StandardStorage<int,double> StorageType;
typedef Sacado::PCE::OrthogPoly<double,StorageType> SGType;
typedef Sacado::ELRCacheFad::DFad<SGType > SGFadType;
typedef Sacado::ETV::Vector<double,StorageType> MPType;
typedef Sacado::ELRCacheFad::DFad<MPType> MPFadType;

// Set up evaluation traits
namespace FEApp {

  // Evaluation types
  struct ResidualType {};
  struct JacobianType {};
  struct TangentType {};

  struct SGResidualType {};
  struct SGJacobianType {};
  struct SGTangentType {};
  
  struct MPResidualType {};
  struct MPJacobianType {};
  struct MPTangentType {};

  // Evaluation traits
  struct EvaluationTraits {
    template <class EvalT> struct apply {};
  };

  template <> struct EvaluationTraits::apply<ResidualType> { 
    typedef RealType type; };
  template <> struct EvaluationTraits::apply<JacobianType> { 
    typedef FadType type; };
  template <> struct EvaluationTraits::apply<TangentType> { 
    typedef FadType type; };
  template <> struct EvaluationTraits::apply<SGResidualType> { 
    typedef SGType type; };
  template <> struct EvaluationTraits::apply<SGJacobianType> { 
    typedef SGFadType  type; };
  template <> struct EvaluationTraits::apply<SGTangentType> { 
    typedef SGFadType type; };
  template <> struct EvaluationTraits::apply<MPResidualType> { 
    typedef MPType type; };
  template <> struct EvaluationTraits::apply<MPJacobianType> { 
    typedef MPFadType  type; };
  template <> struct EvaluationTraits::apply<MPTangentType> { 
    typedef MPFadType type; };
}

// Synonym for the ScalarParameterLibrary/Vector on our traits
typedef Sacado::ScalarParameterLibrary<FEApp::EvaluationTraits> ParamLib;
typedef Sacado::ScalarParameterVector<FEApp::EvaluationTraits> ParamVec;

// Build the MPL vector containing our valid evaluation types
typedef Sacado::mpl::vector<> EvalTypes0;
typedef Sacado::mpl::push_back<EvalTypes0, FEApp::ResidualType>::type EvalTypes1;
typedef Sacado::mpl::push_back<EvalTypes1, FEApp::JacobianType>::type EvalTypes2;
typedef Sacado::mpl::push_back<EvalTypes2, FEApp::TangentType>::type EvalTypes3;
typedef Sacado::mpl::push_back<EvalTypes3, FEApp::SGResidualType>::type EvalTypes4;
typedef Sacado::mpl::push_back<EvalTypes4, FEApp::SGJacobianType>::type EvalTypes5;
typedef Sacado::mpl::push_back<EvalTypes5, FEApp::SGTangentType>::type EvalTypes6;
typedef Sacado::mpl::push_back<EvalTypes6, FEApp::MPResidualType>::type EvalTypes7;
typedef Sacado::mpl::push_back<EvalTypes7, FEApp::MPJacobianType>::type EvalTypes8;
typedef Sacado::mpl::push_back<EvalTypes8, FEApp::MPTangentType>::type EvalTypes9;
typedef EvalTypes9 EvalTypes;

// Turn on/off explicit template instantiation
#define SACADO_ETI

// Define macro for explicit template instantiation
#define INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name) template class name<FEApp::ResidualType>;
#define INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name) template class name<FEApp::JacobianType>;
#define INSTANTIATE_TEMPLATE_CLASS_TANGENT(name) template class name<FEApp::TangentType>;
#define INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL(name) template class name<FEApp::SGResidualType>;
#define INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN(name) template class name<FEApp::SGJacobianType>;
#define INSTANTIATE_TEMPLATE_CLASS_SGTANGENT(name) template class name<FEApp::SGTangentType>;
#define INSTANTIATE_TEMPLATE_CLASS_MPRESIDUAL(name) template class name<FEApp::MPResidualType>;
#define INSTANTIATE_TEMPLATE_CLASS_MPJACOBIAN(name) template class name<FEApp::MPJacobianType>;
#define INSTANTIATE_TEMPLATE_CLASS_MPTANGENT(name) template class name<FEApp::MPTangentType>;

#define INSTANTIATE_TEMPLATE_CLASS(name)      \
  INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name)   \
  INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name)   \
  INSTANTIATE_TEMPLATE_CLASS_TANGENT(name)    \
  INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL(name) \
  INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN(name) \
  INSTANTIATE_TEMPLATE_CLASS_SGTANGENT(name) \
  INSTANTIATE_TEMPLATE_CLASS_MPRESIDUAL(name) \
  INSTANTIATE_TEMPLATE_CLASS_MPJACOBIAN(name) \
  INSTANTIATE_TEMPLATE_CLASS_MPTANGENT(name)

#endif // FEAPP_TEMPLATETYPES_HPP
