/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_PARAMETER_H
#define SUNDANCE_PARAMETER_H


#include "SundanceDefs.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceParameterData.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  /** 
   * A Parameter is an expression subtype used to represent
   * a spatially-constant parameter that may change during the
   * course of a simulation, for example, the time in a transient
   * simulation or a continuation parameter when using a homotopy method.
   * While it is possible to use simple double-precision 
   * constants in expressions, their values are immutable once created.
   * When a constant's value may need to be changed, use a Parameter rather
   * than a simple double.
   */
  class Parameter : public virtual DiscreteFuncElement,
                    public virtual SpatiallyConstantExpr
  {
  public:
    /** */
    Parameter(const double& value, const string& name="");

    /** virtual destructor */
    virtual ~Parameter() {;}

    /** */
    virtual XMLObject toXML() const ;

    /** */
    virtual const double& value() const ;
    
    /** */
    virtual void setValue(const double& value) ;

    //@{
    /** \name Functions for internal use */
    /** */
    Set<MultipleDeriv> 
    internalFindW(int order, const EvalContext& context) const ;
          
    /** */
    Set<MultipleDeriv> 
    internalFindV(int order, const EvalContext& context) const ;
    
    /** */
    Set<MultipleDeriv> 
    internalFindC(int order, const EvalContext& context) const ;

    /** */
    RefCountPtr<Array<Set<MultipleDeriv> > > 
    internalDetermineR(const EvalContext& context,
                       const Array<Set<MultipleDeriv> >& RInput) const ;


    /** */
    virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                       const EvalContext& context) const ;

    /** */
    virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}
    //@}

    /** Write self in text form */
    virtual ostream& toText(ostream& os, bool paren) const 
    {os << "Parameter[" << name() << " = " << value() << "]"; return os;}


  protected:
    /** */
    const ParameterData* data() const ;

    /** */
    ParameterData* data();


  };
}
#endif
