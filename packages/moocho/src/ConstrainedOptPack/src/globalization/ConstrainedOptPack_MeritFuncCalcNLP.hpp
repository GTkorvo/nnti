// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef MERIT_FUNC_CALC_NLP_H
#define MERIT_FUNC_CALC_NLP_H

#include "ConstrainedOptPack_MeritFuncCalc.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "StandardAggregationMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Adds the ability to compute <tt>phi(f(x),c(x),h(x))</tt> at \c x
  * directly instead of having to compute f, c and h first.
  * This class uses an aggregate NLP to perform the computations of \a f(x)
  * \a c(x) and \a h(x).
  */
class MeritFuncCalcNLP : public MeritFuncCalc {
public:

  /** @name Constructors/initializers */
  //@{

  /// <<std aggr>> stereotype members for phi.
  STANDARD_CONST_AGGREGATION_MEMBERS( MeritFuncNLP, phi )

  /// <<std aggr>> stereotype members for nlp.
  STANDARD_CONST_AGGREGATION_MEMBERS( NLP, nlp )

  /** \brief . */
  MeritFuncCalcNLP( const MeritFuncNLP* phi = 0, const NLP* nlp = 0 );

  //@}

  /** @name Overridden from MeritFuncCalc */
  //@{

  /** \brief Return the value of the merit function at x.
   * Here phi(x) is calculated directly using the nlp.
   */
  value_type operator()(const Vector& x) const;

  /// Calls phi().deriv() on phi.
  value_type deriv() const;

  /// Calls <tt>phi().print_merit_func()</tt>.
  void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

};	// end class MeritFuncCalcNLP

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_NLP_H
