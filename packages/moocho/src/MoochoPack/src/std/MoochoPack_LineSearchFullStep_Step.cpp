// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <ostream>

#include "MoochoPack_LineSearchFullStep_Step.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "Teuchos_Assert.hpp"

namespace MoochoPack {

MoochoPack::LineSearchFullStep_Step::LineSearchFullStep_Step(
    const bounds_tester_ptr_t&	bounds_tester
    )
  :
    bounds_tester_(bounds_tester)
{}


bool LineSearchFullStep_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using LinAlgOpPack::V_VpV;

  NLPAlgo        &algo   = rsqp_algo(_algo);
  NLPAlgoState   &s      = algo.rsqp_state();
  NLP            &nlp    = algo.nlp();

  const size_type
    m = nlp.m();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }
  
  // alpha_k = 1.0
  IterQuantityAccess<value_type>
    &alpha_iq   = s.alpha();
  if( !alpha_iq.updated_k(0) )
    alpha_iq.set_k(0) = 1.0;

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\nf_k        = " << s.f().get_k(0);
    if(m)
      out << "\n||c_k||inf = " << s.c().get_k(0).norm_inf();
    out << "\nalpha_k    = " << alpha_iq.get_k(0) << std::endl;
  }

  // x_kp1 = x_k + d_k
  IterQuantityAccess<VectorMutable>  &x_iq = s.x();
  const Vector                       &x_k   = x_iq.get_k(0);
  VectorMutable                      &x_kp1 = x_iq.set_k(+1);
  x_kp1 = x_k;
  Vp_StV( &x_kp1, alpha_iq.get_k(0), s.d().get_k(0) );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||x_kp1||inf   = " << s.x().get_k(+1).norm_inf() << std::endl;
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nx_kp1 =\n" << s.x().get_k(+1);
  }

  if(algo.algo_cntr().check_results()) {
    assert_print_nan_inf(
      x_kp1, "x_kp1",true
      ,int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
      );
    if( nlp.num_bounded_x() ) {
      if(!bounds_tester().check_in_bounds(
          int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
        , int(olevel) >= int(PRINT_VECTORS)					// print_all_warnings
        , int(olevel) >= int(PRINT_ITERATION_QUANTITIES)	// print_vectors
        , nlp.xl(), "xl"
        , nlp.xu(), "xu"
        , x_kp1, "x_kp1"
        ))
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, TestFailed
          ,"LineSearchFullStep_Step::do_step(...) : Error, "
          "the variables bounds xl <= x_k(+1) <= xu where violated!" );
      }
    }
  }

  // Calcuate f and c at the new point.
  nlp.unset_quantities();
  nlp.set_f( &s.f().set_k(+1) );
  if(m) nlp.set_c( &s.c().set_k(+1) );
  nlp.calc_f(x_kp1);
  if(m) nlp.calc_c( x_kp1, false );
  nlp.unset_quantities();

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\nf_kp1        = " << s.f().get_k(+1);
    if(m)
      out << "\n||c_kp1||inf = " << s.c().get_k(+1).norm_inf();
    out << std::endl;
  }

  if( m && static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nc_kp1 =\n" << s.c().get_k(+1); 
  }

  if(algo.algo_cntr().check_results()) {
    assert_print_nan_inf( s.f().get_k(+1), "f(x_kp1)", true, &out );
    if(m)
      assert_print_nan_inf( s.c().get_k(+1), "c(x_kp1)", true, &out );
  }

  return true;
}

void LineSearchFullStep_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "if alpha_k is not updated then\n"
    << L << "    alpha_k = 1.0\n"
    << L << "end\n"
    << L << "x_kp1 = x_k + alpha_k * d_k\n"
    << L << "f_kp1 = f(x_kp1)\n"
    << L << "if m > 0 then c_kp1 = c(x_kp1)\n";
}

} // end namespace MoochoPack
