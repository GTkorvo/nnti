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

#include <assert.h>
#include <math.h>

#include "MoochoPack_LineSearchWatchDog_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 2;

  const char options_group_name[] = "LineSearchWatchDog";

  enum local_EOptions {
    OPT_KKT_ERR_THESHOLD,
    FEAS_KKT_ERR_THESHOLD
  };

  const char* local_SOptions[local_num_options]	= {
    "opt_kkt_err_threshold",
    "feas_kkt_err_threshold"
  };

}

namespace MoochoPack {

LineSearchWatchDog_StepSetOptions::LineSearchWatchDog_StepSetOptions(
      LineSearchWatchDog_Step* target )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        options_group_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      LineSearchWatchDog_Step >( target )
{}

void LineSearchWatchDog_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  switch( (local_EOptions)option_num ) {
    case OPT_KKT_ERR_THESHOLD: {
      target().opt_kkt_err_threshold( ::fabs( ::atof( option_value.c_str() ) ) );
      break;
    }
    case FEAS_KKT_ERR_THESHOLD: {
      target().feas_kkt_err_threshold( ::fabs( ::atof( option_value.c_str() ) ) );
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
