// ////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <cassert>

#include "../../include/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.h"
#include "Misc/include/StringToBool.h"

// Define the options
namespace {

	const int local_num_options = 4;

	const char options_group_name[] = "MeritFuncPenaltyParamUpdate";

	enum local_EOptions {
	    SMALL_MU,
	    MIN_MU_RATIO,
	    MULT_FACTOR,
	    KKT_NEAR_SOL
	};

	const char* local_SOptions[local_num_options]	= {
	    "small_mu",
	    "min_mu_ratio",
	    "mult_factor",
	    "kkt_near_sol"
	};

}

namespace ReducedSpaceSQPPack {

MeritFunc_PenaltyParamUpdate_AddedStepSetOptions::MeritFunc_PenaltyParamUpdate_AddedStepSetOptions(
			MeritFunc_PenaltyParamUpdate_AddedStep* target )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  options_group_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			MeritFunc_PenaltyParamUpdate_AddedStep >( target )
{}

void MeritFunc_PenaltyParamUpdate_AddedStepSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	switch( (local_EOptions)option_num ) {
		case SMALL_MU: {
			target().small_mu( ::atof( option_value.c_str() ) );
			break;
		}
		case MIN_MU_RATIO: {
			target().min_mu_ratio( ::atof( option_value.c_str() ) );
			break;
		}
		case MULT_FACTOR: {
			target().mult_factor( ::atof( option_value.c_str() ) );
			break;
		}
		case KKT_NEAR_SOL: {
			target().kkt_near_sol( ::atof( option_value.c_str() ) );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace ReducedSpaceSQPPack 