// ////////////////////////////////////////////////////////////////
// BasisSystemTesterSetOptions.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <assert.h>
#include <math.h>

#include "AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

	const int local_num_options = 6;

	enum local_EOptions {
		PRINT_TESTS
		,DUMP_ALL
		,TEST_FOR_EXCEPTION
		,NUM_RANDOM_TESTS
		,WARNING_TOL
		,ERROR_TOL
	};

	const char* local_SOptions[local_num_options]	= {
		"print_tests"
		,"dump_all"
		,"throw_exception"
		,"num_random_tests"
		,"warning_tol"
		,"error_tol"
	};

}

namespace AbstractLinAlgPack {

BasisSystemTesterSetOptions::BasisSystemTesterSetOptions(
			  BasisSystemTester* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			BasisSystemTester >( target )
{}

void BasisSystemTesterSetOptions::setOption(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;
	typedef BasisSystemTester target_t;
	switch( (local_EOptions)option_num ) {
		case PRINT_TESTS:
		{
			const std::string &option = option_value.c_str();
			if( option == "PRINT_NONE" )
				target().print_tests( target_t::PRINT_NONE );
		    else if( option == "PRINT_BASIC" )
				target().print_tests( target_t::PRINT_BASIC );
		    else if( option == "PRINT_MORE" )
				target().print_tests( target_t::PRINT_MORE );
		    else if( option == "PRINT_ALL" )
				target().print_tests( target_t::PRINT_ALL );
			else
				TEST_FOR_EXCEPTION(
					true, std::invalid_argument
					,"Error, incorrect value for "
					"\"print_tests\".  Only the options "
					"PRINT_NONE, PRINT_BASIS, PRINT_MORE and PRINT_ALL are allowed" );
			break;
		}
		case DUMP_ALL:
			target().dump_all(
				StringToBool( "dump_all", option_value.c_str() )
				);
			break;
		case TEST_FOR_EXCEPTION:
			target().throw_exception(
				StringToBool( "throw_exception", option_value.c_str() )
				);
			break;
	    case NUM_RANDOM_TESTS:
			target().num_random_tests(::abs(::atoi(option_value.c_str())));
			break;
	    case WARNING_TOL:
			target().warning_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case ERROR_TOL:
			target().error_tol(::fabs(::atof(option_value.c_str())));
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace AbstractLinAlgPack
