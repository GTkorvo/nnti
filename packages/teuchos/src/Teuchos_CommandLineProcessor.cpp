// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// //////////////////////////////////////////////////
// Teuchos_CommandLineProcessor.cpp

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

namespace {

inline int my_max( int a, int b ) { return a > b ? a : b; }

std::string remove_quotes( const std::string& str )
{
	if(str[0] != '\"')
		return str;
	return str.substr(1,str.size()-2);
}

std::string add_quotes( const std::string& str )
{
	if(str[0] == '\"')
		return str;
	return "\"" + str + "\"";
}

} // end namespace

namespace Teuchos {

CommandLineProcessor::CommandLineProcessor(
	bool   throwExceptions
	,bool  recogniseAllOptions
	)
	:throwExceptions_(throwExceptions)
	,recogniseAllOptions_(recogniseAllOptions)
{}

// Set up options

void CommandLineProcessor::setOption(
	const char     option_true[]
	,const char    option_false[]
	,bool          *option_val
	,const char    documentation[]
	)
{
	assert(option_val!=NULL);
	options_list_[std::string(option_true)]
		= opt_val_val_t(OPT_BOOL_TRUE,option_val);
	options_list_[std::string(option_false)]
		= opt_val_val_t(OPT_BOOL_FALSE,option_val);
	options_documentation_list_.push_back(
		opt_doc_t(OPT_BOOL_TRUE,option_true,option_false,std::string(documentation?documentation:""),option_val) );
}

void CommandLineProcessor::setOption(
	const char     option_name[]
	,int           *option_val
	,const char    documentation[]
	)
{
	assert(option_val!=NULL);
	options_list_[std::string(option_name)]
		= opt_val_val_t(OPT_INT,option_val);
	options_documentation_list_.push_back(
		opt_doc_t(OPT_INT,option_name,"",std::string(documentation?documentation:""),option_val) );
}

void CommandLineProcessor::setOption(
	const char     option_name[]
	,double        *option_val
	,const char    documentation[]
	)
{
	assert(option_val!=NULL);
	options_list_[std::string(option_name)]
		= opt_val_val_t(OPT_DOUBLE,option_val);
	options_documentation_list_.push_back(
		opt_doc_t(OPT_DOUBLE,option_name,"",std::string(documentation?documentation:""),option_val) );
}

void CommandLineProcessor::setOption(
	const char     option_name[]
	,std::string   *option_val
	,const char    documentation[]
	)
{
	assert(option_val!=NULL);
	options_list_[std::string(option_name)]
		= opt_val_val_t(OPT_STRING,option_val);
	options_documentation_list_.push_back(
		opt_doc_t(OPT_STRING,option_name,"",std::string(documentation?documentation:""),option_val) );
}

// Parse command line

CommandLineProcessor::EParseCommandLineReturn
CommandLineProcessor::parse(
	int             argc
	,char*          argv[]
	,std::ostream   *errout
	) const
{
	std::string        opt_name;
	std::string        opt_val_str;
	const std::string  help_opt = "help";
	const std::string  pause_opt = "pause-for-debugging";
#ifdef HAVE_MPI
	int procRank = -1;
	MPI_Comm_rank( MPI_COMM_WORLD, &procRank );
#endif
	for( int i = 1; i < argc; ++i ) {
		bool gov_return = get_opt_val( argv[i], &opt_name, &opt_val_str );
		if( !gov_return && recogniseAllOptions() ) {
			print_bad_opt(i,argv,errout);
			return PARSE_UNRECOGNIZED_OPTION;
		}
		if( opt_name == help_opt ) {
			print_help_msg( argc, argv, errout );
			return PARSE_HELP_PRINTED;
		}
		if( opt_name == pause_opt ) {
#ifdef HAVE_MPI
			if(procRank==0) {
#endif
				std::cerr << "\nType 0 and press enter to continue : ";
				int dummy_int = 0;
				std::cin >> dummy_int;
#ifdef HAVE_MPI
			}
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			continue;
		}
		// Lookup the option (we had better find it!)
		options_list_t::const_iterator  itr = options_list_.find(opt_name);
		if( itr == options_list_.end() && recogniseAllOptions() ) {
 			print_bad_opt(i,argv,errout);
			return PARSE_UNRECOGNIZED_OPTION;
		}
		// Changed access to second value of map to not use overloaded arrow operator, 
		// otherwise this code will not compile on Janus (HKT, 12/01/2003) 
		const opt_val_val_t &opt_val_val = (*itr).second;
		switch( opt_val_val.opt_type ) {
			case OPT_BOOL_TRUE:
				*((bool*)opt_val_val.opt_val) = true;
				break;
			case OPT_BOOL_FALSE:
				*((bool*)opt_val_val.opt_val) = false;
				break;
			case OPT_INT:
				*((int*)opt_val_val.opt_val) = ::atoi(opt_val_str.c_str());
				break;
			case OPT_DOUBLE:
				*((double*)opt_val_val.opt_val) = ::atof(opt_val_str.c_str());
				break;
			case OPT_STRING:
				*((std::string*)opt_val_val.opt_val) = remove_quotes(opt_val_str);
				break;
			default:
				assert(0); // Local programming error only
		} 
	}
	return PARSE_SUCCESSFUL;
}

// private

bool CommandLineProcessor::get_opt_val(
	const char     str[]
	,std::string   *opt_name
	,std::string   *opt_val_str
	) const
{
	const int len = strlen(str);
	if( len < 3 )
		return false; // Can't be an option with '--' followed by at least one char
	if( str[0] != '-' || str[1] != '-' )
		return false; // Not a recognised option
	// Find the '='
	int equ_i;
	for( equ_i = 2; equ_i < len && str[equ_i] != '='; ++equ_i );
	// Set opt_name
	opt_name->assign( str + 2, equ_i-2 );
	// Set opt_val_str
	if( equ_i == len ) {
		*opt_val_str = "";
	}
	else {
		opt_val_str->assign( str + equ_i + 1, len - equ_i - 1 );
	}
	return true;
}

void CommandLineProcessor::print_help_msg(
	int             argc
	,char*          argv[]
	,std::ostream   *errout
	) const
{
	using std::setw;
	using std::endl;

	if( !errout )
		return;

	const int opt_type_w = 8;
	const char spc_chars[] = "  ";

	// Get the maximum length of an option name
	int opt_name_w = 19; // For the 'pause-for-debugging' option
	options_documentation_list_t::const_iterator itr;
	for( itr = options_documentation_list_.begin(); itr != options_documentation_list_.end(); ++itr ) {
		opt_name_w = my_max(opt_name_w,itr->opt_name.length());
		if( itr->opt_type )
			opt_name_w = my_max(opt_name_w,itr->opt_name_false.length());
	}
	opt_name_w += 2;

	// Print the help message
	*errout
		<< "Usage: " << argv[0] << " [options]\n"
		<< spc_chars << "options:\n"
		<< spc_chars
		<< "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
		<< std::left << setw(opt_name_w) << "help"
		<< std::left << setw(opt_type_w) << " "
#else
	        << std::setiosflags(std::ios::left) << setw(opt_name_w) << "help"
	        << std::setiosflags(std::ios::left) << setw(opt_type_w) << " "
#endif
		<< "Prints this help message"
		<< endl
		<< spc_chars
		<< "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
		<< std::left << setw(opt_name_w) << "pause-for-debugging"
		<< std::left << setw(opt_type_w) << " "
#else
	        << std::setiosflags(std::ios::left) << setw(opt_name_w) << "pause-for-debugging"
	        << std::setiosflags(std::ios::left) << setw(opt_type_w) << " "
#endif
		<< "Pauses for user input to allow attaching a debugger"
		<< endl;
	for( itr = options_documentation_list_.begin(); itr != options_documentation_list_.end(); ++itr ) {
		*errout
			<< spc_chars
			<< "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
			<< std::left << setw(opt_name_w) << itr->opt_name
			<< std::left << setw(opt_type_w) << opt_type_str(itr->opt_type)
#else
			<< std::setiosflags(std::ios::left) << setw(opt_name_w) << itr->opt_name
			<< std::setiosflags(std::ios::left) << setw(opt_type_w) << opt_type_str(itr->opt_type)
#endif
			<< ( itr->documentation.length() ? itr->documentation.c_str() : "No documentation" )
			<< endl;
		if( itr->opt_type == OPT_BOOL_TRUE ) {
			*errout
				<< spc_chars
				<< "--"
				<< setw(opt_name_w) << itr->opt_name_false;
		}
		else {
			*errout
				<< spc_chars
				<< "  "
				<< setw(opt_name_w) << " ";
		}
		*errout
			<< setw(opt_type_w) << " "
			<< "(default: ";
		switch( itr->opt_type ) {
			case OPT_BOOL_TRUE:
				*errout << "--" << ( (*((bool*)itr->default_val)) ? itr->opt_name : itr->opt_name_false );
				break;
			case OPT_INT:
			case OPT_DOUBLE:
			case OPT_STRING:
				*errout << "--" << itr->opt_name;
				break;
			default:
				assert(0); // Local programming error only
		}		
		switch( itr->opt_type ) {
			case OPT_BOOL_TRUE:
				break;
			case OPT_INT:
				*errout << "=" << (*((int*)itr->default_val));
				break;
			case OPT_DOUBLE:
				*errout <<  "=" << (*((double*)itr->default_val));
				break;
			case OPT_STRING:
				*errout <<  "=" << add_quotes(*((std::string*)itr->default_val));
				break;
			default:
				assert(0); // Local programming error only
		}		
		*errout << ")\n";
	}
	if(throwExceptions_)
		TEST_FOR_EXCEPTION( true, HelpPrinted, "Help message was printed" );
}

void CommandLineProcessor::print_bad_opt(
	int             argv_i
	,char*          argv[]
	,std::ostream   *errout
	) const
{
#   define CLP_ERR_MSG "Error, option " << argv_i-1 << " \'" << argv[argv_i] << "\' is not recognized (use --help)!"
	if(errout)
		*errout << argv[0] << " : " << CLP_ERR_MSG << std::endl;
	if(throwExceptions_)
		TEST_FOR_EXCEPTION( true, UnrecognizedOption, CLP_ERR_MSG );
#   undef CLP_ERR_MSG
}

} // end namespace Teuchos


