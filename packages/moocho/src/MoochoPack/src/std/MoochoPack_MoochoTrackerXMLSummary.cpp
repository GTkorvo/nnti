// ////////////////////////////////////////////////////////////////////////////
// MoochoTrackerXMLSummary.cpp
//
// Copyright (C) 2001
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "MoochoPack/src/std/MoochoTrackerXMLSummary.hpp"
#include "MoochoPack/src/NLPAlgoState.hpp"
#include "MoochoPack/src/moocho_algo_conversion.hpp"
#include "NLPInterfacePack/src/abstract/interfaces/NLPFirstOrder.hpp"
//#include "AbstractLinAlgPack/src/abstract/interfaces/Vector.hpp"
//#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixSymOp.hpp"
#include "Teuchos_dyn_cast.hpp"

using std::endl;
using std::setw;

namespace MoochoPack {

MoochoTrackerXMLSummary::MoochoTrackerXMLSummary(
  const Teuchos::RefCountPtr<std::ostream> &journal_out
  ,const std::string xml_filename
  ,const std::string problem_name
  ,const std::string algorithm_description
  )
	:AlgorithmTracker(journal_out)
  ,obj_value_(0.0)
  ,c_norm_value_(0.0)
  ,xml_filename_(xml_filename)
  ,problem_name_(problem_name)
  ,algorithm_description_(algorithm_description)
{
}	


void MoochoTrackerXMLSummary::output_iteration(const Algorithm& algo) const
{

	using Teuchos::dyn_cast;

	const NLPAlgo         &_algo   = rsqp_algo(algo);
	const NLPAlgoState    &s       = _algo.rsqp_state();
	const NLPObjGrad      &nlp     = dyn_cast<const NLPObjGrad>(_algo.nlp()); 

	if (s.k() == 0) {
		// first iteration...
		// write out a dummy file that will be overwritten if the algorithm completes
		// this way there will be an xml file even if the algorithm exceptions...
		Teuchos::RefCountPtr<std::ofstream> xml_out = Teuchos::rcp( new std::ofstream(xml_filename_.c_str(), std::ios::out));
		std::ofstream &out = *xml_out.get();
		
		// Print the Problem Element
		open_problem_element(out, algo);
		
		char ind[] = "   ";

		out << ind << "<Solution status=\"UNKNOWN_FAILURE\" objective_value=\"?\" constraint_norm=\"?\"/>\n";

		out << ind << "<Iterations overall=\"?\"/>\n";

		close_problem_element(out);

		out.close();
	}
}

void MoochoTrackerXMLSummary::output_final(const Algorithm& algo
	, EAlgoReturn algo_return) const
{
	using Teuchos::dyn_cast;

	const NLPAlgo            &_algo  = rsqp_algo(algo);
	const NLPAlgoState           &s      =_algo.rsqp_state();
	const NLPObjGrad      &nlp    = dyn_cast<const NLPObjGrad>(_algo.nlp()); 
	const NLPFirstOrder  *nlp_foi = dynamic_cast<const NLPFirstOrder*>(&nlp); 

	const size_type
		m = nlp.m();

	Teuchos::RefCountPtr<std::ofstream> xml_out = Teuchos::rcp( new std::ofstream(xml_filename_.c_str(), std::ios::out));
	std::ofstream &out = *xml_out.get();

	// Print the Problem Element
	open_problem_element(out, algo);

	char ind[] = "   ";

	char soln_status[256];
	switch (algo_return) 
		{
		case IterationPack::TERMINATE_TRUE:
			strcpy(soln_status, "SOLVED");
			break;
		case IterationPack::TERMINATE_FALSE:
			strcpy(soln_status, "FAILED");
			break;
		case IterationPack::MAX_ITER_EXCEEDED:
			strcpy(soln_status, "MAX_ITER");
			break;
		case IterationPack::MAX_RUN_TIME_EXCEEDED:
			strcpy(soln_status, "MAX_RUN_TIME");
			break;
		case IterationPack::INTERRUPTED_TERMINATE_TRUE:
			strcpy(soln_status, "INTERRUPTED_SOLVED");
			break;
		case IterationPack::INTERRUPTED_TERMINATE_FALSE:
			strcpy(soln_status, "INTERRUPTED_FAILED");
			break;
		default:
			strcpy(soln_status, "UNKNOWN_STATUS");
			break;
		}

	// Output the solution element
	out << ind << "<Solution status=\"" << soln_status;

	if (s.f().updated_k(0)) {
		out << "\" objective_value=\"" << s.f().get_k(0) << "\"";
	}
	else {
		out << "\" objective_value=\"?\"";		
	}

	if (m && s.c().updated_k(0)) {
		out << " constraint_norm=\"" << s.c().get_k(0).norm_inf() << "\"";
	}
	else {
		out << " constraint_norm=\"?\"";
	}
		
	out	<< "/>\n";

	// Output the Iterations element
	out << ind << "<Iterations overall=\"" << s.k() << "\"/>\n";

	// Output the Evaluations element
	out << ind << "<Evaluations>\n"
		<< ind << ind << "<Objective evaluations=\"" << nlp.num_f_evals() << "\"/>\n"
	    << ind << ind << "<Constraint evaluations=\"" << ( m ? nlp.num_c_evals() : 0 ) << "\"/>\n"
	    << ind << ind << "<Objective_Gradient evaluations=\"" << nlp.num_Gf_evals() << "\"/>\n"
		<< ind << ind << "<Constraint_Gradient evaluations=\"";

	if(m) {
		if (nlp_foi) {
			out << nlp_foi->num_Gc_evals();
		}
		else {
			out << "?";
		}
	}
	else {
		out << 0;
	}

	out << "\"/>\n"
		<< ind << "</Evaluations>\n";

	// Output the Timing element
	/*	out << ind << "<Timing>\n";
	const int n = _algo.num_steps();
	const int final_iter = s.k();
	const int n_iters = final_iter+1;
	std::vector<double>* step_times = new std::vector<double>[n_iters](n+1);
	for (int k=0; k<n_iters; k++) {
		_algo.get_step_times_k(k-(final_iter), &(step_times[k])[0]);
		}

	for (int i=0; i<n+1; i++) {
		if (i != n) {
			out << ind << ind << "<Step name=\"" << _algo.get_step_name(i+1) << "\">\n";
		}
		else {
			// overall step
		out << ind << ind << "<Overall>\n";
			
		}

			
		for (int k=0; k<final_iter; k++) {
			out << ind << ind << ind << "<IterationTime iteration=\"" << k << "\" seconds=\"" << (step_times[k])[i] << "\"/>\n";
		}
		double total, average, min, max, percent;
		_algo.get_final_step_stats(i, &total, &average, &min, &max, &percent);
		
		out << ind << ind << ind << "<TotalTime seconds=\"" << total << "\" percent=\"" << percent*100.0 << "\"/>\n"
			<< ind << ind << ind << "<AverageTime seconds=\"" << average << "\"/>\n"
			<< ind << ind << ind << "<MinTime seconds=\"" << min << "\"/>\n"
			<< ind << ind << ind << "<MaxTime seconds=\"" << max << "\"/>\n";

		if (i != n) {
			out << ind << ind << "</Step>\n";
		}
		else {
			// overall step
		out << ind << ind << "</Overall>\n";
			
		}
	}

	delete [] step_times;
			
	out	<< ind << "</Timing>\n";

	*/
	// Close the problem element
	close_problem_element(out);

	out.close();
}


void MoochoTrackerXMLSummary::output_pre_file() const
{
	Teuchos::RefCountPtr<std::ofstream> xml_out = Teuchos::rcp( new std::ofstream(xml_filename_.c_str(), std::ios::out));
	std::ofstream &out = *xml_out.get();
	
	char ind[] = "   ";
	
	// get a string representation of the current date/time
	time_t current_time = time(NULL);
	char time_str[26];
	strcpy(time_str, ctime(&current_time));
	time_str[24]='\0';
	out << "<Problem name=\"" << problem_name_  << "\" time=\"" << time_str << "\">\n";

	out << ind << "<Dimension n=\"?\" m=\"?\"/>\n";
	
	// get the machine name - NOTE: this is not portable, may need to
	// look into a way to do this on multiple platforms
	char machine_name[256];
	strcpy(machine_name, "UnknownMachine");
	FILE* machine_file = popen("uname -n", "r");
	if (machine_file) {
		fread(machine_name, sizeof(char), 255, machine_file);
		char* end = strchr(machine_name, '\n');
		*end = '\0';
		pclose(machine_file);
	}
	
	out << ind << "<Machine name=\"" << machine_name << "\"/>\n";

	out << ind << "<Solution status=\"UNKNOWN_FAILURE\" objective_value=\"?\" constraint_norm=\"?\"/>\n";
	
	out << ind << "<Algorithm description=\"" << algorithm_description_ << "\"/>\n";

	out << ind << "<Iterations overall=\"?\"/>\n";

	out << "</Problem>\n";
		
	out.close();	

}

void MoochoTrackerXMLSummary::open_problem_element( std::ostream& out, const Algorithm& algo) const
{
	if (out) {
		const NLPAlgo  &_algo  = rsqp_algo(algo);
		const NLP      &nlp    = _algo.nlp(); 

		char ind[] = "   ";

		// get a string representation of the current date/time
		time_t current_time = time(NULL);
		char time_str[26];
		strcpy(time_str, ctime(&current_time));
		time_str[24]='\0';
		out << "<Problem name=\"" << problem_name_  << "\" time=\"" << time_str << "\">\n";

		out << ind << "<Dimension n=\"" << nlp.n() << "\" m=\"" << nlp.m() << "\"/>\n";

		// get the machine name - NOTE: this is not portable, may need to
		// look into a way to do this on multiple platforms
		char machine_name[256];
		strcpy(machine_name, "UnknownMachine");
		FILE* machine_file = popen("uname -n", "r");
		if (machine_file) {
			fread(machine_name, sizeof(char), 255, machine_file);
			char* end = strchr(machine_name, '\n');
			*end = '\0';
			pclose(machine_file);
		}

		out << ind << "<Machine name=\"" << machine_name << "\"/>\n";
		//out << ind << "<File location=\"" << file_location_ << "\"/>\n";
		
		out << ind << "<Algorithm description=\"" << algorithm_description_ << "\">\n";
		// output some description of the algorithm and 
		// its options
		out << ind << "</Algorithm>\n";
	}

}

void MoochoTrackerXMLSummary::close_problem_element( std::ostream& out) const
{
	if (out) {
		out << "</Problem>\n";
	}
}

} // end namespace MoochoPack
