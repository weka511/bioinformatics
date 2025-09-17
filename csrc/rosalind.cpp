/**
 * Copyright (C) 2025 Simon Crase
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 */

#include <iostream>
#include <getopt.h>
#include <stdexcept>
#include <sstream>

#include "rosalind.hpp"
#include "factory.hpp"


using namespace std;

/**
 *  A program written to solve rasalind problems.
 */
int main(int argc, char **argv) {
	try {
		Parameters parameters(argc,argv);
		string problem_name = parameters.get_problem_name();
		if (problem_name.length() ==0) {
			cerr << __FILE__ << " " << __LINE__ << " No problem name specified "<< endl;
			return EXIT_FAILURE;
		}
		ProblemFactory factory;
		shared_ptr<Problem>  problem = factory.create(problem_name);
		FileNameFactory fnf;
		string file_name = fnf.create(problem_name,parameters.get_format(),parameters.get_sequence());
		FileDatasource datasource(file_name);
		FileOutput output(parameters.get_output_name());
		problem->attach(&datasource);
		problem->attach(&output);
		problem->solve();
		output.flush();
		return EXIT_SUCCESS;
	}  catch (const exception& e) {
        cerr << __FILE__ << " " << __LINE__ << " Terminating because of errors: "<< endl;
		cerr  << e.what() << endl;
		return EXIT_FAILURE;
    }
}

struct option long_options[] ={ 
	{"problem", required_argument, NULL, 'p'},
	{"test", required_argument, NULL, 't'},
	{"out", required_argument, NULL, 'o'},
	{"sequence", required_argument, NULL, 'q'},
	{NULL, 0, NULL, 0}
};


Parameters::Parameters(int argc, char **argv){
	char ch;
	while ((ch = getopt_long(argc, argv, "p:to:q:", long_options, NULL)) != -1)
		switch (ch)    {
		 case 'p':
			 _problem_name = optarg; 
			 break;
		  case 't':
			_format = FileNameFactory::Format::TEST;
			break;
		case 'q':
			_sequence = atoi(optarg);
			break;
		default:
			cerr << ch << endl;
			stringstream message;
			message<<__FILE__ <<" " <<__LINE__<<" Error: unrecognized parameter " << ch; 
			throw logic_error(message.str().c_str()); 
	}
}


