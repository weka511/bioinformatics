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
#include "file-adapter.hpp"

using namespace std;


int main(int argc, char **argv) {
	try {
		Parameters parameters(argc,argv);
		ProblemFactory factory;
		shared_ptr<Problem>  problem = factory.create(parameters.get_problem_name());
		FileDatasource datasource("C:\\Users\\Weka\\Downloads\\rosalind_dna_1_dataset.txt");
		FileOutput output("foo.txt");
		problem->attach(&datasource);
		problem->attach(&output);
		problem->solve();
		return EXIT_SUCCESS;
	}  catch (const exception& e) {
        cerr << __FILE__ << " " << __LINE__ << " Terminating because of errors: "<< endl;
		cerr  << e.what() << endl;
		return EXIT_FAILURE;
    }
}

struct option long_options[] ={ 
	{"problem", required_argument, NULL, 'p'},
	{NULL, 0, NULL, 0}
};


Parameters::Parameters(int argc, char **argv){
	char ch;
	while ((ch = getopt_long(argc, argv, "p:", long_options, NULL)) != -1)
		switch (ch)    {
		 case 'p':
			 _problem_name = optarg; 
			 break;
		default:
			cerr << ch << endl;
			stringstream message;
			message<<__FILE__ <<" " <<__LINE__<<" Error: unrecognized parameter " << ch; 
			throw logic_error(message.str().c_str()); 
	}
}


