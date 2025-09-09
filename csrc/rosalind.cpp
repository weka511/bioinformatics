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
 *
 * Implementation of the Barnes Hut algorithm to simulate the evolution of a galaxy.
 */

#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include <stdexcept>
#include <sstream>

#include "rosalind.hpp"

using namespace std;

struct option long_options[] ={ 
	{"problem", required_argument, NULL, 'p'},
	{NULL, 0, NULL, 0}
};

int main(int argc, char **argv) {
	ProblemFactory factory;
	char ch;
	string problem_name = "????";
	while ((ch = getopt_long(argc, argv, "p:", long_options, NULL)) != -1)
	  switch (ch)    {
		 case 'p':
			 problem_name = optarg; 
			 break;
		default:
			cerr << ch << endl;
			exit(EXIT_FAILURE);
		}
	try {
		shared_ptr<Problem>  problem = factory.create(problem_name);
		problem->solve();
		return EXIT_SUCCESS;
	}  catch (const exception& e) {
        cerr << __FILE__ << " " << __LINE__ << " Terminating because of errors: "<< endl;
		cerr  << e.what() << endl;
		return EXIT_FAILURE;
    }
}

shared_ptr<Problem> ProblemFactory::create(string problem_name) {
	if (problem_name == "DNA") 
		return make_shared<DNA>();
	 else {
		stringstream message;
		message<<__FILE__ <<" " <<__LINE__<<" Error: Unable to create problem " << problem_name<<endl; 
		throw logic_error(message.str().c_str()); 
	}
} 