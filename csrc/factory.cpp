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
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cctype>

#include "factory.hpp"
#include "dna.hpp"
#include "rna.hpp"

using namespace std;

ProblemFactory::ProblemFactory(){
	problem_map["DNA"] = &createProblem<DNA>;
	problem_map["RNA"] = &createProblem<RNA>;
}

shared_ptr<Problem> ProblemFactory::create(string problem_name) {
	transform(problem_name.begin(), problem_name.end(), problem_name.begin(),
             [](unsigned char c) { return toupper(c); });
	if (problem_map.find(problem_name) == problem_map.end()) {
		stringstream message;
		message<<__FILE__ <<" " <<__LINE__<<" Error: Unable to create handler for problem " << problem_name<<endl; 
		throw logic_error(message.str().c_str()); 
	} else {
		shared_ptr<Problem> my_ptr(problem_map[problem_name]());
		return my_ptr;
	}
	
} 