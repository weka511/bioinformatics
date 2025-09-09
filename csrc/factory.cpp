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

#include "factory.hpp"

using namespace std;

shared_ptr<Problem> ProblemFactory::create(string problem_name) {
	if (problem_name == "DNA") 
		return make_shared<DNA>();
	 else {
		stringstream message;
		message<<__FILE__ <<" " <<__LINE__<<" Error: Unable to create problem " << problem_name<<endl; 
		throw logic_error(message.str().c_str()); 
	}
} 