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
 * You should   received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 */
 
#ifndef _ROSALIND_HPP
#define _ROSALIND_HPP

#include <string>


using namespace std;

class Parameters{
	string _problem_name = "????";
	
  public:
	Parameters(int argc, char **argv);
	
	string get_problem_name() {return _problem_name;}
};

#endif // _ROSALIND_HPP
