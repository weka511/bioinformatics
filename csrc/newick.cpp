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
 
 #include <iostream>
#include <string> 
#include <cstddef>  
 #include "newick.hpp"
 
 using namespace std;
 
 vector<string> Tokenizer::tokenize(string str) {
	 cout << str << endl;
	 vector<string> tokens;
	 long unsigned int start = 0; // Avoid warning when I compare with `found`
	 auto found = str.find_first_of(_separators);
	 if (found > start)
		 tokens.push_back(str.substr(start,found-start));
	 
	 // `found` points to a token
	 while (found != string::npos) {
		 tokens.push_back(str.substr(found,1));
		 start = found + 1;
		 found = str.find_first_of(_separators,start);
		 if (found > start)
			tokens.push_back(str.substr(start,found-start));
	 }

	 return tokens;
 }