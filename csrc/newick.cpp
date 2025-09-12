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
#include <stack>
#include "newick.hpp"
#include "tokenizer.hpp"

 
 using namespace std;
 
 void Newick::parse(string s){
	Tokenizer tokenizer;
	auto tokens = tokenizer.tokenize(s);
	stack<Node> nodes;
	for (auto token : tokens){
		if (token.is_separator()){}
	}
 }
 
 