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
 
 void Newick::parse(string s){
	Tokenizer tokenizer;
	auto tokens = tokenizer.tokenize(s);
 }
 
 /**
  *  Convert a string to a vector of tokens
  */
 vector<Token> Tokenizer::tokenize(string str) {
	vector<Token> product;
	long unsigned int right = str.length();
	long unsigned int start = 0; // Avoid warning when I compare with `found`
	auto found = str.find_first_of(_separators);
	if (found > start)
		 product.push_back(Token(str.substr(start,found-start),true));
	 
	// `found` points to a token
	while (found != string::npos) {
		 product.push_back(Token(str.substr(found,1),true));
		 start = found + 1;
		 found = str.find_first_of(_separators,start);
		 if (found > start && start < right)
			product.push_back(Token(str.substr(start,found-start),false));
	}

	return _cull_consecutive_spaces(product);
 }
 
 /**
  * Replace consecutive white spaces with a single space
  */
 vector<Token> Tokenizer::_cull_consecutive_spaces(vector<Token> tokens){
	vector<Token> product;
	auto last_token_was_space = true;
	for (auto token : tokens){
		if (!token.is_space())
			product.push_back(token);
		else 
			if (!last_token_was_space)
				product.push_back(token);
		last_token_was_space = token.is_space();
	}

	return product;
 }