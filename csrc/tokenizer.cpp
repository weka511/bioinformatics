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


#include "tokenizer.hpp"
 
 using namespace std;
 
 
	
 /**
  *  Convert a string to a vector of tokens
  */
 vector<Token> Tokenizer::tokenize(string str,bool force_semicolon_at_end) {
	vector<Token> product;
	long unsigned int right = str.length();
	long unsigned int start = 0; // Avoid warning when I compare with `next_token`
	auto next_token = str.find_first_of(_separators);
	if (next_token > start)
		 product.push_back(Token(str.substr(start,next_token-start),false,start,_type_map));
	 
	while (next_token != string::npos) {
		 product.push_back(Token(str.substr(next_token,1),true,next_token,_type_map));
		 start = next_token + 1;
		 next_token = str.find_first_of(_separators,start);
		 if (next_token > start && start < right)
			product.push_back(Token(str.substr(start,next_token-start),false,start,_type_map));
	}

	return _cull_consecutive_spaces(product,force_semicolon_at_end);
 }
 
 /**
  * Replace consecutive white spaces with a single space
  */
 vector<Token> Tokenizer::_cull_consecutive_spaces(vector<Token> tokens,bool force_semicolon_at_end){
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
	/**
	 *  Force semicolon to be last token in vector.
	 */
	if (force_semicolon_at_end)
		while (product.back().get_type() != Token::Type::Semicolon)
			product.pop_back();
	return product;
 }