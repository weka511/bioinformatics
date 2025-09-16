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
 
#ifndef _TOKENIZER_HPP
#define _TOKENIZER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <map>

#include "token.hpp"

using namespace std;



/**
 *  This class converts a string to a vector of tokens
 */
class Tokenizer {
	/**
	 * A token is either a separator (from following string) or a 
	 * string of characters not containing a separator.
	 */
	const string _separators;
	
	/**
	 *  Convert a separator to a type FIXME - _separators is redundant
	 */
	const map<string,Token::Type> _type_map={
			{"(", Token::Type::L},
			{",", Token::Type::Comma},
			{")", Token::Type::R},
			{";", Token::Type::Semicolon},
			{" ", Token::Type::Space},
			{":", Token::Type::Colon},
	};	
	
  public:
	/**
	 * Create Tokenizer, and initialize list of separators.
	 */
	Tokenizer(string separators="(,); :") : _separators(separators){};
  
	/**
	 *  Convert a string to a vector of tokens
	 */
	vector<Token> tokenize(string str);
	
  private:
	/**
	 * Replace consecutive white spaces with a single space
	 */
    vector<Token> _cull_consecutive_spaces(vector<Token> tokens);

};

#endif // _TOKENIZER_HPP
