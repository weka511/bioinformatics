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

using namespace std;



class Token {
  public:
  	/**
	 * Parser needs to know what type of token this is.
	 */
  	enum class Type{
		Undefined,
		L,
		Comma,
		R,
		Semi,
		Space,
		Colon,
		Identifier,
		Number
	};
	
  private:
	/**
	 * Actual text of token.
	 */
	const string _text;
	
	/**
	 *  Whether token matches one of the separators.
	 */
	const bool _is_separator;
	
	/**
	 * Parser needs to know what type of token this is.
	 */
	Type _type;
	
	/**
	 *   Record position of toklen in source line
	 */
	const int _position;
	
  public:

	Token(const string text,const bool is_separator, const int position);
	
	/**
	 * Parser needs to know what type of token this is.
	 */
	Type get_type() {return _type;}
	
	/**
	 * Actual text of token.
	 */
	string get_text() {return _text;};
	
	/**
	 *  Whether token is a space
	 */
	bool is_space() {return _text == " ";};
	
	/**
	 *  Whether token matches one of the separators.
	 */
	bool is_separator() {return _is_separator;}
	
	bool is_numeric() {
		return _type == Type::Number;
	};	
	
	double get_numeric();

	/**
	 *   Position of toklen in source line: for reporting errors
	 */
	int get_position() {return _position;};
};

/**
 *  This class converts a string to a vector of tokens
 */
class Tokenizer {
	string _separators;
	
  public:
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
