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

using namespace std;

class Token {
  public:
  	/**
	 * Parser needs to know what type of token this is.
	 */
  	enum class Type{
		L = 0, 
		Comma = 1,
		R = 2,
		Semicolon = 3,
		Space = 4,
		Colon = 5,
		Identifier = 6,
		Number = 7
	};

	/**
	 * Parser needs to know what type of token this is.
	 */
	Type _type;
	
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
	 *   Record position of toklen in source line
	 */
	const int _position;
	
  public:
	/**
	  *   Create a toekn from a string of text.
	  *
	  *   Parameters:
	  *      text           The text that is used to build the token
	  *      is_separator   Indicates whether token is a separator or an ordinry string
	  *      position       Position of token in input
	  */
	Token(const string text,const bool is_separator, const int position,map<string,Token::Type> type_map);
	/*
	 * Parser needs to know what type of token this is.
	 */
	Type get_type() {
		return _type;
	}
	
	/**
	 * Actual text of token.
	 */
	string get_text() {
		return _text;
	};
	
	/**
	 *  Whether token is a space
	 */
	bool is_space() {
		return _text == " ";
	};
	
	/**
	 *  Whether token matches one of the separators.
	 */
	bool is_separator() {
		return _is_separator;
	}
	
	bool is_numeric() {
		return _type == Type::Number;
	};	
	
	double get_numeric();

	/**.: for reporting errors
	 */
	int get_position() {
		return _position;
	};
	
	friend ostream& operator<<(ostream& os, const Token& token);
};

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
