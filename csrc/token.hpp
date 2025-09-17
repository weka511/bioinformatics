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
 
#ifndef _TOKEN_HPP
#define _TOKEN_HPP

#include <iostream>
#include <string> 
#include <map>

using namespace std;

/**
 * This class is used when parsing a string (e.g. by Newick). 
 * Each Token represents one lexical elemant from string
 */
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
	  *   Create a token from a string of text.
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
	
	/**
	 *  Find numeric value if token is a floating point number
	 */
	double get_numeric();

	/**
	 * Position in original string, for reporting errors
	 */
	int get_position() {
		return _position;
	};
	
	/**
	 *  Used to output token
	 */
	friend ostream& operator<<(ostream& os, const Token& token);
};

#endif // _TOKEN_HPP