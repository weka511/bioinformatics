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
  *   Create a Token from a string of text.
  *
  *   Parameters:
  *      text           The text that is used to build the token
  *      is_separator   Indicates whether token is a separator or an ordinry string
  *      position       POsition of token in input
  */
 Token::Token(const string text,const bool is_separator,const int position,map<string,Token::Type> type_map)
	: _text(text),_is_separator(is_separator),_position(position) {
	if (_is_separator)
		_type = type_map[_text];
	else try{
		stod(_text);
		_type = Type::Number;
	} catch (const invalid_argument& e) {
		_type = Type::Identifier;
	}
}

/**
 *  Find numeric value if token is a floating point number
 */
double Token::get_numeric() {
	return stod(_text);
}	

/**
 *  Used to output token
 */
ostream& operator<<(ostream& os, const Token& token){
      os << (int)token._type << " '" << token._text <<"'"<< " at " << token._position;
      return os;
}