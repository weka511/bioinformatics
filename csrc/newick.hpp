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
 
#ifndef _NEWICK_HPP
#define _NEWICK_HPP

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Newick {
  public:
	class Node{
	  private:
		vector<Node> nodes;
	};
	 
	void parse(string s);
};

class Token {
	const string _text;
	
  public:
	Token(string text) : _text(text) {};
	
	string get_text() {return _text;};
	
	int get_length() {return _text.length();};
	
	bool is_space() {return _text == " ";};
};

/**
 *  This class converts a string to a vector of tokens
 */
class Tokenizer {
	string _separators;
	
  public:
	Tokenizer(string separators="(,); ") : _separators(separators){};
  
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

#endif // _NEWICK_HPP
