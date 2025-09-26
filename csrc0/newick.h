/**
 * Copyright (C) 2023 Simon Crase: simon@greenweavez.nz
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
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 */
 
#ifndef _NEWICK_H
#define _NEWICK_H
#include <vector>
#include <string>

/**
 *  Extract tokens from input string
 */
class Tokenizer {
	std::string  _s;
	
  public:
	enum Token {
			Left, 
			Comma, 
			Right,
			Name,
			Colon,
			Length,
			Semicolon,
			Invalid
		};

	/**
	 *  Iterator
	 *
	 * Allows parser to iterate through Tokens in inout string
	 */
	class Iterator {
		int        _position;
		Token      _current_token;
		Tokenizer& _tokenizer;
		int        _length;
		
	  public:
		Iterator(Tokenizer& tokenizer): _tokenizer(tokenizer),_position(0){;}
		
		void First(){
			_position=0;
			_update_token();
		}
		
		void Next(){
			_position += _length;
			_update_token();
		}
		
		bool isDone() {return _position>=_tokenizer._s.length();}
		
		Token CurrentToken(){return _current_token;}
		
		int GetPosition() {return _position;}
		
		int GetLength() {return _length;}
		
	  private:
		/**
		 * _update_token
		 *
		 * Used to recognize the next token as we move through input sting.
		 */
	  	void _update_token();
		
		/**
		  *  _is_alpha
		  *
		  *  Used to recognize a string of letters, which may include underscore (but not as first charcter)
		  */
		bool _is_alpha(char c, bool extended=false);
		/**
		 *  _is_digit
		 *
		 *  Used to recognize a string of digits, which may include a decomal point (but not as first charcter)
		 */	
		bool _is_digit(char c, bool extended=false);
	};
	
 	Tokenizer(std::string s):_s(s) {}

	std::string GetString(size_t pos, size_t len){
		return _s.substr(pos,len);
	}
};



/**
 *   Parser for trees in Newick format
 */
class Newick : public TreeParser {
	
	std::vector<Clade*> _stack;
	double              _default_branch_length;
	
  public:
  
	Newick(double default_branch_length=1.0) : _default_branch_length(default_branch_length) {}
  
  /**
   *   Parse string in Newick format into a tree
   */
	Clade * parse(std::string s);
		
  private:
	/**
      *   Used when parsing a Newick formatted string when a left parenthesis is encountered
	  */
	void    _parseLeft(Tokenizer::Token previous_token);

	  /**
	  *   Used when parsing a Newick formatted string when a comma is encountered
	  */	
	void    _parseComma(Tokenizer::Token previous_token);
	
	/**
	  *   Used when parsing a Newick formatted string when a right parenthesis is encountered
	  */
	void    _parseRight(Tokenizer::Token previous_token);
	
	/**
	  *   Used when parsing a Newick formatted string when a name is encountered
	  */
	
	void    _parseName(std::string s, Tokenizer::Token previous_token);
	
	/**
	  *   Used when parsing a Newick formatted string when a colon is encountered
	  */	
	void    _parseColon(Tokenizer::Token previous_token);
	
	/**
	  *   Used when parsing a Newick formatted string when a lenght is encountered
	  */
	void    _parseLength(std::string s,Tokenizer::Token previous_token);

	/**
	  *   Used when parsing a Newick formatted string when a semicolon is encountered
	  */	
	void    _parseSemicolon(Tokenizer::Token previous_token);
	
	/**
	 *    Locate the current top of stack
	 */
	Clade * _get_top_of_stack(){return _stack.back();}
	
	/**
	 *     Find depth of the current top of stack
	 */
	int     _get_depth() {return _stack.size();}
	
	/**
	 *    Get the latest clade to be added
	 */
	Clade * _get_latest() {return  _get_top_of_stack()->_children.back();}
	 
	/**
	 *   Establish a parent-child relationship between the current top of the stack and a newly created clade
	 */
	void    _link(Clade * newly_created){
											newly_created->_parent =_get_top_of_stack();
											_get_top_of_stack()->_children.push_back(newly_created);
										}
};

#endif // _NEWICK_H
