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
		Invalid };

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
	  	void _update_token();
		bool _is_alpha(char c, bool extended=false);
		bool _is_digit(char c, bool extended=false);
	};
	
 	Tokenizer(std::string s):_s(s) {}

	std::string GetString(size_t pos, size_t len){
		return _s.substr(pos,len);
	}
};

class Clade{
	friend class Newick;
	
	std::vector<Clade*> _children;
	std::vector<double> _branch_lengths;
	std::string         _name;
	Clade *             _parent;
	
  public:
	Clade() : _name(""),_parent(NULL) {}
	std::vector<Clade*> get_children() {return _children;}
	std::vector<double> get_branch_lengths() {return _branch_lengths;}
	std::string         get_name() {return _name;}
	Clade *             get_parent() {return _parent;}
};

class Newick {
	std::vector<Clade*> _stack;
	double              _default_branch_length;
  public:
	Newick(double default_branch_length=1.0) : _default_branch_length(default_branch_length) {}
  
	Clade * parse(std::string s);
		
  private:
	void    _parseLeft(Tokenizer::Token previous_token);
	
	void    _parseComma(Tokenizer::Token previous_token);
	
	void    _parseRight(Tokenizer::Token previous_token);
	
	void    _parseName(std::string s, Tokenizer::Token previous_token);
	
	void    _parseColon(Tokenizer::Token previous_token);
	
	void    _parseLength(std::string s,Tokenizer::Token previous_token);
	
	void    _parseSemicolon(Tokenizer::Token previous_token);
	
	Clade * _get_top_of_stack(){return _stack.back();}
	int     _get_depth() {return _stack.size();}
	Clade * _get_latest() {return  _get_top_of_stack()->_children.back();}
	void    _link(Clade * newly_created){
											newly_created->_parent =_get_top_of_stack();
											_get_top_of_stack()->_children.push_back(newly_created);
										}
};



#endif
