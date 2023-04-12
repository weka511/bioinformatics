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
	enum Token { Left,  Comma, Right, Name, Colon, Length, Semicolon, Invalid };

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
	
 	Tokenizer(std::string s):_s(s) {
		std::cout << _s << std::endl;
	}

	std::string GetString(size_t pos, size_t len){
		return _s.substr(pos,len);
	}
};

class Clade{
  public:
	Clade() : name(""),parent(NULL) {}
	
	std::vector<Clade*> children;
	std::string         name;
	Clade *             parent;
};

class Newick {
	
	std::vector<Clade*> _stack;
	Clade *             _root;
	public:
		Clade * parse(std::string s);
		
	private:
		void _parseLeft();
		void _parseComma();
		void _parseRight();
		void _parseName(std::string s);
		void _parseColon();
		void _parseLength(std::string s);
		void _parseSemicolon();
};



#endif
