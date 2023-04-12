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
 
#include <iostream>
#include <cstddef> 
#include <string>
#include <cassert>
#include "newick.h"

void  Tokenizer::Iterator::_update_token(){
	_length = 1;
	if (_tokenizer._s[_position]=='(')
		_current_token=Left;
	else if (_tokenizer._s[_position]==',')
		_current_token=Comma;
	else if (_tokenizer._s[_position]==')')
		_current_token=Right;
	else if (_tokenizer._s[_position]==':')
		_current_token=Colon;
	else if (_is_alpha(_tokenizer._s[_position])) {
		_current_token = Name;
		_length        = 0;
		while (_is_alpha(_tokenizer._s[_position+_length],true))
			_length++;
	} else if (_is_digit(_tokenizer._s[_position])) {
		_current_token = Length;
		while (_is_digit(_tokenizer._s[_position+_length],true))
			_length++;
	} else if (_tokenizer._s[_position]==';')
		_current_token=Semicolon;
	else 
		_current_token = Invalid;
}

bool Tokenizer::Iterator::_is_alpha(char c, bool extended){
	bool in_range= 'A'<=c and c<='Z';
	if (extended)
		return in_range || c=='_';
	else
		return in_range;
}
	
bool Tokenizer::Iterator::_is_digit(char c, bool extended){
	bool in_range= '0'<=c and c<='9';
	if (extended)
		return in_range || c=='.';
	else
		return in_range;
}

Clade * Newick::parse(std::string s){
	assert(_stack.size()==0);
	Tokenizer tokenizer(s);
	Tokenizer::Iterator iterator(tokenizer);
	iterator.First();
	while (!iterator.isDone()){
		Tokenizer::Token token = iterator.CurrentToken();
		switch(token) {
			case Tokenizer::Token::Left:
				_parseLeft();
				break;
			case Tokenizer::Token::Comma:
				_parseComma();
				break;
			case Tokenizer::Token::Right:
				_parseRight();
				break;
			case Tokenizer::Token::Name:
				_parseName(tokenizer.GetString(iterator.GetPosition(), iterator.GetLength()));
				break;
			case Tokenizer::Token::Colon:
				_parseColon();
				break;
			case Tokenizer::Token::Length:
				_parseLength(tokenizer.GetString(iterator.GetPosition(), iterator.GetLength()));
				break;
			case Tokenizer::Token::Semicolon:
				iterator.Next();
				assert(iterator.isDone());
				break;
			default:
				std::cout <<"WTF "<<token<< std::endl;
		}
		iterator.Next();
	}
	return _root;
}

void Newick::_parseColon(){
}

void Newick::_parseLeft(){
	Clade * clade = new Clade();
	if (_stack.size()>0)
		_stack.back()->clades.push_back(clade);
	else
		_root = clade;
	_stack.push_back(clade);
}

void Newick::_parseComma(){
	assert(_stack.size()>0);
	Clade * clade = new Clade();
	_stack.back()->clades.push_back(clade);
}

void Newick::_parseRight(){
	assert(_stack.size()>0);
	_stack.pop_back();
}

void Newick::_parseName(std::string s){
	
}

void Newick::_parseLength(std::string s){
	
}



