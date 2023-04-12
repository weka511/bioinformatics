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
		_length        = 1;
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
	_root = NULL;
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
	assert(_stack.size()==0);
	return _root;
}

void Newick::_parseColon(){
}

void Newick::_parseLeft(){
	if (_stack.size()==0){
		 _root = new Clade();
		 _stack.push_back(_root);
	} else {
		Clade * newly_created = new Clade();
		Clade * top_of_stack  = _stack.back();
		newly_created->parent = top_of_stack;
		top_of_stack->children.push_back(newly_created);
		_stack.push_back(newly_created);
	}

	std::cout << "New level "<<_stack.size() << " and node" << std::endl;
}

void Newick::_parseComma(){
	assert(_stack.size()>0);
	Clade * newly_created = new Clade();
	Clade * top_of_stack  = _stack.back();
	newly_created->parent = top_of_stack->parent;
	std::cout << __LINE__ << std::endl;
	top_of_stack->parent->children.push_back(newly_created);
	std::cout << __LINE__ << std::endl;
	std::cout << "node: stack size= " << _stack.size()
		<< ", siblings=" << top_of_stack->children.size() << std::endl;
}

void Newick::_parseRight(){
	assert(_stack.size()>0);
	_stack.pop_back();
	std::cout << "popped" << std::endl;
}

void Newick::_parseName(std::string s){
	std::cout << s << std::endl;
	std::cout << _stack.size() << std::endl;
	
	_stack.back()->children.back()->name = s;
}

void Newick::_parseLength(std::string s){
	
}




