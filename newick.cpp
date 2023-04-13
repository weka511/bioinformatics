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
	bool in_range= ('A'<=c and c<='Z') || ('a'<=c and c<='z');
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
	_stack.push_back(new Clade());     // root
	Tokenizer tokenizer(s);
	Tokenizer::Iterator iterator(tokenizer);
	iterator.First();
	Tokenizer::Token previous_token = Tokenizer::Token::Invalid;
	while (!iterator.isDone()){
		Tokenizer::Token token = iterator.CurrentToken();
		switch(token) {
			case Tokenizer::Token::Left:
				_parseLeft( previous_token);
				break;
			case Tokenizer::Token::Comma:
				_parseComma( previous_token);
				break;
			case Tokenizer::Token::Right:
				_parseRight( previous_token);
				break;
			case Tokenizer::Token::Name:
				_parseName(tokenizer.GetString(iterator.GetPosition(), iterator.GetLength()), previous_token);
				break;
			case Tokenizer::Token::Colon:
				_parseColon( previous_token);
				break;
			case Tokenizer::Token::Length:
				_parseLength(tokenizer.GetString(iterator.GetPosition(), iterator.GetLength()), previous_token);
				break;
			case Tokenizer::Token::Semicolon:
				iterator.Next();
				assert(iterator.isDone());
				break;
			default:
				std::cout <<"WTF "<<token<< std::endl;
		}
		previous_token = token;
		iterator.Next();
	}

	Clade * root = _stack.back();
	_stack.pop_back();
	assert(_stack.size()==0);
	return root->_children.front();
}

void Newick::_parseColon( Tokenizer::Token previous_token){
	assert( previous_token==Tokenizer::Token::Name        || 
			previous_token==Tokenizer::Token::Comma      || 
			previous_token==Tokenizer::Token::Right      || 
			previous_token==Tokenizer::Token::Left);
}

void Newick::_parseLeft( Tokenizer::Token previous_token){
	if (_get_depth()==1){
		Clade * newly_created = new Clade();
		_link(newly_created);
		_stack.push_back(newly_created);
	} else 
		_stack.push_back(_get_latest());
	Clade * first_born = new Clade();
	_get_top_of_stack()->_children.push_back(first_born);
}

void Newick::_parseComma(Tokenizer::Token previous_token){
	assert( previous_token==Tokenizer::Token::Name        || 
			previous_token==Tokenizer::Token::Comma      || 
			previous_token==Tokenizer::Token::Right      || 
			previous_token==Tokenizer::Token::Length      ||
			previous_token==Tokenizer::Token::Left);

	assert(_stack.size()>1);
	
	_link(new Clade());
}

void Newick::_parseRight(Tokenizer::Token previous_token){
	assert( previous_token==Tokenizer::Token::Name        || 
			previous_token==Tokenizer::Token::Comma      || 
			previous_token==Tokenizer::Token::Right      || 
			previous_token==Tokenizer::Token::Length      ||
			previous_token==Tokenizer::Token::Left);
	assert(_stack.size()>0);
	_stack.pop_back();
}

void Newick::_parseName(std::string s, Tokenizer::Token previous_token){
	assert(	previous_token==Tokenizer::Token::Comma      || 
			previous_token==Tokenizer::Token::Right      || 
			previous_token==Tokenizer::Token::Left);
	_get_top_of_stack()->_children.back()->_name = s;
}

void Newick::_parseLength(std::string s, Tokenizer::Token previous_token){
	assert(previous_token==Tokenizer::Token::Colon);
	int deficit = _get_top_of_stack()->_children.size() - _get_top_of_stack()->_branch_lengths.size();
	while (deficit-- >1)
		_get_top_of_stack()->_branch_lengths.push_back(_default_branch_length);
	_get_top_of_stack()->_branch_lengths.push_back(std::stod(s));
}




