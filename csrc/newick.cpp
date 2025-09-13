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

#include <sstream>

#include "newick.hpp"


 
 using namespace std;
 
 void Newick::parse(string s){
	Tokenizer tokenizer;
	auto tokens = tokenizer.tokenize(s);
	explore(tokens,0,tokens.size());
	for (auto token : tokens){
		if (token.is_separator()){}
	}
 }
 
 tuple<shared_ptr<Node>,vector<tuple<int,int>>> Newick::explore(vector<Token> tokens, const int from, const int to){
	 int depth = -1;
	 shared_ptr<Node> node =make_shared<Node>();
	 vector<tuple<int,int>> bounds;
	 for (auto i = from;i < to; i++)
		 switch(tokens[i].get_type()) {
			 
			case Token::Type::Undefined:
				throw _create_error(tokens[i],depth);
			case Token::Type::L:
				depth++;
				break;
			case Token::Type::Comma:
				break;
			case Token::Type::R:
				depth--;
				if (depth < 0) throw _create_error(tokens[i],depth);
				break;
			case Token::Type::Semi:
				break;
			case Token::Type::Space:
				break;
			case Token::Type::Colon:
				break;
			case Token::Type::Identifier:
				break;
			case Token::Type::Number:
				break;
		};
	return make_tuple(node,bounds);
 }

logic_error Newick::_create_error(const Token token, const int depth){
	stringstream message;
	message<<__FILE__ <<" " <<__LINE__<<" Unexpected " <<endl; 
	return logic_error(message.str().c_str()); 
} 
 /**
 switch(tokens[i].get_type()) {
			case Token::Type::Undefined:
			case Token::Type::L:
			case Token::Type::Comma:
			case Token::Type::R:
			case Token::Type::Semi:
			case Token::Type::Space:
			case Token::Type::Colon:
			case Token::Type::Identifier:
			case Token::Type::Number:
			break;
		};
 */