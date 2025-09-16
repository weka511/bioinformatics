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
#include "tokenizer.hpp"

 using namespace std;
 

shared_ptr<Node> Parser::parse(string source){
	Tokenizer tokenizer;
	auto tokens = tokenizer.tokenize(source);	
	return parse_tree(tokens);
}

/**
 * Tree -> Subtree ";"
 */
shared_ptr<Node> Parser::parse_tree(span<Token> tokens) {
	// cout << __FILE__ << " " << __LINE__  << " parse_tree"<<endl;
	// for (auto t : tokens)
		// cout << t << endl;
	span<Token> all_but_last(tokens.begin(), tokens.end() -1); 
	Token token = tokens.back();
	if (token.get_type() == Token::Type::Semicolon)
		return parse_subtree(all_but_last);
	else
		throw _create_error(token, __FILE__, __LINE__);
}

/**
 * Subtree -> Leaf | Internal
 */
shared_ptr<Node> Parser::parse_subtree(span<Token> tokens){
	// cout << __FILE__ << " " << __LINE__  << " parse_subtree"<< endl;
	// for (auto t : tokens)
		// cout << t << endl;
	try {
		return parse_leaf(tokens);
	} catch (const logic_error& e){
		try {
			return parse_internal(tokens);
		} catch (const logic_error& e1) {
			throw e1;
		}
	}		
}

/**
 * Leaf -> Name
 */
shared_ptr<Node> Parser::parse_leaf(span<Token> tokens){
	// cout << __FILE__ << " " << __LINE__  << " parse_leaf "<< tokens.size() <<endl;
	// for (auto t : tokens)
		// cout << t << endl;
	if (tokens.size() == 1 || tokens.size() == 3)
		return parse_name(tokens);
	else
		throw  _create_error(tokens[0], __FILE__, __LINE__); 
}

/**
  * Internal -> "(" BranchSet ")" Name
  */	
shared_ptr<Node> Parser::parse_internal(span<Token> tokens){
	// cout << __FILE__ << " " << __LINE__  << " parse_internal"<< endl;
	// for (auto t : tokens)
		// cout << t << endl;
	if ((int)tokens.size() > 2 &&
		tokens[0].get_type() == Token::Type::L &&
		tokens.back().get_type() == Token::Type::R){
		span<Token> candidate_branchset(tokens.begin()+1, tokens.end() -1); 
		return parse_branchset(candidate_branchset);
	} else
		throw  _create_error(tokens[0], __FILE__, __LINE__);  
}


/**
 * BranchSet -> Branch | Branch "," BranchSet
 */
shared_ptr<Node> Parser::parse_branchset(span<Token> tokens){
	// cout << __FILE__ << " " << __LINE__ << " parse_branchset" << endl;
	// for (auto t : tokens)
		// cout << t << endl;
	const auto first_comma_at_top_level = get_first_comma_at_top_level(tokens);
	if (first_comma_at_top_level < 0)
		return parse_branch(tokens);
	else {
		// cout << __FILE__ << " " << __LINE__ <<  ",=" << first_comma_at_top_level <<endl;
		span<Token> head(tokens.begin(), tokens.begin() + first_comma_at_top_level); 
		parse_branch(head);
		span<Token> tail(tokens.begin() + first_comma_at_top_level +1, tokens.end()); 
		parse_branchset(tail);
	}
	const shared_ptr<Node> node =make_shared<Node>("",0);
	return node;
}

/**
 * Branch -> Subtree Length
 */
shared_ptr<Node> Parser::parse_branch(span<Token> tokens){
	// cout << __FILE__ << " " << __LINE__  << " parse_branch"<< endl;
	// for (auto t : tokens)
		// cout << t << endl;
	parse_subtree(tokens);  // TODO Length
	const shared_ptr<Node> node =make_shared<Node>("",0);
	return node;
}

/**
 * Name -> empty | string
 */
shared_ptr<Node> Parser::parse_name(span<Token> tokens){
	cerr << __FILE__ << " " << __LINE__  << " parse_name"<< endl;
	for (auto t : tokens)
		cerr << t << endl;
	if (tokens.size() == 1){
		const shared_ptr<Node> node =make_shared<Node>(tokens[0].get_text());
		cout << __FILE__ << " " << __LINE__ << " " << *node << endl;
		return node;
	} else {
		span<Token> tail(tokens.begin() + 1, tokens.end()); 
		auto distance = parse_length(tail);
		const shared_ptr<Node> node =make_shared<Node>(tokens[0].get_text(),distance);
		cout << __FILE__ << " " << __LINE__ << " " << *node << endl;
		return node;
	}
}

/**
 * Length -> empty | ":" number  
 */
double Parser::parse_length(span<Token> tokens){
	// cout << __FILE__ << " " << __LINE__  << " parse_length"<< endl;
	if (tokens[0].get_type() ==Token::Type::Colon)
		return tokens[1].get_numeric();
	else
		throw _create_error(tokens[0],__FILE__,__LINE__);

}

logic_error Parser::_create_error(Token token, const string file,const int line){
	stringstream message;
	message<<file <<" " <<line << ": "<<" Unexpected token " << token<<endl; 
	return logic_error(message.str().c_str()); 
} 
	
int Parser::get_first_comma_at_top_level(span<Token> tokens){
	auto depth = 0;
	for (int i=0; i<(int)tokens.size();i++)
		switch(tokens[i].get_type()) {
			case Token::Type::L:
				depth++;
				break;
			case Token::Type::Comma:
				if (depth == 0)
					return i;
				break;
			case Token::Type::R:
				depth--;
				break;
			case Token::Type::Semicolon:
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
	return -1;
}