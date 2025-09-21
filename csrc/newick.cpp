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
 
#include <string> 
#include <sstream>

#include "newick.hpp"
#include "tokenizer.hpp"


 using namespace std;

/**
 * Indicates that parser encountered an error
 */
logic_error Parser::NewickNode::create_error(Token token, const string file,const int line) const {
	stringstream message;
	message<<file <<" " <<line << ": "<<" Unexpected token " << token<<endl; 
	return logic_error(message.str().c_str()); 
} 

/**
 * Used to make a node a child of this one.
 */
void Parser::NewickNode::attach(shared_ptr<Parser::NewickNode> node){
	_children.push_back(node);
}

/**
 *  Used in conjunction with Visitor to perform an operation on every node.
 *  Traversal is by recurseive descent.
 */
void Parser::NewickNode::descend(shared_ptr<Visitor> visitor, const int depth){
	visitor->accept(this,depth);
	for (auto child : _children)
		child->descend(visitor,depth+1);
}

/**
 * Allows us to output node
 */
 ostream& operator<<(ostream& os, const Parser::NewickNode& node){
	os  << node.get_str() ;
    return os;
}

/**
 * Tree -> Subtree ";"
 */  
void Parser::Tree::parse(span<Token> tokens){
	span<Token> all_but_last(tokens.begin(), tokens.end() - 1); 
	Token token = tokens.back();
	if (token.get_type() == Token::Type::Semicolon){
		shared_ptr<Parser::SubTree> sub_tree = make_shared<Parser::SubTree>();
		sub_tree->parse(all_but_last);
		attach(sub_tree);
	} else
		throw create_error(token, __FILE__, __LINE__);
}

/**
 * Subtree -> Leaf | Internal
 */

void Parser::SubTree::parse(span<Token> tokens){
	try {
		shared_ptr<Parser::Leaf> leaf = make_shared<Parser::Leaf>();
		leaf->parse(tokens);
		attach(leaf);
	} catch (const logic_error& e){
		try {
			shared_ptr<Parser::Internal> internal = make_shared<Parser::Internal>();
			internal->parse(tokens);
			attach(internal);
		} catch (const logic_error& e1) {
			throw e1;
		}
	}		
}


/**
 * Leaf -> Name
 */
void Parser::Leaf::parse(span<Token> tokens){
	if (tokens.size() == 1 || tokens.size() == 3){
		shared_ptr<Parser::Name> name = make_shared<Parser::Name>();
		name->parse(tokens);
		attach(name);
	}else
		throw  create_error(tokens[0], __FILE__, __LINE__); 
}


/**
  * Internal -> "(" BranchSet ")" Name
  */	
void Parser::Internal::parse(span<Token> tokens){
	if ((int)tokens.size() > 2 && tokens[0].get_type() == Token::Type::L && tokens.back().get_type() == Token::Type::R){
		span<Token> candidate_branchset(tokens.begin()+1, tokens.end() -1); 
		shared_ptr<Parser::BranchSet> branch_set = make_shared<Parser::BranchSet>();
		branch_set->parse(candidate_branchset);
		attach(branch_set);
	} else
		throw  create_error(tokens[0], __FILE__, __LINE__); 
}

/**
 * BranchSet -> Branch | Branch "," BranchSet
 */
void Parser::BranchSet::parse(span<Token> tokens){
	const auto first_comma_at_top_level = get_first_comma_at_top_level(tokens);
	shared_ptr<Parser::Branch> branch = make_shared<Parser::Branch>();
	if (first_comma_at_top_level == UNDEFINED){
		branch->parse(tokens);
		attach(branch);
	} else {
		span<Token> head(tokens.begin(), tokens.begin() + first_comma_at_top_level); 
		branch->parse(head);
		attach(branch);
		span<Token> tail(tokens.begin() + first_comma_at_top_level +1, tokens.end()); 
		shared_ptr<Parser::BranchSet> branch_set = make_shared<Parser::BranchSet>();
		branch_set->parse(tail);
		attach(branch_set);
	}
}

/**
 * Used to separate Branch from rest of BranchSet
 */		
int Parser::BranchSet::get_first_comma_at_top_level(span<Token> tokens){
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
			default: ;
	};
	return UNDEFINED;
}

/**
 * Branch -> Subtree Length
 */
void Parser::Branch::parse(span<Token> tokens){//TODO length
	shared_ptr<Parser::SubTree> sub_tree = make_shared<Parser::SubTree>();
	sub_tree->parse(tokens);
	attach(sub_tree);
}


/**
 * Name -> empty | string
 */
void Parser::Name::parse(span<Token> tokens){
	_value = tokens[0].get_text();
	if (tokens.size() == 1) return;
	span<Token> tail(tokens.begin() + 1, tokens.end());
	shared_ptr<Parser::Length> length = make_shared<Parser::Length>() ;
	length->parse(tail);
	attach(length);
}

string Parser::Name::get_str() const {
	stringstream str;
	str << get_type() << " " << get_name();
	return str.str();
}
								
/** 
 * Length -> empty | ":" number  
 */
void Parser::Length::parse(span<Token> tokens){
	if (tokens[0].get_type() ==Token::Type::Colon)
		_value = tokens[1].get_numeric();
	else
		throw create_error(tokens[0],__FILE__,__LINE__);
}

string Parser::Length::get_str() const  {
	string repr = to_string(get_length());
	while (repr.length() > 3 && repr[repr.length()-1] == '0' && repr[repr.length()-2] == '0')
		repr.pop_back();
	stringstream str;
	str << get_type() << " " << repr;
	return str.str();
}
	
/**
 *  Convert a string to a vector of tokens,
 *  then attempt to parse to a Tree.
 */
shared_ptr<Parser::Tree> Parser::parse(string source){
	Tokenizer tokenizer;
	shared_ptr<Parser::Tree> tree = make_shared<Parser::Tree>();
	auto tokens = tokenizer.tokenize(source);
	tree->parse(tokens);	
	return tree;
}
