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
 
 int Node::count = 0;
 
 Node::Node(const int depth) :  _id(count++),_depth(depth) {}
 
 ostream& operator<<(ostream& os, const Node& node){
      os  << "Node: " << node._id << ",depth = " << node._depth;
      return os;
}

 void Newick::parse(string s){
	 Node::count = 0;
	 cout <<__FILE__ <<" " <<__LINE__ << ": "<< s << endl;
	Tokenizer tokenizer;
	auto tokens = tokenizer.tokenize(s);
	create_node(tokens,0,tokens.size(),0);
 }
 
 shared_ptr<Node> Newick::create_node(vector<Token> tokens,
									const int from,
									const int to,
									const int depth){
	shared_ptr<Node> node;
	vector<tuple<int,int>> bounds_list;
	tie (node,bounds_list) =  explore(tokens,from,to,depth); 
	for (auto bounds : bounds_list){
		int a,b;
		tie(a,b) = bounds;
		cout <<__FILE__ <<" " <<__LINE__ << ": "<<" [" << tokens[a]<<"] "  << "["<<tokens[b]<< "]"<<endl;
		auto  child =create_node(tokens,a,b,depth+1);
		node->append(child);
	}

	return node;
 }
 
 tuple<shared_ptr<Node>,vector<tuple<int,int>>> Newick::explore(vector<Token> tokens,
																const int from,
																const int to,
																const int depth){
	cout <<__FILE__ <<" " <<__LINE__ << " from="<< from <<", to=" << to << ", depth =" << depth << endl;												
	int working_depth = depth;
	shared_ptr<Node> node = make_shared<Node>(depth);
	vector<tuple<int,int>> bounds;
	auto start = 0;
	for (auto i = from;i < to; i++)
		 switch(tokens[i].get_type()) {
			case Token::Type::L:
				working_depth++;
				start = i + 1;
	
				break;
			case Token::Type::Comma:
				if (working_depth == depth + 1) {
					bounds.push_back(make_tuple(start,i));
					start = i + 1;
				}
				break;
			case Token::Type::R:
				if (working_depth == depth + 1){
					bounds.push_back(make_tuple(start,i));
					start = i + 1;
				}
				working_depth--;
				break;
			case Token::Type::Semicolon:
				if (working_depth > 0) _create_error(tokens[i],depth);
				if (i < to) _create_error(tokens[i],depth);
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

logic_error Newick::_create_error(Token token, const int depth){
	stringstream message;
	message<<__FILE__ <<" " <<__LINE__ << ": "<<" Unexpected token " << token<<endl; 
	return logic_error(message.str().c_str()); 
} 
 