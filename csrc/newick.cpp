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
 
 Node::Node(const int depth, const string name, const double distance) :  _id(count++),_name(name),_depth(depth),_distance(distance) {}
 
 ostream& operator<<(ostream& os, const Node& node){
	os  << "Node " << node._id 	<< " ["<<node._name << "]: "
		<< " depth = " << node._depth
		<< ", distance =  " << node._distance;
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

	// for (auto bounds : bounds_list){
		// int a,b;
		// tie(a,b) = bounds;
		// cout << __FILE__ <<" " <<__LINE__ << ": "<< tokens[a]<<", "  <<tokens[b]<<endl;
		// auto  child =create_node(tokens,a,b,depth+1);
		// node->append(child);
	// }

	return node;
 }
 
 tuple<shared_ptr<Node>,vector<tuple<int,int>>> Newick::explore(vector<Token> tokens,
																const int from,
																const int to, 
																const int depth){
	cout <<__FILE__ <<" " <<__LINE__ << " from="<< from <<", to=" << to << ", depth =" << depth << endl;												
	int working_depth = depth;
	int open = tokens.size();
	int close = -1;
	vector<tuple<int,int>> bounds;
	string name = "";
	auto distance = 1.0;
	auto expect_distance = false;
	auto terminated = false;
	vector <int> splitting_points;
	for (auto i = from;i < to; i++){
		cout <<__FILE__ <<" " <<__LINE__ << " " <<i << tokens[i] <<endl;
		 switch(tokens[i].get_type()) {
			case Token::Type::L:
				if (depth ==working_depth)
					open = i;
				working_depth++;
				break;
				
			case Token::Type::Comma:
				if (working_depth == depth + 1)
					splitting_points.push_back(i+1);
				break;
			case Token::Type::R:
				working_depth--;
				if (depth == working_depth)
					close = i;
				break;
			case Token::Type::Semicolon:
				if (working_depth != 0 || to - i != 1)
					throw _create_error(tokens[i],working_depth,__FILE__,__LINE__);
				terminated = true;
				break;
			case Token::Type::Space:
				break;
			case Token::Type::Colon:
				if (working_depth > depth) break;
				expect_distance = true;
				break;
			case Token::Type::Identifier:
				if (working_depth > depth) break;
				name = tokens[i].get_text();
				break;
			case Token::Type::Number:
				if (working_depth > depth) break;
				if (expect_distance)
					distance = tokens[i].get_numeric();
				else
					throw _create_error(tokens[i],working_depth,__FILE__,__LINE__);
				break;
	}};
	// if (!terminated)
		// throw _create_error(tokens.back(),working_depth,__FILE__,__LINE__);
	const shared_ptr<Node> node = make_shared<Node>(depth,name,distance);
	
	splitting_points.push_back(close);
	for (auto j = 0;j < (int)splitting_points.size();j++)
		if (j==0)
			bounds.push_back(make_tuple(open,splitting_points[j]));
		else
			bounds.push_back(make_tuple(splitting_points[j-1],splitting_points[j]));
		
	return make_tuple(node,bounds);
 }

logic_error Newick::_create_error(Token token, const int depth,const string file,const int line){
	stringstream message;
	message<<file <<" " <<line << ": "<<" Unexpected token " << token<<endl; 
	return logic_error(message.str().c_str()); 
} 
 
 // switch(tokens[i].get_type()) {
			// case Token::Type::L:
				// working_depth++;
				// break;
				
			// case Token::Type::Comma:
				// break;
			// case Token::Type::R:
				// working_depth--;
				// break;
			// case Token::Type::Semicolon:
				// break;
			// case Token::Type::Space:
				// break;
			// case Token::Type::Colon:
				// break;
			// case Token::Type::Identifier:
				// break;
			// case Token::Type::Number:
				// break;
	// };