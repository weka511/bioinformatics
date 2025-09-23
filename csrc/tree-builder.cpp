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
 
#include "tree-builder.hpp"

using namespace std;
 
/**
 * Start processing current node, which is either Internal or a Leaf
 */
void TreeBuilder::accept(Parser::NewickNode* node,const int depth){
	switch (node->get_type()){
		case Parser::Type::Leaf: 
			_stack.back()->append(make_shared<Node>());
			return;
			
		case Parser::Type::Internal:
			if (_stack.size() ==0) {
				_result = make_shared<Node>();
				_stack.push_back(_result);
			} else {
				shared_ptr<Node> new_node = make_shared<Node>();
				_stack.back()->append(new_node);
				_stack.push_back(new_node);
			}	
			return;
			
		case Parser::Type::Name: 
			_stack.back()->get_last_child()->set_name(node->get_name());
			return;
			
		case Parser::Type::Length:
			_stack.back()->get_last_child()->set_distance(node->get_length());
			return;
			
		default:
			return;
	}
}

/**
 * End of current node
 */
void TreeBuilder::farewell(Parser::NewickNode * node, const int depth){
		switch (node->get_type()){
			case Parser::Type::Internal:
				_stack.pop_back();
				return;
				
			default:
				return;
	}
}

/**
 *  Factory method to parse a tree, in Newick format, into a tree of Nodes.
 */
 shared_ptr<Node> TreeBuilder::create(string newick_string){
	Parser parser;
	auto tree = parser.parse(newick_string);
	auto builder = make_shared<TreeBuilder>();
	tree->descend(builder);
	return builder->get_result();
 }
