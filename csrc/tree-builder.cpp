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

#include "tree-builder.hpp"

using namespace std;
 
void TreeBuilder::accept(Parser::NewickNode* node,const int depth){
	switch (node->get_type()){
		case Parser::Type::Leaf: 
			_stack.back()->append(make_shared<Node>());
			if (_needs_comma.back())
				_string << ",";
			else
				_needs_comma.back() = true;
			return;
		case Parser::Type::Internal:
			if (_stack.size() ==0) {
				_result = make_shared<Node>();
				_stack.push_back(_result);
				cout << __FILE__<< " " <<__LINE__ << " " <<*_stack.back() << endl;
			} else {
				shared_ptr<Node> new_node = make_shared<Node>();
				_stack.back()->append(new_node);
				_stack.push_back(new_node);
				cout << __FILE__<< " " <<__LINE__ << " " <<*_stack.back() << endl;
			}	
			_string << "(" ;
			_needs_comma.push_back(false);
			return;
		case Parser::Type::Name: 
			_string << node->get_name();
			_stack.back()->get_last_child()->set_name(node->get_name());
			cout << __FILE__<< " " <<__LINE__ << " " <<*_stack.back() << endl;
			return;
		case Parser::Type::Length:
			_string << ":" << node->get_length();
			_stack.back()->get_last_child()->set_distance(node->get_length());
			return;
		default:
			return;
	}
}

void TreeBuilder::farewell(Parser::NewickNode * node, const int depth){
		switch (node->get_type()){
			case Parser::Type::Internal:
				cout << __FILE__<< " " <<__LINE__ << " " <<*(_stack.back()) << endl;
				_stack.pop_back();
				if (_stack.size()==0)
					cout << __FILE__<< " " <<__LINE__ << " " <<*_result << endl;
				_needs_comma.pop_back();
				_string << ")";
				return;
			default:
				return;
	}
}

string TreeBuilder::get_string(){
	return _string.str();
}