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
		// case Parser::Type::Tree: 
			// return;
		// case Parser::Type::SubTree: 
			// return;
		case Parser::Type::Leaf: 
			if (_needs_comma.back())
				_result << ",";
			else
				_needs_comma.back() = true;
			return;
		case Parser::Type::Internal:
			_result << "(" ;
			_needs_comma.push_back(false);
			return;
		case Parser::Type::BranchSet: 
			return;
		// case Parser::Type::Branch: 
			// return;
		case Parser::Type::Name: 
			_result << node->get_name();
			return;
		// case Parser::Type::Length: 
			// return;
		default:
			return;
	}
}

void TreeBuilder::farewell(Parser::NewickNode * node, const int depth){
		switch (node->get_type()){
		// case Parser::Type::Tree: 
			// return;
		// case Parser::Type::SubTree: 
			// return;
		// case Parser::Type::Leaf: 
			// return;
		case Parser::Type::Internal:
			_needs_comma.pop_back();
			_result << ")";
			return;
		// case Parser::Type::BranchSet: 
			// return;
		// case Parser::Type::Branch: 
			// return;
		// case Parser::Type::Name: 
			// return;
		// case Parser::Type::Length: 
			// return;
		default:
			return;
	}
}

string TreeBuilder::get_result(){
	return _result.str();
}