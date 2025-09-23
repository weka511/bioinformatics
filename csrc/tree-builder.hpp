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
 
#ifndef _TREE_BUILDER_HPP
#define _TREE_BUILDER_HPP

#include <vector>
#include <sstream>
#include <string>

#include "node.hpp"
#include "newick.hpp"

using namespace std;

class TreeBuilder : public Parser::Visitor { 

  private:
	vector<bool> _needs_comma;
	stringstream _string;
	shared_ptr<Node> _result;// = make_shared<Node>(); // adding dummy node doesn't fix segmentation fault
	vector<shared_ptr<Node>> _stack;
	
  public:	
	void accept(Parser::NewickNode* node,const int depth);
	
	void farewell(Parser::NewickNode * node, const int depth);
	
	string get_string();
	
	shared_ptr<Node> get_result() {return _result;}
 };

#endif // _TREE_BUILDER_HPP

