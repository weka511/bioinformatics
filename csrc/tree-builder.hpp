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

#include "node.hpp"
#include "newick.hpp"

using namespace std;

class TreeBuilder : public Parser::Visitor { 

  private:
	/**
	 * This holds a tree that we ae building
	 */
	shared_ptr<Node> _result;
	
	/**
	 *  We keep the most recent nodes, which are still in progress, on the stack.
	 */
	vector<shared_ptr<Node>> _stack;
	
  public:	
  
	/**
	 *  Factory method to parse a tree, in Newick format, into a tree of Nodes.
	 *
	 *  Parameters:
	 *		newick_string            A atring in Newick format
	 *      force_semicolon_at_end   See  #166 Incompatibility between reading Newick file and hand coded Newich string
	 */
	static shared_ptr<Node> create(string newick_string,bool force_semicolon_at_end=false);
	
	/**
	 * Start processing current node, which is either Internal or a Leaf
	 */
	void accept(Parser::NewickNode* node);
	
	/**
	 * End of current node
	 */
	void farewell(Parser::NewickNode * node);
	
	/**
	 * Get the tree that we have built
	 */
	shared_ptr<Node> get_result() {
		return _result;
	}
 };

#endif // _TREE_BUILDER_HPP

