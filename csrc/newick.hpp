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
 
#ifndef _NEWICK_HPP
#define _NEWICK_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "tokenizer.hpp"

using namespace std;

/**
 *  This class represents one node in a tree
 */
class Node{
	
  private:
	/**
	 * The children of this Node
	 */
	vector<shared_ptr<Node>> _children;
	
	/**
	 * The identifier is unique within a tree
	 */
	const int _id;
	
	/**
	 *  Name will be filled in from string representation
	 */
	string _name ="";
	
	/**
	 *   Depth of node in tree. This will be set when tree is built.
	 */
	const int _depth;
	
  public:
	static int count;
	
	Node(const int depth);
	
	void append(shared_ptr<Node> node) {_children.push_back(node);}; 
	
	friend ostream& operator<<(ostream& os, const Node & node);
	
	int get_depth() {return _depth;}; 
};
	
class Newick {
  public:

	void parse(string s);
	
	tuple<shared_ptr<Node>,vector<tuple<int,int>>> explore(vector<Token> tokens, const int from, const int to, const int depth);

	shared_ptr<Node> create_node(vector<Token> tokens,
									const int from,
									const int to,
									const int depth	);
  private:
	logic_error _create_error(const Token token, const int depth);
};



#endif // _NEWICK_HPP
