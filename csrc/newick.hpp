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
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "token.hpp"
#include "node.hpp"

using namespace std;





/**
 *  This class parses a string to a Newick tree
 *
 * Tree -> Subtree ";"
 * Subtree -> Leaf | Internal
 * Leaf -> Name
 * Internal -> "(" BranchSet ")" Name
 * BranchSet -> Branch | Branch "," BranchSet
 * Branch -> Subtree Length
 * Name -> empty | string
 * Length -> empty | ":" number 
 */
class Parser{
  public:
  
	shared_ptr<Node> parse(string source);
	
	/**
	 * Tree -> Subtree ";"
	 */
	shared_ptr<Node> parse_tree(span<Token> tokens);
	
	/**
	 * Subtree -> Leaf | Internal
	 */
	shared_ptr<Node> parse_subtree(span<Token> tokens);

	/**
	 * Leaf -> Name
	 */		
	shared_ptr<Node> parse_leaf(span<Token> tokens);
	
	/**
	 * Internal -> "(" BranchSet ")" Name
	 */	
	shared_ptr<Node> parse_internal(span<Token> tokens);
	
	/**
	 * BranchSet -> Branch | Branch "," BranchSet
	 */
	shared_ptr<Node> parse_branchset(span<Token> tokens);
	
	/**
	 * Branch -> Subtree Length
	 */
	shared_ptr<Node> parse_branch(span<Token> tokens);
	
	/**
	 * Name -> empty | string
	 */
	shared_ptr<Node> parse_name(span<Token> tokens);
	
	/**
	 * Length -> empty | ":" number  
	 */
	double parse_length(span<Token> tokens);
	
	int get_first_comma_at_top_level(span<Token> tokens);
	
  private:
	logic_error _create_error(Token token, const string file,const int line);

};




#endif // _NEWICK_HPP
