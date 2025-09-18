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
class Parser {
  public:
  	class NewickNode {
	  private:
		vector<shared_ptr<Parser::NewickNode>> _children;
		
	  public:
		virtual void parse(span<Token> tokens) = 0;
		
		shared_ptr<Parser::NewickNode> get_child(int index) {return _children[index];}
		
		virtual string get_name() {return "";}
		
		virtual double get_length() {return 1.0;}
		
	  protected:
	    logic_error create_error(Token token, const string file,const int line);
		
		void attach(shared_ptr<Parser::NewickNode> node);
		
	  	/**
		 * Used by get_first_comma_at_top_level to indicate that it failed to find comma
		 */
		const int UNDEFINED = -1;
	};
	
	class Tree : public NewickNode {
	  public:
		void parse(span<Token> tokens);
	};
	
	class SubTree : public NewickNode {
	  public:
		void parse(span<Token> tokens);
	};
	
	class Leaf : public NewickNode {
	  public:
		void parse(span<Token> tokens);
	};
	
	class Internal : public NewickNode {
	  public:
		void parse(span<Token> tokens);
	};
	
	class BranchSet : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		
		int get_first_comma_at_top_level(span<Token> tokens);
	};
	
	class Branch : public NewickNode {
	  public:
		void parse(span<Token> tokens);
	};
	
	class Name : public NewickNode {
	  private:
		string _value;
	  public:
		void parse(span<Token> tokens);
		
		string get_name() {return _value;};
	};
	
	class Length : public NewickNode {
	  private:
		double _value = 1.0;
	  public:
		void parse(span<Token> tokens);
		
		double get_length() {return _value;}
	};
	
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
	
	/**
	 *   Used by parse_branchset() to locate the first comma that is relevant to
	 *   the parse.
	 *
	 *   Returns:
	 *       index of comma, or UNDEFINED
	 */
	int get_first_comma_at_top_level(span<Token> tokens);
	
	/**
	 * Used by get_first_comma_at_top_level to indicate that it failed to find comma
	 */
	const int UNDEFINED = -1;
	
  private:
	/**
	 * Indicates that parser encountered an error
	 */
	logic_error _create_error(Token token, const string file,const int line);

};



#endif // _NEWICK_HPP
