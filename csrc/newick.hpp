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
	class Visitor;
	
  	class NewickNode {
	  private:
		vector<shared_ptr<Parser::NewickNode>> _children;
		
	  public:
		virtual void parse(span<Token> tokens) = 0;
		
		shared_ptr<Parser::NewickNode> get_child(int index) {return _children[index];}
		
		virtual string get_name() {return "";}
		
		virtual double get_length() {return 1.0;}
		
		void descend(shared_ptr<Visitor> visitor);
		
		friend ostream& operator<<(ostream& os, const NewickNode& node);
		
		virtual string get_str() const =0;
	  protected:
	    logic_error create_error(Token token, const string file,const int line);
		
		void attach(shared_ptr<Parser::NewickNode> node);
		
	  	/**
		 * Used by get_first_comma_at_top_level to indicate that it failed to find comma
		 */
		const int UNDEFINED = -1;
	};
	
	/**
	 * Tree -> Subtree ";"
	 */
	class Tree : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return "Tree";}
	};
	
	/**
	 * Subtree -> Leaf | Internal
	 */
	class SubTree : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return "SubTree";}
	};
	
	/**
	 * Leaf -> Name
	 */	
	class Leaf : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return "Leaf";}
	};
	
	/**
	 * Internal -> "(" BranchSet ")" Name
	 */	
	class Internal : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return "Internal";}
	};
	
	/**
	 * BranchSet -> Branch | Branch "," BranchSet
	 */
	class BranchSet : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		
		string get_str() const  {return "BranchSet";}
		
		int get_first_comma_at_top_level(span<Token> tokens);
	};
	
	/**
	 * Branch -> Subtree Length
	 */
	class Branch : public NewickNode {
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return "Branch";}
	};
	
	/**
	 * Name -> empty | string
	 */
	class Name : public NewickNode {
	  private:
		string _value;
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return _value;}
		string get_name() {return _value;};
	};
	
	/**
	 * Length -> empty | ":" number  
	 */
	class Length : public NewickNode {
	  private:
		double _value = 1.0;
	  public:
		void parse(span<Token> tokens);
		string get_str() const  {return to_string(_value);}
		double get_length() {return _value;}
	};
	
	shared_ptr<Tree> parse(string source);
	
	class Visitor {
	  public:
		virtual void accept(NewickNode * node) = 0;
	};
};



#endif // _NEWICK_HPP
