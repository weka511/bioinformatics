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

using namespace std;

/**
 *  This class parses a string to a Newick tree
 *
 */
class Parser {
  public:
	class Visitor;
	
	/*
	 * This class is used to parse Newick format test. It is the parent
	 * for classes that represent the various elements in the definition
	 * of Newick syntax.	 
	 * Tree -> Subtree ";"
	 * Subtree -> Leaf | Internal
	 * Leaf -> Name
	 * Internal -> "(" BranchSet ")" Name
	 * BranchSet -> Branch | Branch "," BranchSet
	 * Branch -> Subtree Length
	 * Name -> empty | string
	 * Length -> empty | ":" number 
	 */
  	class NewickNode {
	  private:
		/**
	     * Nodes that are directly attached to this one.
		 */
		vector<shared_ptr<Parser::NewickNode>> _children;
		
		/**
		 * Identifies role of node in syntac tree
		 */
		string _type;
		
	  public:
	   /**
	    * Used to convert a tree, specified in Newick format, to an actual
		* tree structure.
	    */
		virtual void parse(span<Token> tokens) = 0;
		
		/**
		 * Used to retrieve a particular child of this node
		 */
		shared_ptr<Parser::NewickNode> get_child(int index) {
			return _children[index];
		}
		
		/**
		 * Name (if any) from Newick string
		 */
		virtual string get_name() const {return "-";}
		
		/**
		 * Distance (if any) from Newick string
		 */
		virtual double get_length() const {return 1.0;}
		
		/**
		 *  Used in conjunction with Visitor to perform an operation on every node.
		 *  Traversal is by recurseive descent.
		 */
		void descend(shared_ptr<Visitor> visitor, const int depth=0);
		
		/**
		 *   Used to output Node
		 */
		friend ostream& operator<<(ostream& os, const NewickNode& node);
		
		/**
		 *   Used to output Node
		 */
		virtual string get_str() const {
			return _type;
		}
		
	  protected:
		/**
		 * Used by descendent types to initialize type variable.
		 */
		NewickNode(string type) : _type(type) {};
		
		/**
		 *  Factory method to create exceptions
		 */
		 
	    logic_error create_error(Token token, const string file,const int line) const;
		
		/**
		 * Used to make a node a child of this one.
		 */
		void attach(shared_ptr<Parser::NewickNode> node);
		
		/**
		 * Used to determine type of node
		 */
		string get_type() const {return _type;}
		
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
		Tree() : NewickNode("Tree") {};
		
		/**
		 * Try to perform production:
		 *
		 * Tree -> Subtree ";"
		 */
		void parse(span<Token> tokens);
	};
	
	/**
	 * Subtree -> Leaf | Internal
	 */
	class SubTree : public NewickNode {
		
	  public:
		SubTree() : NewickNode("SubTree") {};
		
		/**
		 * Try to perform production:
		 *
		 * Subtree -> Leaf | Internal
		 */
		void parse(span<Token> tokens);

	};
	
	/**
	 * Leaf -> Name
	 */	
	class Leaf : public NewickNode {
	  public:
	  
		Leaf() : NewickNode("Leaf") {};
		
		/**
		 * Try to perform production:
		 *
		 * Leaf -> Name
		 */
		void parse(span<Token> tokens);

	};
	
	/**
	 * Internal -> "(" BranchSet ")" Name
	 */	
	class Internal : public NewickNode {
	  public:
		Internal() : NewickNode("Internal") {};
		
		/**
		 * Try to perform production:
		 *
		 * Internal -> "(" BranchSet ")" Name
		 */
		void parse(span<Token> tokens);

	};
	
	/**
	 * BranchSet -> Branch | Branch "," BranchSet
	 */
	class BranchSet : public NewickNode {
	  public:
		BranchSet() : NewickNode("BranchSet") {};
		
		/**
		 * Try to perform production:
		 *
		 * BranchSet -> Branch | Branch "," BranchSet
		 */
		void parse(span<Token> tokens);
		
		/**
		 * Used to separate Branch from rest of BranchSet
		 */		 
		int get_first_comma_at_top_level(span<Token> tokens);
	};
	
	/**
	 * Branch -> Subtree Length
	 */
	class Branch : public NewickNode {
	  public:
		Branch() : NewickNode("Branch") {};
		
		/**
		 * Try to perform production:
		 *
		 * Branch -> Subtree Length
		 */
		void parse(span<Token> tokens);

	};
	
	/**
	 * Name -> empty | string
	 */
	class Name : public NewickNode {
	  private:
		string _value;
		
	  public:
		Name() : NewickNode("Name") {};
		
		/**
		 * Try to perform production:
		 *
		 * Name -> empty | string
		 */
		void parse(span<Token> tokens);
		
		string get_str() const;
		
		string get_name() const{return _value;};
	};
	
	/**
	 * Length -> empty | ":" number  
	 */
	class Length : public NewickNode {
	  private:
		double _value = 1.0;
	  public:
		Length() : NewickNode("Length") {};
		
		/**
		 * Try to perform production:
		 *
		 * Length -> empty | ":" number  
		 */
		void parse(span<Token> tokens);
		
		string get_str() const;
		
		double get_length() const {return _value;}
	};
	
	/**
	 *  Convert a string to a vector of tokens,
	 *  then attempt to parse to a Tree.
	 */
	shared_ptr<Tree> parse(string source);
	
	/**
	 *  Used in conjunction with NewickNode::descend() to perform an operation on every node.
	 *  Traversal is by recurseive descent.
	 */
	class Visitor {
	  public:
	    /**
		 * This function is invoked for each node by NewickNode::descend()
		 */
		virtual void accept(NewickNode * node, const int depth) = 0;
	};
};



#endif // _NEWICK_HPP
