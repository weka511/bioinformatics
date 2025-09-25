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
 
#ifndef _NODE_HPP
#define _NODE_HPP

#include <memory>
#include <string>
#include <vector>
#include <iostream>


using namespace std;

/**
 *  This class represents a node in a tree
 */
class Node{
	
  private:
    /**
	 * Used to assign a unique id to each node
	 */
  	static int count;
	
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
	 *  Distance of this node from parent
	 */
	double _distance = 1.0;

  public:
    /**
	 *  This class is used to traverse a tree
	 */
	class Visitor {
	  public:
	   /**
		*  When Node::visit() is called, it executes `accept' for each node.
		*/
		virtual void accept(Node& node,int depth)=0;
	};

	/**
	 * Used in testing to reset Node::count to zero
	 */
	static void reset();
	
	/**
	 *  Create and intialize a node
	 */
	Node();
	
	/**
	 *   Add a child node
	 */
	void append(shared_ptr<Node> node) {
		_children.push_back(node);
	}; 
	
	/**
	 * When the bulider gets a name or a length, it is always for the latest child
	 */
	shared_ptr<Node> get_last_child() const {
		return _children.back();
	}
	
	/**
	 * Allows us to output node
	 */
	friend ostream& operator<<(ostream& os, const Node & node);
	
	/**
	 *  Distance of this node from parent
	 */
	double get_distance() const {
		return _distance;
	}; 
	
	/**
	 *   Assign distance to node
	 */
	void set_distance (const double distance) {
		_distance = distance;
	}
	
	/**
	 *  Name will be filled in from string representation
	 */
	string get_name() const {
		return _name;
	}; 
	
	/**
	 *   Assign name to node
	 */
	void set_name(const string name) {
		_name = name;
	}
	
	/**
	 *  This method is used to traverse a tree,
	 *  invoking the Visitor's `accept` method for each node
	 */
	void visit(Visitor & visitor,int depth=0);
};

#endif // _NODE_HPP
