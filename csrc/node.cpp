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
 
#include "node.hpp"

 using namespace std;
 
 int Node::count = 0;
 
 /**
  * Used in testing to reset Node::count to zero
  */
 void Node::reset(){
	 Node::count = 0;
 }
 
 /**
  *  Create and intialize a node
  */
 Node::Node( const string name, const double distance) : 
	_id(count++),_name(name),_distance(distance) {}
 
 /**
  *  This method is used to traverse a tree,
  *  invoking the Visitor for each node
  */
 void Node::visit(Visitor & visitor,int depth){
	 visitor.accept(*this,depth);
	 for (auto child : _children)
		 child->visit(visitor,depth+1);
 }
 
/**
 * Allows us to output node
 */
 ostream& operator<<(ostream& os, const Node& node){
	os  << "Node " << node._id 	<< " ["<<node._name << "]: "
		<< ", distance =  " << node._distance;
      return os;
}
