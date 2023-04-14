/**
 * Copyright (C) 2023 Simon Crase: simon@greenweavez.nz
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
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 */
 
#ifndef _TREE_H
#define _TREE_H


#include <sstream>
#include <string>
#include <map>
#include <vector>

class Clade{
	friend class Newick;
	
	std::vector<Clade*> _children;
	std::vector<double> _branch_lengths;
	std::string         _name;
	Clade *             _parent;
	
  public:
	class Visitor{
	  public:
		virtual bool visit(Clade *)=0;
	};
	Clade() : _name(""),_parent(NULL) {}
	
	std::vector<Clade*> get_children() {return _children;}
	
	std::vector<double> get_branch_lengths() {return _branch_lengths;}
	
	std::string         get_name() {return _name;}
	
	Clade *             get_parent() {return _parent;}
	
	bool                traverse(Visitor* visitor);
	
	/**
	 *  Destructor deletes children
	 */
	virtual ~Clade() {
		for (std::vector<Clade*>::iterator child=_children.begin(); child!=_children.end(); child++)
			delete *child;
	}
};

class Taxa  {
	friend class ConsistencyChecker;
	std::vector<std::string> _names;
  public:
	Taxa(std::string taxa_string);

};

class ConsistencyChecker : Clade::Visitor{
	std::map<std::string,int> _counts;
	
  public:
	bool visit(Clade *);
	bool is_consistent(Clade *,Taxa& taxa);
};

#endif
