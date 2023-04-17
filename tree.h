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

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

/**
 *  Clade
 *
 *  This class represents one clade in a phylogeny
 */
class Clade{
	friend class Newick;
	
	std::vector<Clade*> _children;
	std::vector<double> _branch_lengths;
	std::string         _name;
	Clade *             _parent;
	std::vector<int>    _taxon_indices;
	
  public:
    /**
	 *  Visitor
	 *
	 *  This abstract class is used by traverse(...) to apply an operation to all noded in phylogeny.
	 */
	class Visitor{
	  public:
		virtual bool visit(Clade *, Clade*)=0;
	};
	
	std::vector<Clade*> get_children() {return _children;}
	
	std::vector<double> get_branch_lengths() {return _branch_lengths;}
	
	std::string         get_name() {return _name;}
	
	void                set_name(std::string name) { _name=name;}
	
	Clade *             get_parent() {return _parent;}
	
	/**
	 *  traverse
	 *
	 *  Used to apply an operation to all noded in phylogeny
	 */
	bool                traverse(Visitor* visitor, Clade * parent=NULL);
	
	void add_taxon(int seq) {
		_taxon_indices.push_back(seq);
		std::cout <<__FILE__ << " " << __LINE__ << ":"<<_taxon_indices.size()<< " " << _taxon_indices.back()<< std::endl;
	}
	
	std::vector<int>    get_taxon_indices() {return _taxon_indices;}
	
	bool is_leaf() {return _children.size()==0;}
	
	/**
	 *  Destructor deletes children
	 */
	virtual ~Clade() {
		for (std::vector<Clade*>::iterator child=_children.begin(); child!=_children.end(); child++)
			delete *child;
	}
};

/**
 * Taxa
 *
 *  This class represent the set of names for leaves in one or more trees
 */
class Taxa  {
	friend class ConsistencyChecker;
	
	/**
	 * Names of leaves
	 */
	std::vector<std::string>  _names;
	
	/**
	 * Determine position of a specified name in list of all names
	 */

	std::map<std::string,int> _positions;
	
  public:
    /**
    * Constructor: create lookup tables from a list of taxon names
    */
	Taxa(std::string taxa_string);
	
	/**
	 *  Accessor to determine sequence number of specifed taxon
	 */
	int get_position(std::string name) {return _positions[name];}
	
	std::set<int> get_complement(std::set<int> leaves);
};

/**
 *  Edge
 *
 *  And edge joins two Clade in phylogeny
 */
 
class Edge : std::tuple<Clade*,Clade*>{
  public:
	Edge(Clade* first,Clade* second) : std::tuple<Clade*,Clade*>(first,second) {;}
};

/**
 *  CladeNamer
 * 
 * Assigned names to internal nodes
 */
class CladeNamer : public Clade::Visitor{
	int seq                   = 0;
	std::map<std::string, Clade*> _clades;
  public:
	bool visit(Clade *, Clade *);
	Clade* get_clade(std::string name) {return _clades[name];};
};

/**
 *  EdgeBuilder
 *
 *  Constructs set of all edges in phylogeny
 */
class EdgeBuilder : public Clade::Visitor{
	std::vector<Edge> & _edges;
	Taxa &              _taxa;
  public:
	EdgeBuilder(std::vector<Edge> &edges, Taxa & taxa) : _edges(edges), _taxa(taxa) {};
	bool visit(Clade *, Clade*);
}; 


/**
 * ConsistencyChecker
 *
 * Used to verify that the phylogeny uses the exact same set of names as the specified Taxa object
 */
class ConsistencyChecker : Clade::Visitor{
	std::map<std::string,int> _counts;
	
  public:
	bool visit(Clade *, Clade*);
	bool is_consistent(Clade *,Taxa& taxa);
};


/**
 * TreeParser
 *
 * Used to parse a tree from some specified string
 */
class TreeParser {
  public:
	virtual Clade * parse(std::string s)=0;
};
	
/**
 * Tree
 *
 * Encapsulate a hierarchu of clades and a consistent collection of taxa
 */
class Tree {
	friend class TreeFactory;
	Clade*            _root;
	std::vector<Edge> _edges;
	Taxa&             _taxa;
	CladeNamer        _clade_namer;
 	Tree(Clade * root, Taxa & taxa);
	
  public:
	std::vector<Edge> get_edges() {return _edges;};
	bool              traverse(Clade::Visitor* visitor) {return _root->traverse(visitor);};
	virtual ~Tree() {delete _root;};
};

/**
 * TreeFactory
 *
 * Create a tree
 */
class TreeFactory {
	TreeParser& _parser;
	Taxa &      _taxa;
  public:
	TreeFactory(TreeParser& parser, Taxa &taxa): _parser(parser),_taxa(taxa){};
	Tree * create(std::string s) {return new Tree(_parser.parse(s),_taxa);};
};

/**
 *  DescendentVisitor
 *
 *  Find the leaves that are descended from each node
 */
class DescendentVisitor : public Clade::Visitor{
	std::map<std::string,std::set<int>> _descendents;
	Taxa&                               _taxa;
  public:
	DescendentVisitor(Taxa & taxa) : _taxa(taxa) {};
	bool visit(Clade *, Clade*);
	std::set<int> get_descendents(std::string name) {return _descendents[name];};
  private:
	void _propagate(Clade * clade,std::set<int> descendents);
}; 

#endif
