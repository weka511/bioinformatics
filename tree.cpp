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
 
 #include "tree.h"
 
 bool Clade::traverse(Visitor* visitor, Clade * parent) {
	for (std::vector<Clade*>::iterator child=_children.begin(); child!=_children.end(); child++)
		if (!(*child)->traverse(visitor, this))
			return false;
	return visitor->visit(this,parent);
}

Taxa::Taxa(std::string taxa_string){
	std::istringstream taxa_stream(taxa_string);
	while (taxa_stream) {
		std::string taxon_name;
		taxa_stream >> taxon_name;
		if (taxon_name.size()>0) {
			_positions[taxon_name] = _names.size();
			_names.push_back(taxon_name);
		}
	}
}

std::set<int> Taxa::get_complement(std::set<int> leaves) {
	std::set<int> result;
	for (int i=0;i<_names.size();i++)
		if (leaves.find(i)==leaves.end())
			result.insert(i);
	return result;
}
	
bool ConsistencyChecker::visit(Clade * clade, Clade * parent) {
	if (clade->get_name().size()==0) return true;
	if (_counts.count(clade->get_name())==1)
		_counts[clade->get_name()]++;
	else
		std::cout <<__FILE__ << " " << __LINE__ << ":"<< "Key \"" << clade->get_name() << "\" in tree is not present in list of taxa" << std::endl;
	return true;
}

bool ConsistencyChecker::is_consistent(Clade * root,Taxa& taxa){
	_counts.clear();
	for (std::vector<std::string>::iterator name=taxa._names.begin(); name!=taxa._names.end(); name++)
		_counts[*name] =0;
	bool result = root->traverse(this);
	for (std::map<std::string, int>::iterator it = _counts.begin(); it != _counts.end(); it++)
		if (it->second!=1){
			std::cout <<__FILE__ << " " << __LINE__ << ":" << "Taxon \"" << it->first << "\" has count " << it->second << " - should be 1." << std::endl;
			result = false;
		}
	return result;
}

Tree::Tree(Clade * root, Taxa & taxa) : _root(root), _taxa(taxa) {
	_root->traverse(&_clade_namer);
	EdgeBuilder edgeBuilder(_edges,taxa);
	_root->traverse(&edgeBuilder); 
	DescendentVisitor descendentVisitor(taxa);
	_root->traverse(&descendentVisitor);
	_descendents = descendentVisitor._descendents;
}



bool CladeNamer::visit(Clade * clade, Clade * parent) {
	if ( clade->get_name().size()==0)
		clade->set_name( std::to_string(seq++));
	_clades[clade->get_name()] = clade;
	return true;
}

/**
 *  EdgeBuilder::visit(...)
 *
 *  Used to construct set of all edges in phylogeny
 */
bool EdgeBuilder::visit(Clade * clade, Clade * parent) {
	if (parent!=NULL){
		Edge edge(parent,clade);
		_edges.push_back(edge);
	} 

	return true;
}

bool DescendentVisitor::visit(Clade * clade, Clade* parent) {
	if (clade->is_leaf()) {
		int leaf_index = _taxa.get_position(clade->get_name());
		_descendents[clade->get_name()].insert(leaf_index);
	}
	while (parent!=NULL){
		_propagate(parent,_descendents[clade->get_name()]);
		parent = parent->get_parent();
	}

	return true;
}

void DescendentVisitor::_propagate(Clade * clade,std::set<int> descendents){
	for (std::set<int>::iterator descendent=descendents.begin();descendent!=descendents.end();descendent++)
		_descendents[clade->get_name()].insert(*descendent);
}
