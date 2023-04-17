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
 
#include <iostream>
#include "tree.h"
#include "newick.h"
#include "qrtd.h"



int  QuartetDistanceCalculator::get_distance(std::string T1, std::string T2) {
	Newick      newick;
	TreeFactory factory(newick,_taxa);
	Tree *      tree1 = factory.create(T1);
	Tree *      tree2 = factory.create(T2);
	int         result = _get_distance(tree1,tree2);
	delete tree2;
	delete tree1;
	return result;
}

int QuartetDistanceCalculator::_get_distance(Tree* T1, Tree* T2){
	prepare(T1);
	prepare(T2);
	return -1;
}

void QuartetDistanceCalculator::prepare(Tree* T){
	for (std::vector<Edge>::iterator edge = T->get_edges().begin(); edge != T->get_edges().end();++edge)
		if (edge->is_internal()) {
			Clade* first = edge->first;
			Clade* second = edge->second;
			std::cout <<__FILE__ << " " << __LINE__ <<" "<< first->get_name() << "," << second->get_name() << std::endl;
			std::set<int> B = T->get_descendents(second->get_name());
			std::set<int> A =  _taxa.get_complement(B);
			for (std::set<int>::iterator a=A.begin();a!=A.end();++a)
				for (std::set<int>::iterator b=B.begin();b!=B.end();++b)
					std::cout <<__FILE__ << " " << __LINE__ <<" "<< *a << "->" << *b << std::endl;
	}
}