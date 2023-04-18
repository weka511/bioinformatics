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

#include <algorithm>
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
	std::set<uint64_t> Q1, Q2;
	_prepare(T1,&Q1);
	_prepare(T2, &Q2);
	std::vector<uint64_t> Q_intersection;
    std::set_intersection(Q1.begin(), Q1.end(), Q2.begin(), Q2.end(),
                          std::back_inserter(Q_intersection));
	return Q1.size() + Q2.size() - 2*Q_intersection.size();
}

void QuartetDistanceCalculator::_prepare(Tree* T,std::set<uint64_t> * quartets){
	for (std::vector<Edge>::iterator edge = T->get_edges_begin(); edge != T->get_edges_end();++edge) {
		if (!edge->is_internal()) continue;
		Clade* first = edge->first;
		Clade* second = edge->second;
		std::set<int> B = T->get_descendents(second->get_name());
		std::set<int> A = _taxa.get_complement(B);
		for (std::set<int>::iterator a1=A.begin();a1!=A.end();++a1)
			for (std::set<int>::iterator a2=A.begin();a2!=A.end();++a2)
				for (std::set<int>::iterator b1=B.begin();b1!=B.end();++b1)
					for (std::set<int>::iterator b2=B.begin();b2!=B.end();++b2)
						if (*a1<*a2 && *b1<*b2)
							quartets->insert(_create_quartet(*a1,*a2,*b1,*b2));		 
	}
}

uint64_t QuartetDistanceCalculator::_create_quartet(int a1,int a2,int b1,int b2) {
 	if (a1>b1)
		return _create_quartet(b1,b2,a1,a2);
	uint64_t result=a1;
	result <<=QuartetDistanceCalculator::shift;
	result += a2;
	result <<=QuartetDistanceCalculator::shift;
	result += b1;
	result <<=QuartetDistanceCalculator::shift;
	result += b2;
	return result;
}