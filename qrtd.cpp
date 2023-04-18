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
	_prepare(T1);
//	std::set<uint64_t> q2=prepare(T2);
	return -1;
}

void QuartetDistanceCalculator::_prepare(Tree* T){
	std::vector<uint64_t> result;
	for (std::vector<Edge>::iterator edge = T->get_edges_begin(); edge != T->get_edges_end();++edge) {
		std::cout <<__FILE__ << " " << __LINE__  << std::endl;
		if (edge->is_internal()) {
			std::cout <<__FILE__ << " " << __LINE__  << std::endl;
			Clade* first = edge->first;
			Clade* second = edge->second;
			std::cout <<__FILE__ << " " << __LINE__ <<" "<< first->get_name() << "," << second->get_name() << std::endl;
			std::set<int> B = T->get_descendents(second->get_name());
			std::set<int> A = _taxa.get_complement(B);
			std::cout <<__FILE__ << " " << __LINE__ <<" "<<A.size() << "," <<B.size() << std::endl;
			for (std::set<int>::iterator a1=A.begin();a1!=A.end();++a1)
				for (std::set<int>::iterator a2=A.begin();a2!=A.end();++a2)
					for (std::set<int>::iterator b1=B.begin();b1!=B.end();++b1)
						for (std::set<int>::iterator b2=B.begin();b2!=B.end();++b2)
							if (*a1<*a2 && *b1<*b2){
								uint64_t quartet=_create_quartet(*a1,*a2,*b1,*b2);
								std::cout <<__FILE__ << " " << __LINE__ <<" "<<std::hex << quartet <<std::dec << std::endl;
								result.push_back(quartet);
							}	
			std::cout <<__FILE__ << " " << __LINE__ <<" "<<result.size() << std::endl;	
			
			std::cout <<__FILE__ << " " << __LINE__ <<" "<< first->get_name() << "," << second->get_name() << std::endl;
		} else 
			std::cout <<__FILE__ << " " << __LINE__  << std::endl;
	}
	std::cout <<__FILE__ << " " << __LINE__ <<" "<<result.size() << std::endl;	
//	return result;
}

uint64_t QuartetDistanceCalculator::_create_quartet(int a1,int a2,int b1,int b2) {
/* 	if (a1>b1)
		return _create_quartet(b1,b2,a1,a2); */
	uint64_t result=a1;
	result <<=QuartetDistanceCalculator::shift;
	result += a2;
	result <<=QuartetDistanceCalculator::shift;
	result += b1;
	result <<=QuartetDistanceCalculator::shift;
	result += b2;
	return result;
}