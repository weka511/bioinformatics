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
#include "catch.hpp"
#include "newick.h"

TEST_CASE( "Newick tests", "[newick]" ) {
	
	SECTION("empty tree"){
		Newick newick;
		Clade * root = newick.parse("();");
		REQUIRE (root->children.size()==0);
	} 
	
  	SECTION("no nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(,,(,));");
		REQUIRE (root->children.size()==3);
		Clade * child_0 = root->children[0];
		REQUIRE(child_0->children.size()==0);
		Clade * child_1 = root->children[1];
		REQUIRE(child_1->children.size()==0);
		Clade * child_2 = root->children[2];
		REQUIRE(child_2->children.size()==2);
	}
	
/* 	SECTION("leaf nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(A,B,(C,D));");
		REQUIRE (root->clades.size()==3);
		Clade * c1 = root->clades[0];
		std::cout<<c1->name<<std::endl;
//		REQUIRE (c1->name=='A');
	//	Clade * c2 = root->clades[2];
	//	REQUIRE (c2->clades.size()==2);
	} */
	
/* 		SECTION("all nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(A,B,(C,D)E)F;");
		REQUIRE (root->clades.size()==3);
	//	Clade * c2 = root->clades[2];
	//	REQUIRE (c2->clades.size()==2);
	} */
}