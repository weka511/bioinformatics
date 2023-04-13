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
#include <string>

TEST_CASE( "Newick tests", "[newick]" ) {
	
  	SECTION("Trivial trees"){
		Newick newick;
		Clade * tree0 = newick.parse("();");
		REQUIRE (tree0->children.size()==1);
		Clade * tree5 = newick.parse("(,,,,);");
		REQUIRE (tree5->children.size()==5);
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
	
 	SECTION("leaf nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(Alice,Bob,(Charlie,David));");
		REQUIRE (root->children.size()==3);
		REQUIRE (root->children[0]->name=="Alice");
		REQUIRE (root->children[1]->name=="Bob");
		Clade * child_0 = root->children[0];
	 	REQUIRE(child_0->children.size()==0);
		Clade * child_1 = root->children[1];
		REQUIRE(child_1->children.size()==0);
		Clade * child_2 = root->children[2];
		REQUIRE(child_2->children.size()==2);
		REQUIRE (child_2->children[0]->name=="Charlie");
		REQUIRE (child_2->children[1]->name=="David");
	} 
	
 		SECTION("all nodes are named"){
			Newick newick;
			Clade * root = newick.parse("(Alice,Bob,(Charlie,David)Ermintrude)Frank;");
			REQUIRE (root->children.size()==3);
			REQUIRE (root->children[0]->name=="Alice");
			REQUIRE (root->children[1]->name=="Bob");
			REQUIRE (root->children[2]->name=="Ermintrude");
			Clade * child_0 = root->children[0];
			REQUIRE(child_0->children.size()==0);
			Clade * child_1 = root->children[1];
			REQUIRE(child_1->children.size()==0);
			Clade * child_2 = root->children[2];
			REQUIRE(child_2->children.size()==2);
			REQUIRE (child_2->children[0]->name=="Charlie");
			REQUIRE (child_2->children[1]->name=="David");
			REQUIRE (root->name=="Frank");
	} 
}