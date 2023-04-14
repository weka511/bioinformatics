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
 

#include "catch.hpp"
#include "newick.h"

TEST_CASE( "Newick tests", "[newick]" ) {
	
  	SECTION("Trivial trees"){
		Newick newick;
		Clade * tree0 = newick.parse("();");
		REQUIRE (tree0->get_children().size()==1);
		delete tree0;
		Clade * tree5 = newick.parse("(,,,,);");
		REQUIRE (tree5->get_children().size()==5);
		delete tree5;
	} 
	
    SECTION("no nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(,,(,));");
		REQUIRE (root->get_children().size()==3);
		Clade * child_0 = root->get_children()[0];
	 	REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		delete root;
	} 
	
  	SECTION("leaf nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(Alice,Bob,(Charlie,David));");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_children()[0]->get_name()=="Alice");
		REQUIRE (root->get_children()[1]->get_name()=="Bob");
		Clade * child_0 = root->get_children()[0];
	 	REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		REQUIRE (child_2->get_children()[0]->get_name()=="Charlie");
		REQUIRE (child_2->get_children()[1]->get_name()=="David");
		delete root;
	}  
	
	SECTION("leaf nodes are named - 3 levels"){
		Newick newick;
		Clade * root = newick.parse("(A,C,((B,D),E));");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_children()[0]->get_name()=="A");   // 1st level -- (A,C,...)
		REQUIRE (root->get_children()[1]->get_name()=="C");
		Clade * child_2 = root->get_children()[2];            // 2nd level -- (...,E)
		REQUIRE(child_2->get_children().size()==2);
		Clade * child_20 = child_2->get_children()[0];
		REQUIRE(child_20->get_children().size()==2);
		REQUIRE (child_20->get_children()[0]->get_name()=="B");
		REQUIRE (child_20->get_children()[1]->get_name()=="D");
		Clade * child_21 = child_2->get_children()[1];
		REQUIRE (child_21->get_name()=="E");
		delete root;
	}  
	
 	SECTION("all nodes are named"){
		Newick newick;
		Clade * root = newick.parse("(Alice,Bob,(Charlie,David)Ermintrude)Frank;");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_children()[0]->get_name()=="Alice");
		REQUIRE (root->get_children()[1]->get_name()=="Bob");
		REQUIRE (root->get_children()[2]->get_name()=="Ermintrude");
		Clade * child_0 = root->get_children()[0];
		REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		REQUIRE (child_2->get_children()[0]->get_name()=="Charlie");
		REQUIRE (child_2->get_children()[1]->get_name()=="David");
		REQUIRE (root->get_name()=="Frank");
		delete root;
	} 
	
    SECTION("Branch lengths"){
		Newick newick;
		Clade * root = newick.parse("(:0.1,:0.2,(:0.3,:0.4):0.5);");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_branch_lengths()[0]==0.1);
		REQUIRE (root->get_branch_lengths()[1]==0.2);
		REQUIRE (root->get_branch_lengths()[2]==0.5);
		Clade * child_0 = root->get_children()[0];
	 	REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		REQUIRE (child_2->get_branch_lengths()[0]==0.3);
		REQUIRE (child_2->get_branch_lengths()[1]==0.4);
		delete root;
	} 
	
	SECTION("Branch lengths - one missing"){
		Newick newick(42);
		Clade * root = newick.parse("(:0.1,,(:0.3,:0.4):0.5);");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_branch_lengths()[0]==0.1);
		REQUIRE (root->get_branch_lengths()[1]==42.0);
		REQUIRE (root->get_branch_lengths()[2]==0.5);
		Clade * child_0 = root->get_children()[0];
	 	REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		REQUIRE (child_2->get_branch_lengths()[0]==0.3);
		REQUIRE (child_2->get_branch_lengths()[1]==0.4);
	} 
	
 	SECTION("Branch lengths -- all"){
		Newick newick;
		Clade * root = newick.parse("(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_branch_lengths()[0]==0.1);
		REQUIRE (root->get_branch_lengths()[1]==0.2);
		REQUIRE (root->get_branch_lengths()[2]==0.5);
		Clade * child_0 = root->get_children()[0];
	 	REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		REQUIRE (child_2->get_branch_lengths()[0]==0.3);
		REQUIRE (child_2->get_branch_lengths()[1]==0.4);
		REQUIRE(root->get_parent()->get_branch_lengths()[0]==0);
		delete root;
	} 

 	SECTION("leaf nodes are named, and we have lengths"){
		Newick newick;
		Clade * root = newick.parse("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);");
		REQUIRE (root->get_children().size()==3);
		REQUIRE (root->get_children()[0]->get_name()=="A");
		REQUIRE (root->get_children()[1]->get_name()=="B");
		REQUIRE (root->get_branch_lengths()[0]==0.1);
		REQUIRE (root->get_branch_lengths()[1]==0.2);
		REQUIRE (root->get_branch_lengths()[2]==0.5);
		Clade * child_0 = root->get_children()[0];
	 	REQUIRE(child_0->get_children().size()==0);
		Clade * child_1 = root->get_children()[1];
		REQUIRE(child_1->get_children().size()==0);
		Clade * child_2 = root->get_children()[2];
		REQUIRE(child_2->get_children().size()==2);
		REQUIRE (child_2->get_children()[0]->get_name()=="C");
		REQUIRE (child_2->get_children()[1]->get_name()=="D");
		REQUIRE (child_2->get_branch_lengths()[0]==0.3);
		REQUIRE (child_2->get_branch_lengths()[1]==0.4);
		delete root;
	} 	
}