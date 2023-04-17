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
#include "tree.h"
#include "newick.h"


TEST_CASE( "Consistency tests", "[consist]" ) {	
 	SECTION("Consistency"){
		Newick             newick;
		Clade *            root = newick.parse("(A,C,((B,D),E));");
		Taxa               taxa("A B C D E");
		ConsistencyChecker consistency_checker;
		REQUIRE(consistency_checker.is_consistent(root,taxa));
		delete root;
	}
	
	SECTION("Inconsistency"){
		Newick             newick;
		Clade *            root = newick.parse("(A,C,((B,D),F));");
		Taxa               taxa("A B C D E");
		ConsistencyChecker consistency_checker;
		REQUIRE_FALSE(consistency_checker.is_consistent(root,taxa));
		delete root;
	}
	
	SECTION("Edges"){
		Newick      newick;
		Taxa        taxa("A B C D E");
		TreeFactory factory(newick,taxa);
		Tree *      tree = factory.create("(A,C,((B,D),E));");
		REQUIRE (tree->get_edges().size()==7);
		delete tree;
	}
	
 	SECTION("Descendents"){
		Newick      newick;
		Taxa        taxa("A B C D E");
		TreeFactory factory(newick,taxa);
		Tree *      tree = factory.create("(A,C,((B,D),E));");
		std::set<int> d0 = tree->get_descendents("0");
		REQUIRE(d0.size()==2);
		REQUIRE(d0.find(1)!=d0.end());
		REQUIRE(d0.find(3)!=d0.end());
		std::set<int> d1 = tree->get_descendents("1");
		REQUIRE(d1.size()==3);
		REQUIRE(d1.find(1)!=d1.end());
		REQUIRE(d1.find(3)!=d1.end());
		REQUIRE(d1.find(4)!=d1.end());
		std::set<int> d2 = tree->get_descendents("2");
		REQUIRE(d2.size()==5);
		std::set<int> d1c = taxa.get_complement(d1);
		REQUIRE(d1c.size()==2);
		REQUIRE(d1c.find(0)!=d1.end());
		REQUIRE(d1c.find(5)!=d1.end());
		delete tree;
	} 
	
}
