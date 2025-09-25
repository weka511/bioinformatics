/**
 * Copyright (C) 2025 Simon Crase
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
 *
 * (,,(,));                               no nodes are named
 * (A,B,(C,D));                           leaf nodes are named
 * (A,B,(C,D)E)F;                         all nodes are named
 * (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
 * (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
 * (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
 * (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
 * ((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;    a tree rooted on a leaf node (rare)
 */
 
 
#include "catch.hpp"
#include "newick.hpp"
#include "tokenizer.hpp"
 
  
TEST_CASE( "Newick Tests", "[newick]" ) {
	Parser parser;
	Tokenizer tokenizer;
	
	SECTION("get_first_comma_at_top_level") {
		auto tokens = tokenizer.tokenize("A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5");	
		Parser::BranchSet branch_set;
		REQUIRE(branch_set.get_first_comma_at_top_level(tokens) == 3);
	}
	
	SECTION("get_first_comma_at_top_level(1)") {
		auto tokens = tokenizer.tokenize("(C:0.3,D:0.4)E:0.5,A:0.1,B:0.2");	
		Parser::BranchSet branch_set;
		REQUIRE(branch_set.get_first_comma_at_top_level(tokens) == 12);
	}

	SECTION("a naked leaf node ") {
		shared_ptr<Parser::Tree> tree = parser.parse("SomeLeaf;");
		shared_ptr<Parser::NewickNode> subtree = tree->get_child(0);
		shared_ptr<Parser::NewickNode> leaf = subtree->get_child(0);
		shared_ptr<Parser::NewickNode> name = leaf->get_child(0);
		REQUIRE(name->get_name() == "SomeLeaf");
	}
	
	SECTION("a naked leaf node with length") {
		auto tokens = tokenizer.tokenize("SomeLeaf:3.1415;");	
		Parser::Tree tree;
		tree.parse(tokens);
		shared_ptr<Parser::NewickNode> subtree = tree.get_child(0);
		shared_ptr<Parser::NewickNode> leaf = subtree->get_child(0);
		shared_ptr<Parser::NewickNode> name = leaf->get_child(0);
		REQUIRE(name->get_name() == "SomeLeaf");
		shared_ptr<Parser::NewickNode> length = name->get_child(0);
		REQUIRE(length->get_length() == 3.1415);
	}
	
	SECTION("Test parser with labels") {
		auto tree = parser.parse("(A,B,(C,D));");
		auto subtree = tree->get_child(0);
		REQUIRE(subtree->get_type() == Parser::Type::SubTree);
		auto internal = subtree->get_child(0);
		REQUIRE(internal->get_type() == Parser::Type::Internal);
		auto branchset1 = internal->get_child(0);
		REQUIRE(branchset1->get_type() == Parser::Type::BranchSet);
		auto branch1 = branchset1->get_child(0);
		REQUIRE(branch1->get_type() == Parser::Type::Branch);
		auto subtree1 = branch1->get_child(0);
		REQUIRE(subtree1->get_type() == Parser::Type::SubTree);	
		auto leaf1 = subtree1->get_child(0);
		REQUIRE(leaf1->get_type() == Parser::Type::Leaf);
		auto name1 = leaf1->get_child(0);
		REQUIRE(name1->get_type() == Parser::Type::Name);
		REQUIRE(name1->get_name() == "A");
		auto branchset2 = branchset1->get_child(1);
		REQUIRE(branchset2->get_type() == Parser::Type::BranchSet);
		auto branch2 = branchset2->get_child(0);
		REQUIRE(branch2->get_type() == Parser::Type::Branch);
		auto subtree2 = branch2->get_child(0);
		REQUIRE(subtree2->get_type() == Parser::Type::SubTree);	
		auto leaf2 = subtree2->get_child(0);
		REQUIRE(leaf2->get_type() == Parser::Type::Leaf);
		auto name2 = leaf2->get_child(0);
		REQUIRE(name2->get_type() == Parser::Type::Name);
		REQUIRE(name2->get_name() == "B");
		// auto branchset3 = branchset2->get_child(1);
		// REQUIRE(branchset3->get_type() == Parser::Type::BranchSet);
		// auto branch3 = branchset3->get_child(0);
		// REQUIRE(branch3->get_type() == Parser::Type::Branch);
		// auto subtree3 = branch3->get_child(0);
		// REQUIRE(subtree3->get_type() == Parser::Type::SubTree);	
		// auto leaf3 = subtree3->get_child(0);
		// REQUIRE(leaf3->get_type() == Parser::Type::Leaf);
		// auto name3 = leaf3->get_child(0);
		// REQUIRE(name3->get_type() == Parser::Type::Name);
		// REQUIRE(name3->get_name() == "C");
	}
	
	SECTION("Test parser with labels and lengths(1)") {
		auto tree = parser.parse("(A:2.1,B:1.2,(C:3.0,D:0.2));");
		auto subtree = tree->get_child(0);
		REQUIRE(subtree->get_type() == Parser::Type::SubTree);
		auto internal = subtree->get_child(0);
		REQUIRE(internal->get_type() == Parser::Type::Internal);
		auto branchset1 = internal->get_child(0);
		REQUIRE(branchset1->get_type() == Parser::Type::BranchSet);
		auto branch1 = branchset1->get_child(0);
		REQUIRE(branch1->get_type() == Parser::Type::Branch);
		auto subtree1 = branch1->get_child(0);
		REQUIRE(subtree1->get_type() == Parser::Type::SubTree);	
		auto leaf1 = subtree1->get_child(0);
		REQUIRE(leaf1->get_type() == Parser::Type::Leaf);
		auto name1 = leaf1->get_child(0);
		REQUIRE(name1->get_type() == Parser::Type::Name);
		REQUIRE(name1->get_name() == "A");
		// auto distance1 = leaf1->get_child(1);   //segfault!
		// REQUIRE(distance1->get_type() == Parser::Type::Length);
		// REQUIRE(distance1->get_length() == 2.1);
		auto branchset2 = branchset1->get_child(1);
		REQUIRE(branchset2->get_type() == Parser::Type::BranchSet);
		auto branch2 = branchset2->get_child(0);
		REQUIRE(branch2->get_type() == Parser::Type::Branch);
		auto subtree2 = branch2->get_child(0);
		REQUIRE(subtree2->get_type() == Parser::Type::SubTree);	
		auto leaf2 = subtree2->get_child(0);
		REQUIRE(leaf2->get_type() == Parser::Type::Leaf);
		auto name2 = leaf2->get_child(0);
		REQUIRE(name2->get_type() == Parser::Type::Name);
		REQUIRE(name2->get_name() == "B");
	}
	

	
 }