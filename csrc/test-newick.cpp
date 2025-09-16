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
 
 class Displayer: public Node::Visitor {
	 void accept(Node& node,int depth){
		 cout << __FILE__ << " " << __LINE__  <<node << " " << depth << endl;
	 }
 };
 
TEST_CASE( "Newick Tests", "[newick]" ) {
	Parser parser;
	Node::reset();
	Tokenizer tokenizer;
	Displayer displayer;
	
	SECTION("get_first_comma_at_top_level") {
		auto tokens = tokenizer.tokenize("A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5");	
		REQUIRE(parser.get_first_comma_at_top_level(tokens) == 3);
	}
	
	SECTION("get_first_comma_at_top_level(1)") {
		auto tokens = tokenizer.tokenize("(C:0.3,D:0.4)E:0.5,A:0.1,B:0.2");	
		REQUIRE(parser.get_first_comma_at_top_level(tokens) == 12);
	}

	
	SECTION("a naked leaf node ") {
		auto tokens = tokenizer.tokenize("SomeLeaf;");	
		shared_ptr<Node> node = parser.parse_tree(tokens);
		REQUIRE(node->get_name() == "SomeLeaf");
	}
	
	
	SECTION("Test parser with labels") {
		auto tokens = tokenizer.tokenize("(A,B,(C,D));");
		shared_ptr<Node> node = parser.parse_tree(tokens);
		node->visit(displayer);
	}
	
	SECTION("Test parser with labels and lengths") {
		auto tokens = tokenizer.tokenize("(A:2.1,B,(C,D));");
		shared_ptr<Node> node = parser.parse_tree(tokens);
		node->visit(displayer);
	}
	
	SECTION("Test parser with labels and lengths(1)") {
		shared_ptr<Node> node = parser.parse("(A:2.1,B:1.2,(C:3.0,D:0.2));");
		node->visit(displayer);
	}
	
 }