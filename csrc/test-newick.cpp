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
  
 
TEST_CASE( "Newick Tests", "[newick]" ) {
	Newick newick;
	Node::count = 0;
	Tokenizer tokenizer;
	
	SECTION("Test parser with label and top level distance") {
		auto tokens = tokenizer.tokenize("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F:6.7;");	
		REQUIRE(tokens.size() == 26);
		shared_ptr<Node> node;
		vector<tuple<int,int>> bounds;
		tie(node,bounds) = newick.explore(tokens, 0, tokens.size(), 0);
		REQUIRE(node->get_depth() == 0);
		REQUIRE(node->get_distance() == 6.7);
		REQUIRE(node->get_name() == "F");
		for (auto bound : bounds){
			cout <<__FILE__ <<" " <<__LINE__ << ": " << get<0>(bound) << "-" << get<1>(bound) << endl;
			shared_ptr<Node> child;
			vector<tuple<int,int>> bounds1;
			tie(child,bounds1) = newick.explore(tokens, get<0>(bound), get<1>(bound), 1);
			cout <<__FILE__ <<" " <<__LINE__ << ": " << *child << endl;
		}

	}
	
	// SECTION("Test parser with label without top level distance") {
		// auto tokens = tokenizer.tokenize("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;");	
		// shared_ptr<Node> node;
		// vector<tuple<int,int>> bounds;
		// tie(node,bounds) = newick.explore(tokens, 0, tokens.size(), 0);
		// REQUIRE(node->get_depth() == 0);
		// REQUIRE(node->get_distance() == 1.0);
		// REQUIRE(node->get_name() == "F");
		// for (auto bound : bounds)
			// cout << get<0>(bound) << "-" << get<1>(bound) << endl;
	// }
	
	// SECTION("a tree rooted on a leaf node (rare)") {
		// auto tokens = tokenizer.tokenize("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;");	
		// shared_ptr<Node> node;
		// vector<tuple<int,int>> bounds;
		// tie(node,bounds) = newick.explore(tokens, 0, tokens.size(), 0);
		// REQUIRE(node->get_depth() == 0);
		// REQUIRE(node->get_distance() == 1.0);
		// REQUIRE(node->get_name() == "A");
		// for (auto bound : bounds)
			// cout << get<0>(bound) << "-" << get<1>(bound) << endl;
	// }
	
	// SECTION("a naked leaf node ") {
		// auto tokens = tokenizer.tokenize("B:0.2,");	
		// shared_ptr<Node> node;
		// vector<tuple<int,int>> bounds;
		// tie(node,bounds) = newick.explore(tokens, 0, tokens.size(), 1);
		// REQUIRE(node->get_depth() == 0);
		// REQUIRE(node->get_distance() == 0.2);
		// REQUIRE(node->get_name() == "B");
	// }
	
	// SECTION("Test parser with no labels") {
		// Newick newick;
		// newick.parse("(,,(,));");
	// }
	
	// SECTION("Test parser with labels") {
		// Newick newick;
		// newick.parse("(A,B,(C,D));");
	// }
 }