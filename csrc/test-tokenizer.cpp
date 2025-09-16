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
 *((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;    a tree rooted on a leaf node (rare)
 */
 
 #include "catch.hpp"
 #include "tokenizer.hpp"
  
 TEST_CASE( "Tokenizer Tests", "[tokenizer]" ) {
	 
	 SECTION("Test tokenizer newick") {
		Tokenizer tokenizer;
		auto tokens = tokenizer.tokenize("(C,(B,D),(A,E));");
		REQUIRE(tokens[0].get_text() == "(");
		REQUIRE(tokens[0].is_separator());
		REQUIRE(tokens[1].get_text() == "C");
		REQUIRE(!tokens[1].is_separator());
		REQUIRE(tokens[2].get_text() == ",");
		REQUIRE(tokens[2].is_separator());
		REQUIRE(tokens[3].get_text() == "(");
		REQUIRE(tokens[3].is_separator());
		REQUIRE(tokens[4].get_text() == "B");
		REQUIRE(!tokens[4].is_separator());
		REQUIRE(tokens[5].get_text() == ",");
		REQUIRE(tokens[5].is_separator());
		REQUIRE(tokens[6].get_text() == "D");
		REQUIRE(!tokens[6].is_separator());
		REQUIRE(tokens[7].get_text() == ")");
		REQUIRE(tokens[7].is_separator());
		REQUIRE(tokens[8].get_text() == ",");
		REQUIRE(tokens[8].is_separator());
		REQUIRE(tokens[9].get_text() == "(");
		REQUIRE(tokens[9].is_separator());
		REQUIRE(tokens[10].get_text() == "A");
		REQUIRE(!tokens[10].is_separator());
		REQUIRE(tokens[11].get_text() == ",");
		REQUIRE(tokens[11].is_separator());
		REQUIRE(tokens[12].get_text() == "E");
		REQUIRE(!tokens[12].is_separator());
		REQUIRE(tokens[13].get_text() == ")");
		REQUIRE(tokens[13].is_separator());
		REQUIRE(tokens[14].get_text() == ")");
		REQUIRE(tokens[14].is_separator());
		REQUIRE(tokens[15].get_text() == ";");
		REQUIRE(tokens[15].is_separator());
		REQUIRE(tokens.size() == 16);		
	}
	 
	SECTION("Test tokenizer nodes") {
		Tokenizer tokenizer;
		auto tokens = tokenizer.tokenize("A B  C      D E");
		REQUIRE(tokens[0].get_text() == "A");
		REQUIRE(tokens[1].get_text() == " ");
		REQUIRE(tokens[2].get_text() == "B");
		REQUIRE(tokens[3].get_text() == " ");
		REQUIRE(tokens[4].get_text() == "C");
		REQUIRE(tokens[5].get_text() == " ");
		REQUIRE(tokens[6].get_text() == "D");
		REQUIRE(tokens[7].get_text() == " ");
		REQUIRE(tokens[8].get_text() == "E");
		REQUIRE(tokens.size() == 9);
	}
	

 }