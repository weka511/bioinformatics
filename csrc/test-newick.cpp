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
 */
 
 #include "catch.hpp"
 #include "newick.hpp"
  
 TEST_CASE( "Newick Tests", "[newick]" ) {
	 
	 SECTION("Test tokenizer newick") {
		Tokenizer tokenizer;
		auto tokens = tokenizer.tokenize("(C,(B,D),(A,E));");
		REQUIRE(tokens[0] == "(");
		REQUIRE(tokens[1] == "C");
		REQUIRE(tokens[2] == ",");
		REQUIRE(tokens[3] == "(");
		REQUIRE(tokens[4] == "B");
		REQUIRE(tokens[5] == ",");
		REQUIRE(tokens[6] == "D");
		REQUIRE(tokens[7] == ")");
		REQUIRE(tokens[8] == ",");
		REQUIRE(tokens[9] == "(");
		REQUIRE(tokens[10] == "A");
		REQUIRE(tokens[11] == ",");
		REQUIRE(tokens[12] == "E");
		REQUIRE(tokens[13] == ")");
		REQUIRE(tokens[14] == ")");
		REQUIRE(tokens[15] == ";");
		REQUIRE(tokens.size() == 16);		
	 }
	 
	SECTION("Test tokenizer nodes") {
		Tokenizer tokenizer;
		auto tokens = tokenizer.tokenize("A B C D E");
		REQUIRE(tokens[0] == "A");
		REQUIRE(tokens[1] == " ");
		REQUIRE(tokens[2] == "B");
		REQUIRE(tokens[3] == " ");
		REQUIRE(tokens[4] == "C");
		REQUIRE(tokens[5] == " ");
		REQUIRE(tokens[6] == "D");
		REQUIRE(tokens[7] == " ");
		REQUIRE(tokens[8] == "E");
		REQUIRE(tokens.size() == 9);
	 }
 }