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
#include "qrtd.h"

TEST_CASE( "QRTD tests", "[qrtd]" ) {
	
 /*  	SECTION("Sample"){
		REQUIRE(get_qrtd("A B C D E",
			"(A,C,((B,D),E));",
			"(C,(B,D),(A,E));")==4);
	}  */
}

TEST_CASE( "Consitency tests", "[consist]" ) {	
 	SECTION("Consistency"){
		Newick             newick;
		Clade *            root = newick.parse("(A,C,((B,D),E));");
		Taxa               taxa("A B C D E");
		ConsistencyChecker consistency_checker;
		REQUIRE(consistency_checker.is_consistent(root,taxa));
		delete root;
	}
	
	SECTION("Consistency"){
		Newick             newick;
		Clade *            root = newick.parse("(A,C,((B,D),F));");
		Taxa               taxa("A B C D E");
		ConsistencyChecker consistency_checker;
		REQUIRE_FALSE(consistency_checker.is_consistent(root,taxa));
		delete root;
	}
	
}
