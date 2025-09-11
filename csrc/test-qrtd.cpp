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
 #include "qrtd.hpp"
 #include "test-adapter.hpp"
 
 TEST_CASE( "QRTD Tests", "[qrtd]" ) {
	 
	 SECTION("quartet distance") {
		QRTD qrtd;
		TestDatasource datasource;
		datasource.push_back("A B C D E");
		datasource.push_back("(A,C,((B,D),E));");
		datasource.push_back("(C,(B,D),(A,E));");
		qrtd.attach(&datasource);
		TestOutput output;
		qrtd.attach(&output);
		qrtd.solve();
		REQUIRE(output.get(0) == "4");
	 }
 }