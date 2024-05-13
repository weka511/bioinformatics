/**
 * Copyright (C) 2024 Simon Crase: simon@greenweavez.nz
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
#include "kmin.h"

TEST_CASE( "KMIN tests", "[kmin]" ) {
	
   	SECTION("Sample"){
		EditDistanceCalculator calc;
		std::vector< std::pair<int, int>> Result = calc.get_distance(2,"ACGTAG","ACGGATCGGCATCGT");
		REQUIRE(Result.size()==3);
	} 
}

