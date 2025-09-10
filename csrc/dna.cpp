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

#include <map>

#include "dna.hpp"

using namespace std;

void DNA::solve() {
	map<char,int> counts;
	counts['A'] = 0;
	counts['C'] = 0;
	counts['G'] = 0;
	counts['T'] = 0;
	auto dna_string= get_input(0);
	for (auto base : dna_string)
		counts[base] ++;
	
	vector<int> c;
	c.push_back(counts['A']);
	c.push_back(counts['C']);
	c.push_back(counts['G']);
	c.push_back(counts['T']);
	append(c);
}
