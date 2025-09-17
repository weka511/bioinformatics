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
 * You should   received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 *
 * Classes that allow Rosalind problems to be accessed from unit tests
 */
 
 #include <sstream>
 #include "test-adapter.hpp"
 
 using namespace std;
 
 void TestOutput::append(string text) {
	 push_back(text);
}
	


string TestOutput::get(int i) {
	return (*this)[i];
}