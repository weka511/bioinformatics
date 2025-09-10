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
 */
 
#ifndef _TEST_ADAPTER_HPP
#define _TEST_ADAPTER_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "problem.hpp"

using namespace std;



class TestOutput : public OutputAdapter, public vector<string>{
   public:
	void append(string text);
	void append(vector<int> counts);
	string get(int i);
};




#endif //_TEST_ADAPTER_HPP
