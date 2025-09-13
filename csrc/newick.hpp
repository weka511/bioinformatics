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
 
#ifndef _NEWICK_HPP
#define _NEWICK_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "tokenizer.hpp"

using namespace std;

class Node{
  private:
	vector<shared_ptr<Node>> nodes;
};
	
class Newick {
  public:

	 
	void parse(string s);
	
	tuple<shared_ptr<Node>,vector<tuple<int,int>>> explore(vector<Token> tokens, const int from, const int to);
	
  private:
	logic_error _create_error(const Token token, const int depth);
};



#endif // _NEWICK_HPP
