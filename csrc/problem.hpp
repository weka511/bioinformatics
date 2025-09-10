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
 
#ifndef _PROBLEM_HPP
#define _PROBLEM_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;

class Datasource : public vector<string>{
   public:
	
};

class Output : public vector<string>{
   public:
	void append(string text) {push_back(text);};
	void append(vector<int> counts);
	string get(int i) {return (*this)[i];};
};

/**
 * Abstract class that serves as parent for solutions to problems.
 */
class Problem{
  private:
    Datasource* _datasource;
	Output* _output;
	string _name;
	
  protected:
    Problem(string name) : _name(name) {};
	
	string get_input(int i) { 
		return(*_datasource)[i];
	}
	
	void append(vector<int> counts){
		_output->append(counts);
	}
	
	void append(string text){
		_output->append(text);
	}
  public:
	/**
	 * Solve a problem
     */
    virtual void solve()=0;
	
	void attach(Datasource* datasource) { _datasource = datasource;}
	
	void attach(Output* output) { _output = output;}
	
	string get(int i) {return _output->get(i);};
};

/**
 * Used by factory to instantiate problems
 */
template <class T>
Problem* createProblem() {
	return new T;
};


#endif //_PROBLEM_HPP
