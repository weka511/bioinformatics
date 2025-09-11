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
#include <string>
#include <vector>

using namespace std;

class Datasource {
   public:
	virtual string get_input(int i)=0;
};

class OutputAdapter  {
   public:
	virtual void append(string text)=0;
	
	void append(vector<int> counts);
	
	virtual void flush() {;};
};

/**
 * Abstract class that serves as parent for solutions to problems.
 */
class Problem{
  private:
    Datasource* _datasource;
	OutputAdapter* _output;
	string _name;
	
  protected:
    Problem(string name) : _name(name) {};
	
	string get_input(int i) { 
		return _datasource->get_input(i);
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
    virtual void solve() = 0;
	
	void attach(Datasource* datasource) { _datasource = datasource;}
	
	void attach(OutputAdapter* output) { _output = output;}
	
};

/**
 * Used by factory to instantiate problems
 */
template <class T>
Problem* createProblem() {
	return new T;
};


#endif //_PROBLEM_HPP
