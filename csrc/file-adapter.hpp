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
 
#ifndef _FILE_ADAPTER_HPP
#define _FILE_ADAPTER_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "problem.hpp"

using namespace std;



class FileNameFactory {
  private:
	string path="C:\\Users\\Weka\\Downloads\\";
	
  public:
	enum Format{SUBMIT,TEST};
		
	string create(string problem_name,Format format,int seq);
};

class FileDatasource : public Datasource{
	private: 
	
   public:
	FileDatasource(string file_name);
	string get_input(int i) { 
		return "";
	}
};

class FileOutput : public OutputAdapter{
  private:
 
	
  public:
	FileOutput(string file_name);
   
	void append(string text);
	void append(vector<int> counts);
};




#endif //_FILE_ADAPTER_HPP
