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
 * Classes that allow Roslind problems to interact with files.
 */
 
#ifndef _FILE_ADAPTER_HPP
#define _FILE_ADAPTER_HPP

#include <string>
#include <vector>

#include "problem.hpp"

using namespace std;


/**
 * This class creates a file name for a problem
 */
class FileNameFactory {
  private:
	string _path="C:\\Users\\Weka\\Downloads\\";
	
  public:
	enum Format{
		SUBMIT, 	// rosalind_ba11f
		TEST 		// rosalind_dna_1_dataset
	};
	
	/**
	 * Create a file name for a problem
	 */
	string create(string problem_name, const Format format,int const seq=0);
};

/**
 *  This class is responsible for providing input data for a problem from a file.
 */
class FileDatasource : public Datasource{
  private: 
	vector<string> _lines;
		
  public:
	FileDatasource(string file_name);
	
	string get_input(int i);
};

/**
 *  This class is responsible for storing results data from a problem
 *  into a file.
 */
class FileOutput : public OutputAdapter{
  private:
	vector<string> _lines;
	string _file_name;
	
  public:
	FileOutput(string file_name);
   
	void append(string text);

	void flush();

};

#endif //_FILE_ADAPTER_HPP
