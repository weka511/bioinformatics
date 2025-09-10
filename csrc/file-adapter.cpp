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
 
 // rosalind_ba11f
 // rosalind_dna_1_dataset
 
#include <fstream>
#include <stdexcept>
#include <sstream>
#include "file-adapter.hpp"
 
 using namespace std;
 
 string FileNameFactory::create(string problem_name,Format format,int seq=0) {
		switch (format) {
		case SUBMIT:
			break;
		case TEST:
			break;
		}  
		return "";
};
	
 FileDatasource::FileDatasource(string file_name){
	 ifstream file(file_name);
	 if (file.is_open()) {
		 string line;
		while (getline(file, line))
			std::cout << line << std::endl;
	} else {
        stringstream message;
		message<<__FILE__ <<" " <<__LINE__<<" Error: Unable to copen file  " << file_name<<endl; 
		throw logic_error(message.str().c_str()); 
    }
 }
 
 FileOutput::FileOutput(string file_name){
 }
 
 void FileOutput::append(string text) {
	
}
	
 void FileOutput::append(vector<int> counts){

}
