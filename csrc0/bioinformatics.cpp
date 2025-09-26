/**
 * Copyright (C) 2023 Simon Crase: simon@greenweavez.nz
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
#include <fstream> 
#include "tree.h"
#include "newick.h"
#include "qrtd.h"



int main(){
	std::string text1, text2, text3;
	std::ifstream data_file("data/rosalind_qrtd.txt");
	getline (data_file, text1);
	Taxa taxa(text1);
	std::cout <<__FILE__ << " " <<__LINE__ << " " << taxa.size() << std::endl;
	getline (data_file, text2);
	getline (data_file, text3);
	data_file.close();
	QuartetDistanceCalculator calculator(taxa);
	int dist = calculator.get_distance(text2,text3);
	std::cout <<__FILE__ << " " <<__LINE__ << " " << dist << std::endl;
	return 0;
}
 