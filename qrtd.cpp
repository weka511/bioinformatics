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

#include "qrtd.h"

Taxa::Taxa(std::string taxa_string){
	std::istringstream taxa_stream(taxa_string);
	while (taxa_stream) {
		std::string taxon_name;
		taxa_stream >> taxon_name;
		if (taxon_name.size()>0)
			_names.push_back(taxon_name);
	}
}

bool ConsistencyChecker::visit(Clade * clade) {
	if (clade->get_name().size()==0) return true;
	if (_counts.count(clade->get_name())==1)
		_counts[clade->get_name()]++;
	else
		std::cout<< "Key \"" << clade->get_name() << "\" in tree is not present in list of taxa" << std::endl;
	return true;
}

bool ConsistencyChecker::is_consistent(Clade * root,Taxa& taxa){
	_counts.clear();
	for (std::vector<std::string>::iterator name=taxa._names.begin(); name!=taxa._names.end(); name++)
		_counts[*name] =0;
	bool result = root->depth_first_search(this);
	for (std::map<std::string, int>::iterator it = _counts.begin(); it != _counts.end(); it++)
		if (it->second!=1){
			std::cout << "Taxon \"" << it->first << "\" has count " << it->second << " - should be 1." << std::endl;
			result = false;
		}
	return result;
}

int get_qrtd(std::string taxa_string, std::string T1, std::string T2) {
	Taxa taxa(taxa_string);

	return -1;
}
