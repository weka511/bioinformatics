/**
 * Copyright (C) 2024 Simon Crase: simon@greenweavez.nz
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
 
#ifndef _KMIN_H
#define _KMIN_H

#include <string>
#include <vector>


class EditDistanceCalculator {
	
  public:
	EditDistanceCalculator() {};
	std::vector< std::pair<int, int>>  get_distance(int k, std::string s, std::string t) ;
};
#endif //_KMIN_H
