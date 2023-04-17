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
 
#ifndef _QRTD_H
#define _QRTD_H

#include <string>
#include "tree.h"

class QuartetDistanceCalculator {
  public:
	QuartetDistanceCalculator(Taxa & taxa) : _taxa(taxa) {}
	int get_distance(std::string T1, std::string T2);
	
  private:
	Taxa & _taxa;
    int _get_distance(Tree* T1, Tree* T2);
	void prepare( Tree* T);
};
#endif //_QRTD_H
