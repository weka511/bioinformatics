#!/usr/bin/env python

#    Copyright (C) 2019 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
#    BA11A Construct the Graph of a Spectrum

from spectrum import SpectrumGraph

if __name__=='__main__':
    graph = SpectrumGraph([87, 137, 200, 208, 287, 355, 400, 483, 529, 612, 685, 725, 832, 856, 959, 960, 1063, 1072, 1120, 1159, 1221, 1260, 1308, 1317, 1420, 1421, 1524, 1548, 1655, 1695, 1768, 1851, 1897, 1980, 2025, 2093, 2172, 2180, 2243, 2293, 2380])
    for k in sorted(graph.keys()):
        for w,p in graph[k]:
            print ('{0}->{1}:{2}'.format(k,w,p))
