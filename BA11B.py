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
#    BA11B Construct the Graph of a Spectrum

from spectrum import DecodeIdealSpectrum

if __name__=='__main__':
    ss =DecodeIdealSpectrum([71, 115, 228, 234, 305, 408, 414, 564, 600, 711, 729, 825, 843, 956, 971, 1055, 1084, 1158, 1247, 1261, 1348, 1375, 1461, 1503, 1574, 1618, 1689, 1731, 1817, 1844, 1931, 1945, 2034, 2108, 2137, 2221, 2236, 2349, 2367, 2463, 2481, 2592, 2628, 2778, 2784, 2887, 2958, 2964, 3077, 3121, 3192])
    print (ss)

