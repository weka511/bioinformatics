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
    ss =DecodeIdealSpectrum([103, 131, 259, 287, 387, 390, 489, 490, 577, 636, 690, 693, 761, 840, 892, 941, 1020, 1070, 1176, 1198, 1247, 1295, 1334, 1462, 1481, 1580, 1599, 1743, 1762, 1842, 1861, 2005, 2024, 2123, 2142, 2270, 2309, 2357, 2406, 2428, 2534, 2584, 2663, 2712, 2764, 2843, 2911, 2914, 2968, 3027, 3114, 3115, 3214, 3217, 3317, 3345, 3473, 3501, 3604])
    print (ss)
    # Should be:
    # CRQCSLAMQRASQHYVYVWPQETFGFVCRM
