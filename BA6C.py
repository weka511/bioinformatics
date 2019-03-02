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
#    BA6C Compute the 2-Break Distance Between a Pair of Genomes

def ba6c(a,b):
    def blocks():
        return 0
    def cycles():
        return 0
    return blocks()-cycles()

if __name__=='__main__':
    print (ba6c(
        [[+1, +2, +3, +4, +5, +6]],
        [[+1, -3, -6, -5],[+2, -4]]
    ))
