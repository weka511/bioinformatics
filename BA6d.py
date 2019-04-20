# Copyright (C) 2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# BA6D Find a Shortest Transformation of One Genome into Another by 2-Breaks 

from fragile import d2break

def FindShortestTransformation(s,t):
     return d2break(s,t)

if __name__=='__main__':
     #print (d2break([[+1, +2, +3, +4, +5, +6]],
                    #[[+1, -3, -6, -5],[+2, -4]]))
     print (d2break([[+1, -2, -3, +4]],[[+1, +2, -4, -3]]))     
     #FindShortestTransformation([+1, -2, -3, +4],[+1, +2, -4, -3])