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
# EDIT 	Edit Distance http://rosalind.info/problems/edit/



if __name__=='__main__':
    with open(r'C:\Users\Simon\Downloads\rosalind_ba5g(2).txt','r') as f:
        strings = []
        for line in f:
            strings.append(line.strip())
            print (line.strip())
        d,_=edit(strings[0],strings[1])
        print(d)
            
        