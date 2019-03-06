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

from fragile import d2break


if __name__=='__main__':
    def conv(xx):
        return [int(s) for s in xx]
    def parse(line):
        ff=line[1:-1].split(')(')
        result=[conv(fff.split(' ')) for fff in ff]
        return result
    with open('/Users/Simon/Downloads/rosalind_ba6c.txt') as f:
 
        print (
            d2break(
               parse(f.readline().strip()),
                parse(f.readline().strip())
        ))
