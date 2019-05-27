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
#    BA4M 	Solve the Turnpike Problem   


import math,numpy
from rosalind import read_list,write_list,RosalindException
from spectrum import Turnpike


if __name__=='__main__':
    try:
        write_list (Turnpike(read_list('/Users/Simon/Downloads/rosalind_ba4m_5_dataset.txt'),
                             check=True),
                    out='BA4M.txt')
    except RosalindException as e:
        print (e.args[0])