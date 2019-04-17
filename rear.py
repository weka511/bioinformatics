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
#    REAR 	Reversal Distance

import time
from fragile import rear

if __name__=='__main__':
    def parse(line):
        return [int(c) for c in line.strip().split()]
    
    start_time = time.time()
    i=0
    original=''


    data = [
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [3, 1, 5, 2, 7, 4, 9, 6, 10, 8],
        
        [3, 10, 8, 2, 5, 4, 7, 1, 6, 9],
        [5, 2, 3, 1, 7, 4, 10, 8, 6, 9],
        
        [8, 6, 7, 9, 4, 1, 3, 10, 2, 5],
        [8, 2, 7, 6, 9, 1, 5, 3, 10, 4],
        
        [3, 9, 10, 4, 1, 8, 6, 7, 5, 2],
        [2, 9, 8, 5, 1, 7, 3, 4, 6, 10],
        
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]    
    ]
    
    #it = iter(data)
    #result = []
    #for a in it:
        #print ('----')
        #result.append(rear(a,next(it)))
    #print (result)
    with open (r'C:\Users\Simon\Downloads\rosalind_rear(3).txt') as f:
        result = []
        for line in f:
            if i%3==0:
                original=parse(line)
            if i%3==1:
                result.append(rear(original,parse(line)))
                print("--- {0} seconds ---".format(time.time() - start_time))
            i+=1
    
        print(result)
        