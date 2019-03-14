# Copyright (C) 2017 Greenweaves Software Pty Ltd

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



from graphs import ddeg

#n=5
#m=4
#A=[[1, 2],
   #[2, 3],
   #[4, 3],
   #[2, 4]]

#print(ddeg(n,m,A))

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_ddeg.txt') as f:
        A=[]
        for line in f:
            text=line.strip()
            pair=text.split(' ')
            print (pair)
            A.append((int(pair[0]),int(pair[1])))
        (n,m)=A[0]

        print (ddeg(n,m,A[1:]))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))