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
#    HEA buiding a heap

def hea(n,A):
    def next_free(T):
        pass
    def insert(a,T,current):
        if len(T)==0:
            T.append((-1,1,2,a))
            return 0
        else:
            parent,left,right,_=T[current]
            if len(T)<=left:
                T.append((current,2*left+1,2*left+2,a))
                return current
            if len(T)<=right:
                T.append((current,2*right+1,2*right+2,a))
                return next_free(T)            
    def linearize(T):
        return []
    T   = []
    ptr = -1
    
    for a in A:
        ptr=insert(a,T,ptr)
    return linearize(T)

if __name__=='__main__':
    print (hea(5,[1, 3, 5, 7, 20]))