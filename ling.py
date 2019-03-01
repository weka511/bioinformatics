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

def ling(s,a=4):
    def sub():
        def subk(k):
            kmers=set()
            for i in range(len(s)-k+1):
                print (s[i:i+k])
                kmers.add(s[i:i+k])
            return len(kmers)
        return sum(subk(k) for k in range(1,len(s)+1))
    def m(n):
        def mm(k):
            return a if k==1 else n-k+1
        return sum(mm(k) for k in range(1,n+1))

    return sub()/m(len(s))

if __name__=='__main__':
    print (ling('ATTTGGATT'))