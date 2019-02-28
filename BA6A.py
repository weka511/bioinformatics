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
# BA6A Implement GreedySorting to Sort a Permutation by Reversals

def GreedySorting(P):
    def kReverse(k):
        pos = P.index(k if k in P else -k)
        return P[0:k-1] + [-P[j] for j in range(pos,k-2,-1)] + P[pos+1:]
    def format():
        def f(p):
            return str(p) if p<0 else '+' + str(p)
        return '(' + ' '.join(f(p) for p in P) + ')'
    reversalDistance = 0
    for k in range(1,len(P)+1):
        if k!=P[k-1]:
            P=kReverse(k)
            reversalDistance+=1
            print (format())
            if P[k-1]==-k:
                P[k-1]=k
                reversalDistance+=1
                print (format())
    return reversalDistance

if __name__=='__main__':
    print ( GreedySorting([-3, +4, +1, +5, -2]))
