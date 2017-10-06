#    Copyright (C) 2017 Greenweaves Software Pty Ltd
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
#    #BA4M 	Solve the Turnpike Problem   


import math,numpy
from rosalind import read_list,write_list,RosalindException

def Turnpike(D,check=False):
    '''
    Solve the Turnpike Problem
    Parameters:
       D     The differences
       check Indicates whether differences are to be checked 
    
    Based Mark Weiss's treatment
    https://users.cs.fiu.edu/~weiss/cop3337_f99/assignments/turnpike.pdf
    '''
    
    def find_remaining_points(X,D,first,last):
        '''
        Extend a partial solution by adding more points.
          Parameters
              X         Partial solution - typically contains NaNs for unassigned points,
                                           with known values filled in from the ends
              D
              first     Index of end of knownn values at beginning X
              last      Index of start of know value at end of X
        '''
        def get_set_diffs(diffs):
            '''
            Ascertain whether a some partial set of differences is a subset of D
            '''
            diffs.sort()
            set_diffs=[]
            i=0
            for diff in diffs:
                found=False
                while i<len(D):
                    if  D[i]<diff:
                        set_diffs.append(D[i])
                    elif D[i]==diff:
                        found=True
                        i+=1
                        break
                    i+=1
                if not found:
                    return None
            while i<len(D):
                set_diffs.append(D[i])
                i+=1 
            return set_diffs
        
        def explore(candidate,X,first,last):
            '''
            Explore lower levels of tree
            '''
            # Constuct set of differences between candidate and known members of X
            diffs=[abs(candidate-x) for x in X if not numpy.isnan(x) and x!=candidate]
            set_diffs=get_set_diffs(diffs)
            if set_diffs==None:
                return None
            elif len(set_diffs)==0:
                return X
            else:
                return find_remaining_points(X,set_diffs,first,last)
        
        # There are two cases to consider: either the largest remaining unprocessed value in D
        # is part of the solution, or it isn't. We will explore these two cases separately.
        # We maintain a tree of data structures, so we can explore the tree of solutions
        # see https://users.cs.fiu.edu/~weiss/cop3337_f99/assignments/turnpike.pdf for details
        x_max=D[-1]
        XX=X[:]           # Clone this so the 
        XX[last-1]=x_max  # Add candidate at end
        trial_solution=explore(x_max,XX,first,last-1) #process level below - one fewer unknown at end
        if trial_solution==None:  # largest remaining unprocessed value was a false lead
            XX=X[:]
            XX[first+1]=X[-1]-x_max # Added new candiate at beginning
            return explore(X[-1]-x_max,XX,first+1,last) #process level below - one fewer unknown at start
        else:
            return trial_solution
    
    def check_diffs(reconstruction):
        '''
        Verify that a particular reconstruction does indeed give rise to original differences
        '''
        diffs=[a-b for a in reconstruction for b in reconstruction]
        diffs.sort()
        if len(diffs)!=len(D):
            raise RosalindException ('Length of reconstructed diffs ({0}) does not match length of D ({1}) '.format(len(diffs), len(D)))
        mismatches=0
        for a,b in zip(D,diffs):
            if a!=b:
                mismatches+=1
                print (a,b)
       
        if mismatches>0:
            raise RosalindException ('Found {0} mismatches'.format(mismatches))
        return diffs
     
    # Start by initializing array of points. We know that its length
    # must be the square root of the array of differences.   
    len_D=len (D)
    len_X=int(math.sqrt(len_D))
    X=[float('nan')]*len_X    #We fill in all values as "unknown" 
    X[0]=0                    #Actually we are given the first value, zero
    X[-1]=D[-1]               # We also know that the last point must match the last difference.
    reconstruction= find_remaining_points(X,[d for d in D[:-1] if d>0],0,-1)
    if check:
        check_diffs(reconstruction)
    return reconstruction

if __name__=='__main__':
    try:
        write_list (Turnpike(read_list('data/rosalind_ba4m(2).txt'),check=True),out='BA4M.txt')
    except RosalindException as e:
        print (e.args[0])