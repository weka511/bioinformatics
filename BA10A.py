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

def ProbabilityHiddenPath(path,Transition):
    def prob():
        result=1
        for i in range(1,len(path)):
            result *= Transition[(path[i-1],path[i])]
 
        return result
    return 0.5*prob()

if __name__=='__main__':
    Transition = {
        ('A','A'): 0.194, ('A','B'):   0.806,
        ('B','A'): 0.273, ('B','B'):   0.727,
    }
    print (ProbabilityHiddenPath('AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB',Transition))
    
    Transition1 = {
        ('A','A'): 0.497, ('A','B'):   0.503,
        ('B','A'): 0.263, ('B','B'):   0.737,
    }
    print (ProbabilityHiddenPath('BBABBAAABBAABBBBABABBBBAAABBBBBBAABAABAAABABAABBBA',Transition1))    