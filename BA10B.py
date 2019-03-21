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

def ProbabilityOutcomeGivenHiddenPath(string,path,Emission):
    result = 1
    for i in range(0,len(path)):
        result *= Emission[(path[i],string[i])]
    return result  

if __name__=='__main__':
    Emission = {
        ('A','x'): 0.419, ('A','y'):   0.321, ('A','z') : 0.26,
        ('B','x'): 0.185, ('B','y'):   0.551, ('B','z') : 0.263,
    }
    print (ProbabilityOutcomeGivenHiddenPath('xzyyxyzyxzyzxxzyzzxxzyxyxyxzxzyyzyzxzxzyxyzyyyzzxz',
                                             'BBBAABAAABBABABBABBAAAAAABAAABAABABBBABBAABABAABAB',
                                             Emission))
    
  