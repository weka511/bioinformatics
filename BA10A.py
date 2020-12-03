# Copyright (C) 2019-2020 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

# BA10A Compute the Probability of a Hidden Path

from hmm import ProbabilityHiddenPath

if __name__=='__main__':
    print (ProbabilityHiddenPath('AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB',
                                 {                                #sample
                                                                  ('A','A'): 0.194, ('A','B'):   0.806,
                                                                  ('B','A'): 0.273, ('B','B'):   0.727,
                                                                  }))
    
    print (ProbabilityHiddenPath('BBABBBABBAABABABBBAABBBBAAABABABAAAABBBBBAABBABABB',
                                 {                               #extra
                                                                 ('A','A'): 0.863, ('A','B'):   0.137,
                                                                 ('B','A'): 0.511, ('B','B'):   0.489,
                                                                 }))    