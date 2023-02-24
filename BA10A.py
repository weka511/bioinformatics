#!/usr/bin/env python

# Copyright (C) 2019-2023 Simon Crase

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

'''
BA10A Compute the Probability of a Hidden Path
'''

import numpy as np

from hmm import ProbabilityHiddenPath

if __name__=='__main__':
    print (ProbabilityHiddenPath('AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB',
                                 'AB',
                                 np.array([[ 0.194,   0.806],[ 0.273,   0.727]])))


    print (ProbabilityHiddenPath('BBABBBABBAABABABBBAABBBBAAABABABAAAABBBBBAABBABABB',
                                 'AB',
                                 np.array([[ 0.863,   0.137],[0.511,   0.489]])))

    print (ProbabilityHiddenPath('BBABBAAABBAABBBBABABBBBAAABBBBBBAABAABAAABABAABBBA',
                                 'AB',
                                 np.array([[ 0.497,	0.503],[0.263,	0.737]])))

    print (ProbabilityHiddenPath('ABBBAABBAABABBABBAABBBABABABAABAABBBBAAAAAAABABAAB',
                                 'AB',
                                 np.array([[0.849,	0.151],[0.43,	0.57]])))
