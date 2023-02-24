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
# along with this pogram.  If not, see <http://www.gnu.org/licenses/>

'''
   BA10B Compute the Probability of an Outcome Given a Hidden Path
'''

import numpy as np

from hmm import ProbabilityOutcomeGivenHiddenPath

if __name__=='__main__':

    print (ProbabilityOutcomeGivenHiddenPath('xzyyxyzyxzyzxxzyzzxxzyxyxyxzxzyyzyzxzxzyxyzyyyzzxz', 'xyz',
                                             'BBBAABAAABBABABBABBAAAAAABAAABAABABBBABBAABABAABAB', 'AB',
                                             np.array([[0.419,0.321, 0.26],
                                                      [0.185, 0.551,0.263]])))


    print (ProbabilityOutcomeGivenHiddenPath('zyyyxzxzyyzxyxxyyzyzzxyxyxxxxzxzxzxxzyzzzzyyxzxxxy', 'xyz',
                                             'BAABBAABAABAAABAABBABBAAABBBABBAAAABAAAABBAAABABAA', 'AB',
                                                                                          np.array([[ 0.093, 0.581, 0.325],
                                                                                                   [0.77, 0.21,0.02]])))

