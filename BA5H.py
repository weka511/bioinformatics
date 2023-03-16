#!/usr/bin/env python

#    Copyright (C) 2019-2023 Greenweaves Software Limited
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

''' BA5H Find a Highest-Scoring Fitting Alignment of Two Strings '''

from timeit import default_timer

import numpy as np

from reference_tables import createSimpleDNASubst
from align            import FindHighestScoringFittingAlignment
from helpers          import create_strings


if __name__=='__main__':
      start_time = default_timer()
      strings    = create_strings(path='data')
      d,s1,t1    = FindHighestScoringFittingAlignment(strings[0],strings[1])
      print ('Score = {0}'.format(d))
      print (s1)
      print (t1)
      print ('Elapsed Time = {0}'.format(default_timer() - start_time))
