#!/usr/bin/env python

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
#    BA6K Implement 2-BreakOnGenome


from fragile import perform_2_BreakOnGenome



if __name__=='__main__':
     #print (perform_2_BreakOnGenome([+1, -2, -4, +3], 1, 6, 3, 8))
     def format(genome):
          def ff(g):
               return '+'+str(g) if g>0 else str(g)
          def f(g):
               return '('+ ' '.join([ff(g0) for g0 in g]) + ')'
          return ' '.join([f(g) for g in genome])

     print (
          format(
               perform_2_BreakOnGenome(
                    [-1, +2, +3, -4, -5, +6, +7, -8, -9, +10, +11, -12, -13, +14, +15, +16, +17, -18, +19, +20, -21, -22, -23, -24, +25, -26, -27, +28, +29, -30, -31, -32, +33, -34, -35, +36, +37, -38, +39, +40, -41, +42, -43, -44, -45, +46, +47, +48, -49, -50, -51, -52, +53, -54, -55, -56, -57, -58, -59, -60, -61, -62, +63, +64, +65, -66, +67],
                    17, 19, 67, 70
               )))
