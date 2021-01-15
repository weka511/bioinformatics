#    Copyright (C) 2019-2021 Greenweaves Software Limited
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#    BA6A Implement GreedySorting to Sort a Permutation by Reversals

from   fragile import GreedySorting
import argparse
import os
import time
from   helpers import read_strings

def get_permutation(Input):
    return [int(i) for i in Input[1:-1].split()]

def format(P,signed=False):
    def f(p):
        return str(p) if not signed or p<0 else '+' + str(p)
    return '(' + ' '.join(f(p) for p in P) + ')'

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('Find a Shortest Transformation of One Genome into Another by 2-Breaks')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        for permutation in GreedySorting(  [
            -132, -42, +143, +38, -12, -79, -33, -48, -5, -19, -9, -1, -100, -92, +65, -70, +109, -123, +18, -118, -7, -78,
            -145, -6, +31, +36, -102, -58, -46, -53, -138, +73, +75, +113, -13, -124, -43, -24, +3, +114, +20, +99, -96, -15,
            -121, -98, +134, +104, +125, +22, -64, -68, +94, -30, +110, +103, +57, +116, -85, +130, -44, -131, +142, -50, 
            +107, -52, +89, +63, +144, -47, +14, -133, +8, -72, +4, -37, +29, +10, -90, -97, -83, -106, -105, -60,
            -108, +119, +122, -2, +32, -112, -61, -111, +141, +84, -126, +82, +88, -95, +40, -129, -137, -67, +17,
            -54, +66, +74, -16, -93, -140, +86, -115, +80, +35, +139, -91, +39, +41, -26, +76, -21, +62, -87, -27,
            -117, -71, -28, -77, -59, -11, -45, +69, -120, +128, +51, +101, -25, +34, +136, -49, -55, +81, +56, -127,
            +135, -23
            ],
            signed=True):
            print (permutation)

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for permutation in GreedySorting(get_permutation(Input[0])):
                print (format(permutation))
                f.write(f'{format(permutation)}\n')
        
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')  