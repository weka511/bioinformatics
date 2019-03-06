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

# BA6G Implement Cycle to Chromosome 

from fragile import CycleToChromosome,ChromosomeToCycle


if __name__=='__main__':
    def disp(i):
        return str(i) if i < 0 else '+' + str(i)
    cycle = [1, 2, 4, 3, 5, 6, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 17, 18, 19, 20, 22, 21, 23, 24, 25, 26, 28, 27, 30, 29, 32, 31, 33, 34, 36, 35, 37, 38, 39, 40, 42, 41, 44, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 55, 57, 58, 59, 60, 62, 61, 64, 63, 65, 66, 68, 67, 69, 70, 72, 71, 73, 74, 76, 75, 77, 78, 79, 80, 82, 81, 83, 84, 85, 86, 88, 87, 89, 90, 92, 91, 93, 94, 95, 96, 98, 97, 100, 99, 102, 101, 103, 104, 106, 105, 108, 107, 109, 110, 112, 111, 114, 113, 116, 115, 117, 118, 120, 119, 121, 122, 123, 124, 126, 125, 128, 127, 129, 130, 131, 132, 134, 133, 135, 136, 138, 137, 140, 139]
    Chromosome =CycleToChromosome(cycle)
    
    # Test - invert transformation and compare
    cycle1 = ChromosomeToCycle(Chromosome)
    for i in range(max(len(cycle),len(cycle1))):
        if cycle[i] != cycle1[i]:
            print (i, cycle[i],cycle1[i])
    digits = ' '.join(disp(c) for c in Chromosome)
    print ('(' + digits + ')')