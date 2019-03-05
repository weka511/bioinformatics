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
#    BA6J Implement 2-BreakOnGenomeGraph

def get2BreakOnGenomeGraph(graph,i0,i1,j0,j1):
    def eq(x,y):
        u,v=x
        w,z=y
        return (u ==w and v==z) or (w==v and u==z)
    removed = [x for x in graph if not eq(x,(i0,i1)) and not eq(x,(j0,j1))]
    return removed +[(i0,j0)] + [(i1,j1)]

if __name__=='__main__':
    #print (get2BreakOnGenomeGraph([(2, 4), (3, 8), (7, 5), (6, 1)],1, 6, 3, 8))
    #print (get2BreakOnGenomeGraph([(2, 4), (3, 5), (6, 8), (7, 10), (9, 12), (11, 13), (14, 15), (16, 17), (18, 19), (20, 21), (22, 24), (23, 25), (26, 28), (27, 30), (29, 32), (31, 33), (34, 35), (36, 37), (38, 40), (39, 42), (41, 44), (43, 45), (46, 47), (48, 49), (50, 51), (52, 53), (54, 56), (55, 58), (57, 59), (60, 61), (62, 64), (63, 66), (65, 68), (67, 69), (70, 72), (71, 73), (74, 75), (76, 77), (78, 80), (79, 81), (82, 83), (84, 86), (85, 88), (87, 90), (89, 92), (91, 94), (93, 96), (95, 98), (97, 99), (100, 102), (101, 104), (103, 106), (105, 107), (108, 109), (110, 111), (112, 113), (114, 116), (115, 117), (118, 120), (119, 121), (122, 124), (123, 126), (125, 128), (127, 129), (130, 132), (131, 134), (133, 136), (135, 1)],
                                  #87, 90, 74, 75))
    print (get2BreakOnGenomeGraph(
        [(2, 3), (4, 6), (5, 8), (7, 10), (9, 11), (12, 14), (13, 15), (16, 18), (17, 20), 
         (19, 22), (21, 23), (24, 25), (26, 27), (28, 29), (30, 32), (31, 34), (33, 36), (35, 37), (38, 40), (39, 42), (41, 44), 
         (43, 45), (46, 48), (47, 50), (49, 51), (52, 53), (54, 55), (56, 57), (58, 59), (60, 61), (62, 64), (63, 66), (65, 68), (67, 69), (70, 71), (72, 73), (74, 76), (75, 78), (77, 79), (80, 82), (81, 83), (84, 85), (86, 87), (88, 90), (89, 92), (91, 94), (93, 96), (95, 98), (97, 100), (99, 102), (101, 103), (104, 105), (106, 107), (108, 109), (110, 112), (111, 113), (114, 115), (116, 117), (118, 119), (120, 121), (122, 123), (124, 126), (125, 128), (127, 129), (130, 132), (131, 1)
],13, 15, 86, 87
    
    ))