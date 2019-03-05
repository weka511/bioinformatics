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
# BA6I Implement Graph to Genome

from BA6G import CycleToChromosome

def GraphToGenome(GenomeGraph):
    def diff(a,b):
        _,x=a
        y,_=b
        return abs(x-y)
    def build_cycle(pair,dcycles):
        result = [pair]
        while pair in dcycles:
            pair = dcycles[pair]
            result.append(pair)
        return result
    def black_edges(cycle):
        result=[]
        for i in range(len(cycle)):
            a,_ = cycle[i]
            _,b = cycle[i-1]
            result.append((b,a))
        return result
    #for  Nodes in GenomeGraph:
        #P.append( CycleToChromosome(Nodes))
    extract = [(a,b,diff(a,b)) for (a,b) in zip([GenomeGraph[-1]]+GenomeGraph[0:],GenomeGraph)]
    gaps    = [(a,b) for (a,b,diff) in extract if diff>1]
    cycles  = [(a,b) for (a,b,diff) in extract if diff==1]
    dcycles = {}
    for (a,b) in cycles:
        dcycles[a]=b
    P       = [build_cycle(pair,dcycles) for _,pair in gaps]
    Q       = [black_edges(p) for p in P]
    next_node = 1
    R = []
    for q in Q:
        r = []
        for a,b in q:
            r.append(next_node if a<b else -next_node)
            next_node+=1
        R.append(r)
    return R

if __name__=='__main__':
    print (GraphToGenome([(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]))
    
    print (GraphToGenome([(1, 4), (3, 6), (5, 7), (8, 9), (10, 11), (12, 14), (13, 15), (16, 18), (17, 19), (20, 22), (21, 24), (23, 26), (25, 27), (28, 29), (30, 32), (31, 33), (34, 36), (35, 38), (37, 39), (40, 42), (41, 44), (43, 2), (45, 48), (47, 50), (49, 51), (52, 54), (53, 56), (55, 57), (58, 60), (59, 62), (61, 63), (64, 66), (65, 68), (67, 69), (70, 72), (71, 74), (73, 76), (75, 78), (77, 79), (80, 81), (82, 84), (83, 86), (85, 87), (88, 90), (89, 91), (92, 94), (93, 46), (96, 98), (97, 99), (100, 102), (101, 103), (104, 106), (105, 107), (108, 109), (110, 111), (112, 114), (113, 115), (116, 117), (118, 119), (120, 122), (121, 123), (124, 125), (126, 127), (128, 130), (129, 131), (132, 134), (133, 135), (136, 138), (137, 139), (140, 141), (142, 144), (143, 145), (146, 147), (148, 149), (150, 152), (151, 95), (154, 155), (156, 158), (157, 159), (160, 161), (162, 163), (164, 166), (165, 168), (167, 169), (170, 171), (172, 174), (173, 175), (176, 178), (177, 180), (179, 182), (181, 184), (183, 186), (185, 187), (188, 190), (189, 191), (192, 194), (193, 196), (195, 153), (198, 200), (199, 201), (202, 204), (203, 205), (206, 207), (208, 210), (209, 211), (212, 214), (213, 215), (216, 217), (218, 219), (220, 222), (221, 223), (224, 226), (225, 227), (228, 230), (229, 232), (231, 234), (233, 235), (236, 238), (237, 239), (240, 241), (242, 244), (243, 245), (246, 197), (248, 250), (249, 252), (251, 253), (254, 255), (256, 258), (257, 260), (259, 262), (261, 263), (264, 265), (266, 267), (268, 270), (269, 272), (271, 273), (274, 276), (275, 278), (277, 280), (279, 281), (282, 283), (284, 286), (285, 247), (288, 290), (289, 292), (291, 293), (294, 295), (296, 297), (298, 299), (300, 302), (301, 303), (304, 306), (305, 307), (308, 310), (309, 312), (311, 314), (313, 316), (315, 318), (317, 319), (320, 322), (321, 323), (324, 325), (326, 327), (328, 330), (329, 332), (331, 333), (334, 336), (335, 338), (337, 339), (340, 287)]))