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
    #print (GraphToGenome([
        #(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)
    
    #]))
    
    print (GraphToGenome([
        
    (1, 3), (4, 5), (6, 8), (7, 9), (10, 12), (11, 14), (13, 15), (16, 18), (17, 19), (20, 21), (22, 23), (24, 25), (26, 28), (27, 30), (29, 31), (32, 34), (33, 36), (35, 38), (37, 39), (40, 41), (42, 44), (43, 45), (46, 47), (48, 50), (49, 52), (51, 53), (54, 55), (56, 2), (57, 60), (59, 62), (61, 63), (64, 65), (66, 68), (67, 70), (69, 71), (72, 73), (74, 76), (75, 78), (77, 80), (79, 81), (82, 83), (84, 86), (85, 88), (87, 90), (89, 91), (92, 94), (93, 95), (96, 98), (97, 99), (100, 102), (101, 104), (103, 105), (106, 107), (108, 110), (109, 58), (112, 114), (113, 115), (116, 118), (117, 119), (120, 121), (122, 123), (124, 126), (125, 128), (127, 130), (129, 131), (132, 134), (133, 136), (135, 138), (137, 139), (140, 142), (141, 144), (143, 146), (145, 148), (147, 149), (150, 151), (152, 154), (153, 155), (156, 158), (157, 160), (159, 162), (161, 163), (164, 165), (166, 167), (168, 111), (169, 171), (172, 173), (174, 176), (175, 177), (178, 179), (180, 182), (181, 184), (183, 185), (186, 188), (187, 189), (190, 191), (192, 194), (193, 196), (195, 197), (198, 200), (199, 201), (202, 204), (203, 206), (205, 208), (207, 210), (209, 212), (211, 214), (213, 216), (215, 218), (217, 219), (220, 170), (221, 224), (223, 226), (225, 228), (227, 230), (229, 231), (232, 233), (234, 236), (235, 238), (237, 239), (240, 242), (241, 243), (244, 246), (245, 248), (247, 249), (250, 252), (251, 254), (253, 256), (255, 258), (257, 259), (260, 261), (262, 264), (263, 265), (266, 268), (267, 270), (269, 222), (272, 274), (273, 275), (276, 278), (277, 280), (279, 281), (282, 284), (283, 286), (285, 287), (288, 290), (289, 292), (291, 293), (294, 296), (295, 297), (298, 299), (300, 301), (302, 304), (303, 305), (306, 308), (307, 309), (310, 311), (312, 313), (314, 315), (316, 317), (318, 320), (319, 321), (322, 324), (323, 271), (326, 327), (328, 329), (330, 332), (331, 334), (333, 335), (336, 338), (337, 340), (339, 341), (342, 344), (343, 345), (346, 347), (348, 349), (350, 352), (351, 353), (354, 355), (356, 358), (357, 360), (359, 362), (361, 364), (363, 325), (365, 368), (367, 369), (370, 372), (371, 373), (374, 376), (375, 378), (377, 380), (379, 381), (382, 384), (383, 386), (385, 388), (387, 389), (390, 391), (392, 393), (394, 396), (395, 397), (398, 399), (400, 401), (402, 403), (404, 405), (406, 407), (408, 410), (409, 412), (411, 413), (414, 416), (415, 417), (418, 420), (419, 422), (421, 366), (423, 425), (426, 428), (427, 429), (430, 431), (432, 433), (434, 436), (435, 438), (437, 440), (439, 441), (442, 443), (444, 446), (445, 447), (448, 450), (449, 451), (452, 454), (453, 455), (456, 458), (457, 459), (460, 462), (461, 464), (463, 466), (465, 467), (468, 469), (470, 471), (472, 473), (474, 476), (475, 478), (477, 479), (480, 482), (481, 424)
    ]))    
    
 