# Copyright (C) 2017 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# BA7C Implement Additive Phylogeny
from BA7B import ComputeLimbLength

class Tree(object):
    def __init__(self,N=-1):
        self.nodes=list(range(N))
        self.edges={}
    def link(self,start,end,weight=1): 
        print ('Linking {0} to {1} ({2})'.format(start,end,weight))
        self.half_link(start,end,weight)
        self.half_link(end,start,weight)
    def unlink(self,i,k):
        try:
            self.half_unlink(i,k)
            self.half_unlink(k,i)
        except KeyError:
            print ('Could not unlink {0} from {1}'.format(i,k))
            self.print()        
        
    def half_link(self,a,b,weight):
        if not a in self.nodes:
            self.nodes.append(a)        
        if a in self.edges:
            self.edges[a].append((b,weight))
        else:
            self.edges[a]=[(b,weight)]
    
    def half_unlink(self,a,b):
        links=[(e,w) for (e,w) in self.edges[a] if e != b]
        if len(links)<len(self.edges[a]):
            self.edges[a]=links
            print ('Unlinked {0} from {1}'.format(a,b))
        else:
            print ('Could not unlink {0} from {1}'.format(a,b))
            self.print()
    
    def are_linked(self,a,b):
        return len([e for (e,w) in self.edges[a] if e == b])>0
        
    def print(self,includeNodes=False):
        print('-----------------')
        self.nodes.sort()
        if includeNodes:
            print (self.nodes)
        for node in self.nodes:
            if node in self.edges:
                for edge in self.edges[node]:
                    end,weight=edge
                    print ('{0}->{1},{2}'.format(node,end,weight))
                    
    def next_node(self):
        return len(self.nodes)
    
    def traverse(self,i,k,depth=0):
        #print ('Traverse',i,k,depth)
        if not i in self.edges: return ([],[])
        path=[i]
        weights=[0]

        for j,w in self.edges[i]:
            if j in path: continue
            path.append(j)
            weights.append(w)
            if j==k:
                return list(zip(path,weights))
            else:
                if depth>25: return ([],[])      #FIXME
                return self.traverse(j,k,depth+1)
            

def DPrint(D):
    print ('=====================')
    for drow in D:
        print (', '.join([str(d) for d in drow]))
    
def AdditivePhylogeny(D,n,N=-1):
    def find_ikn(D):
        '''
        Find three leaves such that Di,k = Di,n + Dn,k
        '''
        ik=[]
        for i in range(n):
            for k in range(n):
                #print('i={0},k={1},D_ik={2},D_in={3},d_nk={4}'.format(i,k,D[i][k],D[i][n-1],D[n-1][k]))
                if D[i][k]==D[i][n-1]+D[n-1][k] and i!=k:
                    #print (i,k,n)
                    ik.append([i,k,n-1])
        return ik[0]
    
    if N==-1:
        N=n
        
    if n==2:
        T=Tree(N)
        T.link(0,1,D[0][1])
        return T
    else:
        limbLength=ComputeLimbLength(n,n-1,D)
        
        D_bald=[d_row[:] for d_row in D]
        for j in range(n-1):   
            D_bald[n-1][j]-=limbLength
            D_bald[j][n-1]=D_bald[n-1][j]  
         
        i,k,node=find_ikn(D_bald)
        x=D_bald[i][n-1]
        
        D_Trimmed=[D_bald[l][:-1] for l in range(n-1)]
        
        T=AdditivePhylogeny(D_Trimmed,n-1,N)
        print('i={0},n={1},k={2},x={3},limb length={4}'.format(i,n,k,x,limbLength))
        # v= the (potentially new) node in T at distance x from i on the path between i and k
        path,weights=T.traverse(i,k)
        v=None
        if len(path)>=2 and path[0]==i and path[-1]==k:
            print ("OK",pathe,weights)
        else:
            v=T.next_node()
            weight_i=ComputeLimbLength(n,i,D)
            weight_k=ComputeLimbLength(n,k,D)
            print ('New Node {0} between {1}({3}) and {2}({4})'.format(v,i,k,weight_i,weight_k))
            T.unlink(i,k)
            T.link(v,i,weight_i)
            T.link(v,k,weight_k)
        # add leaf n back to T by creating a limb (v, n) of length limbLength
        T.link(node,v,limbLength)
        T.print()
        return T
    
if __name__=='__main__':
    n=4
    D=[[0,   13,  21,  22],
        [13,  0 ,  12,  13],
        [21,  12 , 0 ,  13],
        [22 , 13 , 13,  0]]
    #for j in range(4):
        #print ('{0}: {1}'.format(j,ComputeLimbLength(n,j,D)))
    #for j in range(3):
        #print ('{0}: {1}'.format(j,ComputeLimbLength(3,j,D)))        
    #n=29
    #D=[
        #[0,3036,4777,1541,2766,6656,2401,4119,7488,4929,5344,3516,1485,6392,2066,3216,7008,7206,1187,6491,3379,6262,6153,4927,6670,4997,9010,5793,9032],
        #[3036,0,6323,3087,4312,8202,1619,5665,9034,2205,6890,966,3031,7938,3612,492,4284,8752,2023,8037,4925,3538,3429,6473,3946,2273,10556,7339,10578],
        #[4777,6323,0,4054,2571,3455,5688,1972,4287,8216,2143,6803,4126,3191,3183,6503,10295,4005,4474,3290,2964,9549,9440,1726,9957,8284,5809,2592,5831],
        #[1541,3087,4054,0,2043,5933,2452,3396,6765,4980,4621,3567,890,5669,1343,3267,7059,6483,1238,5768,2656,6313,6204,4204,6721,5048,8287,5070,8309],
        #[2766,4312,2571,2043,0,4450,3677,1913,5282,6205,3138,4792,2115,4186,1172,4492,8284,5000,2463,4285,1173,7538,7429,2721,7946,6273,6804,3587,6826],
        #[6656,8202,3455,5933,4450,0,7567,3851,1082,10095,1872,8682,6005,1992,5062,8382,12174,800,6353,1047,4843,11428,11319,3053,11836,10163,2604,2545,2626],
        #[2401,1619,5688,2452,3677,7567,0,5030,8399,3512,6255,2099,2396,7303,2977,1799,5591,8117,1388,7402,4290,4845,4736,5838,5253,3580,9921,6704,9943],
        #[4119,5665,1972,3396,1913,3851,5030,0,4683,7558,2539,6145,3468,3587,2525,5845,9637,4401,3816,3686,2306,8891,8782,2122,9299,7626,6205,2988,6227],
        #[7488,9034,4287,6765,5282,1082,8399,4683,0,10927,2704,9514,6837,2824,5894,9214,13006,1172,7185,1879,5675,12260,12151,3885,12668,10995,2144,3377,2166],
        #[4929,2205,8216,4980,6205,10095,3512,7558,10927,0,8783,2859,4924,9831,5505,1891,3719,10645,3916,9930,6818,2973,2864,8366,3381,1708,12449,9232,12471],
        #[5344,6890,2143,4621,3138,1872,6255,2539,2704,8783,0,7370,4693,1608,3750,7070,10862,2422,5041,1707,3531,10116,10007,1741,10524,8851,4226,1233,4248],
        #[3516,966,6803,3567,4792,8682,2099,6145,9514,2859,7370,0,3511,8418,4092,1146,4938,9232,2503,8517,5405,4192,4083,6953,4600,2927,11036,7819,11058],
        #[1485,3031,4126,890,2115,6005,2396,3468,6837,4924,4693,3511,0,5741,1415,3211,7003,6555,1182,5840,2728,6257,6148,4276,6665,4992,8359,5142,8381],
        #[6392,7938,3191,5669,4186,1992,7303,3587,2824,9831,1608,8418,5741,0,4798,8118,11910,2542,6089,1827,4579,11164,11055,2789,11572,9899,4346,2281,4368],
        #[2066,3612,3183,1343,1172,5062,2977,2525,5894,5505,3750,4092,1415,4798,0,3792,7584,5612,1763,4897,1785,6838,6729,3333,7246,5573,7416,4199,7438],
        #[3216,492,6503,3267,4492,8382,1799,5845,9214,1891,7070,1146,3211,8118,3792,0,3970,8932,2203,8217,5105,3224,3115,6653,3632,1959,10736,7519,10758],
        #[7008,4284,10295,7059,8284,12174,5591,9637,13006,3719,10862,4938,7003,11910,7584,3970,0,12724,5995,12009,8897,1442,2699,10445,1088,3447,14528,11311,14550],
        #[7206,8752,4005,6483,5000,800,8117,4401,1172,10645,2422,9232,6555,2542,5612,8932,12724,0,6903,1597,5393,11978,11869,3603,12386,10713,2694,3095,2716],
        #[1187,2023,4474,1238,2463,6353,1388,3816,7185,3916,5041,2503,1182,6089,1763,2203,5995,6903,0,6188,3076,5249,5140,4624,5657,3984,8707,5490,8729],
        #[6491,8037,3290,5768,4285,1047,7402,3686,1879,9930,1707,8517,5840,1827,4897,8217,12009,1597,6188,0,4678,11263,11154,2888,11671,9998,3401,2380,3423],
        #[3379,4925,2964,2656,1173,4843,4290,2306,5675,6818,3531,5405,2728,4579,1785,5105,8897,5393,3076,4678,0,8151,8042,3114,8559,6886,7197,3980,7219],
        #[6262,3538,9549,6313,7538,11428,4845,8891,12260,2973,10116,4192,6257,11164,6838,3224,1442,11978,5249,11263,8151,0,1953,9699,1104,2701,13782,10565,13804],
        #[6153,3429,9440,6204,7429,11319,4736,8782,12151,2864,10007,4083,6148,11055,6729,3115,2699,11869,5140,11154,8042,1953,0,9590,2361,2592,13673,10456,13695],
        #[4927,6473,1726,4204,2721,3053,5838,2122,3885,8366,1741,6953,4276,2789,3333,6653,10445,3603,4624,2888,3114,9699,9590,0,10107,8434,5407,2190,5429],
        #[6670,3946,9957,6721,7946,11836,5253,9299,12668,3381,10524,4600,6665,11572,7246,3632,1088,12386,5657,11671,8559,1104,2361,10107,0,3109,14190,10973,14212],
        #[4997,2273,8284,5048,6273,10163,3580,7626,10995,1708,8851,2927,4992,9899,5573,1959,3447,10713,3984,9998,6886,2701,2592,8434,3109,0,12517,9300,12539],
        #[9010,10556,5809,8287,6804,2604,9921,6205,2144,12449,4226,11036,8359,4346,7416,10736,14528,2694,8707,3401,7197,13782,13673,5407,14190,12517,0,4899,1758],
        #[5793,7339,2592,5070,3587,2545,6704,2988,3377,9232,1233,7819,5142,2281,4199,7519,11311,3095,5490,2380,3980,10565,10456,2190,10973,9300,4899,0,4921],
        #[9032,10578,5831,8309,6826,2626,9943,6227,2166,12471,4248,11058,8381,4368,7438,10758,14550,2716,8729,3423,7219,13804,13695,5429,14212,12539,1758,4921,0]]
    AdditivePhylogeny(D,n).print()
    
    #Output
    #0->55:745
    #1->48:156
    #2->52:788
    #3->54:409
    #4->53:280
    #5->49:125
    #6->51:492
    #7->50:657
    #8->31:311
    #9->41:820
    #10->47:280
    #11->46:723
    #12->45:417
    #13->44:864
    #14->43:236    
 
    def AdditivePhylogenyOld(D,n,T=None):
        
        def find_ik(n):
            '''
            Find three leaves such that Di,k = Di,n + Dn,k
            '''
            for i in range(n):
                for k in range(n):
                    print('i={0},k={1},D_ik={2},D_in={3},d_nk={4}'.format(i,k,D[i][k],D[i][n-1],D[n-1][k]))
                    if D[i][k]==D[i][n-1]+D[n-1][k]:
                        return (i,k)
            print ('not found')
         
        def  explore(i,k,path=[]):
            path.append(i)
            if i==k:
                return path
            else:
                for a,_ in T.edges[i]:
                    if a in path:
                        continue
                    test=explore(a,k,path)
                    if test!=None:
                        return test
            
        def Node(x,i,k):
            '''
            the (potentially new) node in T at distance x from i on the path between i and k
            '''
            print('Node: x={0},i={1},k={2}'.format(x,i,k))
            #for drow in D:
                #print (', '.join([str(d) for d in drow])) 
            print ('Path from i[{0}] to k[{1}] has length {2}'.format(i,k,D[i][k]))
            path=explore(i,k)
            print ("Path is",path)
            if path!=None:
                length=0
                for index in range(len(path)-1):
                    print (index,path[index],T.edges[path[index]],path[index+1])
                    for l,w  in T.edges[path[index]]:
                        print(l,w,path[index+1])
                        if path[index+1]==l:
                            length+=w
                            print(l,w,path[index+1],length)
                            if length==x:
                                print ("Found",l)
                                return l,False
                            if length>x:
                                print ("Overshot",l,length)
    
            return (T.next_node(),True)
        
        
        print ('\nAdditivePhylogeny n={0}'.format(n)) 
        if T==None:
            T=Tree(n)
       
        if n==2:
            T.link(0,1,D[0][1])
            #T.print(includeNodes=True)
            return T
        else:
            limbLength=ComputeLimbLength(n,n-1,D)
            for drow in D:
                print (', '.join([str(d) for d in drow]))        
            print('limbLength of {0} is {1}'.format(n-1,limbLength))
            for j in range(n-1):           #D_bald
                D[n-1][j]-=limbLength
                D[j][n-1]=D[n-1][j]
            i,k=find_ik(n)
            #print ('D{2},{0}[{3}]=D{0},{2}[{4}]+D{2},{1}[{5}]'.format(i,k,n-1,D[i][k], D[i][n-1], D[n-1][k])) # Di,k = Di,n + Dn,k
            x=D[i][n-1]
            #remove row n and column n from D
            D_Trimmed=[D[l][:-1] for l in range(n-1)]
            T=AdditivePhylogeny(D_Trimmed,n-1,T)
            print('i={0},n={1},k={2},x={3},limb length={4}'.format(i,n,k,x,limbLength))
            #for drow in D:
                #print (', '.join([str(d) for d in drow]))
            v,is_new=Node(x,i,k)
            print('Addnode {0}, at length {2} from {1} {3}'.format(n,v,limbLength,is_new))
    
            if is_new:
                T.unlink(i,k)
                T.link(i,v,ComputeLimbLength(n,i,D))
                T.link(k,v,ComputeLimbLength(n,k,D))
            T.link(n,v,limbLength)
            return T    