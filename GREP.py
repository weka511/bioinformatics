#!/usr/bin/env python
'''
 Copyright (C) 2017 Greenweaves Software Pty Ltd

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

'''

import copy

from rosalind import dbru,read_strings

def count_kmers(S):
    '''
    Construct a list of counts for kmers, so we can make sure
    each cycle has the mataching frequencies.
    '''
    counts={}
    for s in S:
        if s in counts:
            counts[s]+=1
        else:
            counts[s]=1
    return counts

def create_lookup(B,E):
    '''
    Turn edges into a loolup table for successors
    '''
    F={}
    for b in B:
        F[b]=[f for (e,f) in E if e==b]
    return F

def remove_unused_kmer(counts):
    '''Get rid of zero counts'''
    removes=[]
    for key,value in counts.items():
        if value==0:
            removes.append(key)
    for key in removes:
        del counts[key]
    return counts

def format(r):
    return ''.join([rr[0] for rr in r] + [r[-1][-1]])

def grep(S):
    '''
    GREP Genome Assembly with Perfect Coverage and Repeats
    '''
    counts=count_kmers(S)
    B,E=dbru(S,include_revc=False)
    F=create_lookup(B,E)
    # We are going to build a list of runs, whicg are candidates for cycles.
    # Each time we have a choice of successor, we generate another run, so all
    # candidates are expanded successively.
    # If a candiate cannot be extended, jusr leave it in the list and remove
    # it at the end.
    Runs=[[S[0]]]
    counts[S[0]]-=1
    counts=remove_unused_kmer(counts)

    # We maintain a separate list of counts for each run, so we
    # know which symbols are available.

    CountsForRuns=[counts]

    for n in range(len(S)-1): # Grow runs by one base each iteration

        NewRuns=[]  # Stick extra runs on the end when we are done

        for i in range(len(Runs)):
            run=Runs[i]
            counts=CountsForRuns[i]
            last=run[-1][1:]
            succ=F[last]
            counts_old=copy.deepcopy(counts)
            j=0
            added=False
            while j<len(succ):
                kmer=last+succ[j][-1]
                if kmer in counts_old:
                    if added:
                        new_counts=copy.deepcopy(counts_old)
                        new_counts[kmer]-=1
                        new_counts=remove_unused_kmer(new_counts)
                        new_run=copy.deepcopy(run[:-1])
                        new_run.append(kmer)
                        CountsForRuns.append(new_counts)
                        NewRuns.append(new_run)
                    else:
                        counts[kmer]-=1
                        counts=remove_unused_kmer(counts)
                        run.append(kmer)
                        added=True
                j+=1
        Runs = Runs + NewRuns

    return [format(r)[:-1] for r in Runs if len(r)==len(S)]

if __name__=='__main__':

    for s in grep(read_strings('c:/Users/Weka/Downloads/rosalind_grep(1).txt')):
        print (s)
