# Copyright (C) 2015 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.

# Exhaustively search for k-mers and compute statistics

import os, random, matplotlib.pyplot as plt

# Generate DNA at random
def random_dna(length):
    def choose_one():
        return 'ATGC'[random.choice(range(4))]
    return ''.join([choose_one() for i in range(length)])

# Count kmers that appear in dna

def generate_stats(dna,k):
    kmers={}
    for i in range(len(dna)-k):
        kmer=dna[i:i+k]
        if kmer in kmers:
            kmers[kmer]+=1
        else:
            kmers[kmer]=1
    return kmers

# Lookup table for reverse complements

dna_pairing={
 'T':'A',
 'A':'T',
 'C':'G',
 'G':'C'
}

# Computer reverse complement

def reverse_complement(string):
    return ''.join([dna_pairing[x] for x in string[::-1]])

# Find d-neighbours of kmer
# FIXME: the repeated joins are ugle
def neighbours(kmer,k,d):
    def join(xs):
        result=[]
        for x in xs:
            result=result+x
        return result
    # Generate a list of d-lists of indices that can be changed: e.g. for a d neighbourhood,
    # determine all the pairs that could be changed
    def indices(d):
        def incr(used):
            return [[i]+used for i in range(k) if i not in used]
        if d==0:
            return []
        elif d==1:
            return incr([])
        else:
            return join([incr(used) for used in indices(d-1)])
    # Generate all possible substitutions at a specific potion in a kmer
    def substitute(kmer,pos):
        def sub(letter):
            return kmer[0:pos] + letter + kmer[pos+1:]
        return [sub(letter) for letter in ['T','G','C','A']]
    # Generare alternative versions of a kmer given d-lists of indices
    def generate(iset):
        result=[kmer]
        for i in iset:
            temp=[]
            for rr in result:
                temp.append(substitute(rr,i))
            result=result+temp
        return result[1:]
    return join(join([generate(i) for i in indices(d)]))
   
# Given a table of kmers and counts, explore d-neighbourhood of each kmer and 
# update count to include counts of neighbours
def explore_neighbourhood(stats,k,d):
    if d==0:
        return stats    
    accumulated={}
    for key in stats.keys():
        accumulated[key]=sum([stats[kmer] for kmer in neighbours(key,k,d) if kmer in stats])    
    return accumulated

# include reverse complement in couts

def accumulate(stats):
    accumulated={}
    for key in stats.keys():
        rc=reverse_complement(key)
        if not rc in accumulated.keys():
            if rc!=key and rc in stats:
                accumulated[key]=stats[key]+stats[rc]
            else:
                accumulated[key]=stats[key]
    return accumulated

# Compute counts for kmers (inclusing reverse complements and neightbours)
# and plot them
def plot_stats(folder,name,k,d,m,fig):
    print (name,fig)
    with open (os.path.join(folder,name),'r') as file:
        dna=''.join(file.read().split('\n'))
        
        stats=accumulate(explore_neighbourhood(generate_stats(dna,k),k,d))
        counts=sorted(stats.values(),reverse=True)
        stats_random=accumulate(explore_neighbourhood(generate_stats(random_dna(len(dna)),k),k,d))
        counts_random=sorted(stats_random.values(),reverse=True)
        plt.figure(fig)
        plt.plot(counts,label=name)  
        plt.plot(counts_random,label="random")
        base=name.replace('.txt','')
        title = '%(k)d-mers in %(d)d-neighbourhood of %(base)s'%locals()
        plt.title(title)
        plt.legend(loc='upper right')
        with open("kmers.log","a") as out:
            out.write('%(title)s\n'%locals())
            for i in range(m):
                out.write('%(count)d,%(count_random)d,%(ratio)f\n'%{
                    'count': counts[i],
                    'count_random':counts_random[i],
                    'ratio' : counts[i]/max(counts_random[i],1)})
    plt.savefig('%(k)d-mers-%(d)d-neighbourhood-%(base)s.png'%locals())
   
if __name__=='__main__':
    fig=1
    #plot_stats('data','E-coli.txt',9,1,25,fig)
    for k in [7,8,9,10,11,13,17,19,29,39]:
        for name in [
            'E-coli.txt',
            'Salmonella_enterica.txt',
            'Vibrio_cholerae.txt'
            ]:
            plot_stats('data',name,k,0,25,fig)
            fig+=1
     
    plt.show()
