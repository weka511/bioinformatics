# Copyright (C) 2019 Greenweaves Software Limited

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

# BA10C Implement the Viterbi Algorithm 

from hmm import Viterbi

if __name__=='__main__':
    from helpers import create_strings
    
    strings    = create_strings(ext=3)
    xs         = strings[0]
    alphabet   = strings[2].split()
    States     = strings[4].split()
    Transition = {}
    i          = 0
    while i<len(States):
        items = strings[7+i].split()
        for j in range(1,len(items)):
            Transition[(items[0],States[j-1])] = float(items[j])
        i+=1
    i+=9
    
    Emission   = {}
    while i<len(strings):
        items = strings[i].split()
        for j in range(1,len(items)):        
            Emission[(items[0],alphabet[j-1])] = float(items[j])
        i+=1  
    print ( Viterbi(xs,alphabet,States,Transition,Emission))
    print ()
    Transition = {
        ('A','A') : 0.641, ('A','B') :   0.359,
        ('B','A') : 0.729, ('B','B') :   0.271,
    }
    Emission = {
        ('A','x') : 0.117, ('A','y') :   0.691, ('A','z') : 0.192,
        ('B','x') : 0.097, ('B','y') :   0.42,  ('B','z') : 0.483,
    }    
    print (Viterbi('xyxzzxyxyy','xyz','AB',Transition,Emission))
    
    Transition1 = {
        ('A','A') : 0.634, ('A','B') :   0.366,
        ('B','A') : 0.387, ('B','B') :   0.613,
    }
    Emission1 = {
        ('A','x') : 0.532, ('A','y') :   0.226, ('A','z') : 0.241,
        ('B','x') : 0.457, ('B','y') :   0.192, ('B','z') : 0.351,
    }    
    print (Viterbi('zxxxxyzzxyxyxyzxzzxzzzyzzxxxzxxyyyzxyxzyxyxyzyyyyzzyyyyzzxzxzyzzzzyxzxxxyxxxxyyzyyzyyyxzzzzyzxyzzyyy',
                   'xyz','AB',Transition1,Emission1))
    # AAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBAAA
    # AAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBAAA

    
    

   
    
