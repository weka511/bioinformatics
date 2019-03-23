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
#    Hidden Markov Models

# BA10A Compute the Probability of a Hidden Path

from numpy import argmax

def ProbabilityHiddenPath(path,Transition):
    result = 0.5
    for i in range(1,len(path)):
        result *= Transition[(path[i-1],path[i])]
    return result

# BA10B Compute the Probability of an Outcome Given a Hidden Path 

def ProbabilityOutcomeGivenHiddenPath(string,path,Emission):
    result = 1
    for i in range(0,len(path)):
        result *= Emission[(path[i],string[i])]
    return result 

# BA10C Implement the Viterbi Algorithm 

def Viterbi(xs,alphabet,States,Transition,Emission):
    
    def calculateproduct_weights(s_source = 1):
        def product_weight(k,x):
            return max([s[-1][l] * Transition[(States[l],k)] * Emission[(k,x)] for l in range(len(States))])
        
        s = []  # first index is position, 2nd state
    
        s.append([s_source * (1/len(States)) * Emission[(k,xs[0])] for k in States])
        for x in xs[1:]:
            s.append([product_weight(k,x) for k in States])
        return s
    
     
    def backtrack(s):
        n     = len(s) - 1
        state = argmax(s[n])
        path  = [States[state]]
        while True:
            ps = [s[n-1][l] * Transition[(States[l],States[state])]  for l in range(len(States))]
            state = argmax(ps)
            path.append(States[state])
            n-=1
            if n<=0: return path[::-1]
        #for step in s[::-1]:
            #state = States[0]
            #best  = step[0]
            #for i in range(1,len(States)):
                #if step[i]>best:
                    #state = States[i]
                    #best  = step[i]
            #path.append(state)
        
    return ''.join(backtrack(calculateproduct_weights()))