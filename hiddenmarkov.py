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