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
#    spec Inferring Protein from Spectrum

from reference_tables import amino_acids
from bisect import bisect
from helpers import create_list

def create_lookup(amino_acids=amino_acids):
    pairs = sorted([(abbrev,value.mon_mass) for abbrev,value in amino_acids.items()],
                   key =lambda x:x[1])
    pairs.append(('?',999999999999))
    masses = [mass for (_,mass) in pairs]
    return masses,pairs
    
def spectrum2protein(ms):
    masses,pairs = create_lookup()

    diffs = [m1-m0 for m0,m1 in zip(ms[:-1],ms[1:])]
    def get_abbrev(diff):
        index = bisect(masses,diff)
        m1 = masses[index]
        m0 = masses[(index-1) if index>0 else 0]
        if index>0 and diff-m0 < m1-diff:
            index-=1
        abbrev,_ = pairs[index]

        return abbrev
    return ''.join([get_abbrev(diff) for diff in diffs])
  

if __name__=='__main__':

    pp=spectrum2protein([
        2355.76105643,
        2426.79816643,
        2563.85707643,
        2664.90475643,
        2792.99971643,
        2894.04739643,
        3025.08788643,
        3126.13556643,
        3255.17815643,
        3356.22583643,
        3443.25786643,
        3546.26705643,
        3617.30416643,
        3773.40527643,
        3910.46418643,
        4073.52751643,
        4229.62862643,
        4357.72358643,
        4513.82469643,
        4676.88802643,
        4839.95135643,
        4953.99428643,
        5025.03139643,
        5140.05833643,
        5268.11691643,
        5415.18532643,
        5562.25373643,
        5725.31706643,
        5824.38547643,
        5895.42258643,
        6081.50189643,
        6168.53392643,
        6281.61798643,
        6384.62717643,
        6481.67993643,
        6596.70687643,
        6759.77020643,
        6830.80731643,
        6958.90227643,
        7071.98633643,
        7129.00779643,
        7257.06637643,
        7314.08783643,
        7443.13042643,
        7571.22538643,
        7685.26831643,
        7786.31599643,
        7917.35648643,
        8016.42489643,
        8113.47765643,
        8216.48684643,
        8315.55525643,
        8414.62366643,
        8511.67642643,
        8612.72410643,
        8725.80816643,
        8853.90312643,
        8984.94361643,
        9099.97055643,
        9214.01348643,
        9343.05607643,
        9472.09866643,
        9571.16707643,
        9708.22598643,
        9836.28456643,
        9983.35297643,
        10096.4370364,
        10252.5381464,
        10323.5752564,
        10422.6436664,
        10479.6651264,
        10642.7284564,
        10828.8077664,
        10956.9027264,
        11084.9976864,
        11200.0246264,
        11337.0835364,
        11408.1206464,
        11545.1795564,
        11648.1887464,
        11795.2571564,
        11898.2663464,
        12001.2755364,
        12102.3232164,
        12233.3637064,
        12346.4477664,
        12532.5270764,
        12647.5540164,
        12748.6016964,
        12876.6966564,
        13004.7552364,
        13133.7978264,
        13262.8404164,
        13393.8809064,
        13556.9442364,
        13657.9919164,
        13745.0239464,
        13848.0331364,
        13905.0545964
       
])
    with open('spec.txt','w') as o:
        o.write('{0}\n'.format(pp))    