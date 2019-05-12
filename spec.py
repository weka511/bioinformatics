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
    pairs = sorted([(abbrev,1000*value.mon_mass) for abbrev,value in amino_acids.items()],
                   key =lambda x:x[1])
    masses = [mass for (_,mass) in pairs]
    return masses,pairs
    
def spectrum2protein(ms):
    masses,pairs = create_lookup()
    diffs = [1000*(m1-m0) for m0,m1 in zip(ms[:-1],ms[1:])]
    def get_abbrev(diff):
        index = bisect(masses,diff)
        m1 = masses[index]
        m0 = masses[(index-1) if index>0 else 0]
        if diff-m0 < m1-diff:
            index-=1
        abbrev,_ = pairs[index]
        return abbrev
    return ''.join([get_abbrev(diff) for diff in diffs])
  

if __name__=='__main__':
    #for key,value in amino_acids.items():
        #print (key,value.mon_mass)
    #print (spectrum2protein([3524.8542,
                     #3710.9335,
                     #3841.974,
                     #3970.0326,
                     #4057.0646]))
    pp=spectrum2protein([3721.5396671,
3834.6237271,
4020.7030371,
4134.7459671,
4263.7885571,
4391.8835171,
4504.9675771,
4652.0359871,
4723.0730971,
4851.1316771,
4938.1637071,
5052.2066371,
5165.2906971,
5236.3278071,
5364.3863871,
5451.4184171,
5554.4276071,
5667.5116671,
5780.5957271,
5936.6968371,
6007.7339471,
6094.7659771,
6195.8136571,
6298.8228471,
6385.8548771,
6522.9137871,
6637.9407271,
6694.9621871,
6796.0098671,
6943.0782771,
7000.0997371,
7057.1211971,
7170.2052571,
7283.2893171,
7397.3322471,
7500.3414371,
7663.4047671,
7750.4367971,
7821.4739071,
7949.5688671,
8048.6372771,
8234.7165871,
8337.7257771,
8468.7662671,
8624.8673771,
8780.9684871,
8894.0525471,
8995.1002271,
9052.1216871,
9165.2057471,
9351.2850571,
9537.3643671,
9640.3735571,
9755.4004971,
9869.4434271,
9984.4703671,
10121.5292771,
10178.5507371,
10235.5721971,
10363.6307771,
10434.6678871,
10531.7206471,
10717.7999571,
10831.8428871,
10962.8833771,
11065.8925671,
11136.9296771,
11250.0137371,
11436.0930471,
11592.1941571,
11778.2734671,
11941.3367971,
12078.3957071,
12193.4226471,
12306.5067071,
12492.5860171,
12605.6700771,
12718.7541371,
12815.8068971,
12952.8658071,
13108.9669171,
13222.0509771,
13337.0779171,
13450.1619771,
13564.2049071,
13727.2682371,
13828.3159171,
13942.3588471,
14079.4177571,
14208.4603471,
14321.5444071,
14434.6284671,
14491.6499271
])
    with open('spec.txt','w') as o:
        o.write('{0}\n'.format(pp))    