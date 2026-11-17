#!/usr/bin/env python3
"""Finite character and sector verification for the non-CY sign-type cubic.

Does not establish the graded Morita equivalence or orbifold HKR theorem.
"""
import json
from itertools import product
from pathlib import Path

H = [(a,b,(a+b)%2) for a,b in product(range(2),repeat=2)]
e=(0,0,0); A=(0,1,1); B=(1,0,1); C=(1,1,0)
assert set(H)=={e,A,B,C}

def torsion(g,h):
    return (g[0]*h[1]+g[1]*h[0]) % 2

assert all(torsion(g,g)==0 for g in H)
assert all(any(torsion(g,h) for h in H) for g in H if g!=e)
assert all(torsion(B,h)==h[1] and torsion(C,h)==h[2] for h in H)
euler={e:9,A:5,B:1,C:1}
trace={g:euler[g]-2 for g in H}
assert [trace[g] for g in (e,A,B,C)]==[7,3,-1,-1]
rank=sum(trace.values())//len(H)
assert rank==2
counts={'-1':2,'0':rank+2,'1':2}
assert counts=={'-1':2,'0':4,'1':2}
sectors=[{'g':list(g),'fixed':('surface' if g==e else 'line_and_3_points' if g==A else 'elliptic_and_point'),
          'euler_characteristic':euler[g], 'trace_H2':trace[g],
          'HH_minus1':int(g in (B,C)), 'HH_0':4 if g==e else 0,
          'HH_1':int(g in (B,C))} for g in (e,A,B,C)]
cert={'group':'(Z/2)^2', 'torsion_exponent':'a_g*b_h+b_g*a_h',
      'invariant_picard_rank':rank,'sectors':sectors,'hh_dimensions':counts}
Path(__file__).with_name('sign_family_sector_certificate.json').write_text(
    json.dumps(cert,indent=2,ensure_ascii=False)+'\n')
print('PASS',counts)
