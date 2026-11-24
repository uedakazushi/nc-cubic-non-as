#!/usr/bin/env python3
"""Six-zero support exclusions through path length eight.
Finite-field ranks only skip impossible generic-excess witnesses. Every
exclusion is checked over the full characteristic-zero fraction field.
"""
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import os
import subprocess
from sympy import GF
from verify_ou_zero_supports import I, PAIRS, TRIPLES, normalize, support, path_dim

RESIDUAL = [(117358318,2),(117370606,2),(125759214,2),(131005679,2),
            (131037935,3),(133099246,2),(134145774,2)]
EXPECTED = {(2,4):51,(3,4):2,(2,5):76,(2,6):30,(2,7):23,(2,8):12,
            (2,0):6,(3,0):1}


def obstruction(mask, coeff, K, free):
    F = GF(1009)
    sup = set(support(mask))
    sample = {a:F.convert(2+free.index(a) if a in free else int(a in sup))
              for a in TRIPLES}
    for m in range(3,9):
        for r in I:
            for i,j in PAIRS:
                expected = {3:3+int(i==j),4:5,5:7,6:9+int(i==j),7:12,8:15}[m]
                if path_dim(sample,m,r,i,j,F) <= expected:
                    continue
                dim = path_dim(coeff,m,r,i,j,K)
                if dim > expected:
                    return (m,r,i,j,dim,expected)
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--first',type=int,default=0)
    p.add_argument('--last',type=int,default=201)
    args = p.parse_args()
    assert 0 <= args.first <= args.last <= 201
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'support_orbits'
        src = Path(__file__).with_name('ou_zero_supports.cpp')
        subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
        masks = [int(m) for z,m in map(str.split,subprocess.check_output([str(exe)],text=True).splitlines()) if z=='6']
    assert len(masks)==201
    counts,residual=Counter(),[]
    for mask in masks[args.first:args.last]:
        coeff,K,free=normalize(mask)
        w=obstruction(mask,coeff,K,free)
        counts[len(free),w[0] if w else 0]+=1
        if w is None:
            residual.append((mask,len(free)))
        print('CHECK',mask,len(free),w,flush=True)
    if args.first==0 and args.last==201:
        assert counts==EXPECTED and sorted(residual)==RESIDUAL
    print('COUNTS',sorted(counts.items()),'RESIDUAL',sorted(residual),flush=True)

if __name__=='__main__':
    main()
