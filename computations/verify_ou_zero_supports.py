#!/usr/bin/env python3
"""Exact support-orbit and path-dimension exclusions for z=7,...,12.
Requires C++17 and SymPy. Rational samples only skip impossible universal
rank-excess witnesses; every exclusion is verified over the full field.
"""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import subprocess
from sympy import Matrix, QQ, symbols
from sympy.polys.matrices import DomainMatrix

I = range(3)
TRIPLES = list(product(I, repeat=3))
PAIRS = list(product(I, repeat=2))
EXPECTED = {7:344, 8:387, 9:217, 10:53, 11:4, 12:1}
RESIDUAL7 = [116309742,116834030,117368558,125234926,125531767,
             128957047,129972855,133118583]


def support(mask):
    return [t for n,t in enumerate(TRIPLES) if mask >> n & 1]


def incidence(sup):
    return Matrix([[int(ij == (i,j)) for ij in PAIRS]
                   + [int(jk == (j,k)) for jk in PAIRS]
                   + [int(ki == (k,i)) for ki in PAIRS] for i,j,k in sup])


def path_dim(coeff, length, rotation, i, j, K):
    words = list(product(I, repeat=length-1))
    index = {w:n for n,w in enumerate(words)}
    rows = []
    for pos in range(length-1):
        r = (rotation+pos) % 3
        for w in words:
            if w[pos]:
                continue
            row = {}
            for a in I:
                v = w[:pos]+(a,)+w[pos+1:]
                full = (i,)+v+(j,)
                t = full[pos:pos+3]
                q = coeff[t[-r:]+t[:-r]]
                if q:
                    row[index[v]] = q
            if row:
                rows.append(row)
    M = DomainMatrix(dict(enumerate(rows)), (len(rows),len(words)), K)
    return len(words)-M.rank()


def normalize(mask):
    sup = support(mask)
    _, pivots = incidence(sup).T.rref()
    free = [t for i,t in enumerate(sup) if i not in pivots]
    variables = symbols('u:'+str(len(free)))
    K = QQ.frac_field(*variables) if variables else QQ
    coeff = {t: K.one if t in sup else K.zero for t in TRIPLES}
    for t,v in zip(free,variables):
        coeff[t] = K.from_sympy(v)
    return coeff,K,free


def obstruction(mask, coeff, K, free, lengths):
    # Any rational specialization has rank <= generic rank. If its quotient
    # dimension is <= the Euler requirement, generic excess is impossible.
    # If it exceeds the requirement, test the original fraction-field matrix.
    sample_values = dict(zip(free,range(2,2+len(free))))
    sup = set(support(mask))
    sample = {t:QQ(sample_values.get(t,1) if t in sup else 0) for t in TRIPLES}
    for length in lengths:
        for rot in I:
            for i,j in PAIRS:
                expected = {3:3+int(i==j),4:5,5:7,6:9+int(i==j),7:12}[length]
                d0 = path_dim(sample,length,rot,i,j,QQ)
                if d0 <= expected:
                    continue
                dim = path_dim(coeff,length,rot,i,j,K) if free else d0
                if dim > expected:
                    return (length,rot,i,j,dim,expected)
    return None


def periodic_support(mask):
    # Exact support test, allowing independent relabellings of each axis.
    T = set(support(mask))
    repeated = {t for t in TRIPLES if len(set(t)) == 2}
    diag = {(i,i,i) for i in I}
    cycles = [{(i,(i+s)%3,(i+2*s)%3) for i in I} for s in (1,2)]
    for a,b,c in product(permutations(I),repeat=3):
        U = {(a[i],b[j],c[k]) for i,j,k in T}
        if repeated <= U and not diag & U and any(not cyc & U for cyc in cycles):
            return True
    return False


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--zeros', type=int, nargs='+', default=list(EXPECTED))
    args = p.parse_args()
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'support_orbits'
        src = Path(__file__).with_name('ou_zero_supports.cpp')
        subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
        data = [tuple(map(int,line.split())) for line in subprocess.check_output([str(exe)],text=True).splitlines()]
    for z in args.zeros:
        masks = [m for n,m in data if n == z]
        assert len(masks) == EXPECTED[z]
        counts, residual = Counter(), []
        for mask in masks:
            coeff,K,free = normalize(mask)
            lengths = (3,4,5,6,7) if z in (7,8) else (3,4,5)
            w = obstruction(mask,coeff,K,free,lengths)
            counts[len(free),w[0] if w else 0] += 1
            if not w:
                residual.append((mask,len(free)))
        expected_counts = {
            7:{(1,4):177,(2,4):15,(1,5):102,(2,5):1,(1,6):26,(1,7):15,(1,0):7,(2,0):1},
            8:{(0,4):237,(1,4):82,(0,5):49,(1,5):6,(0,6):10,(0,7):2,(1,0):1},
            9:{(0,4):201,(0,5):10,(1,4):5,(1,0):1},
            10:{(0,4):52,(0,5):1},11:{(0,4):4},12:{(0,4):1},
        }
        assert counts == expected_counts[z], (z,counts)
        expected_residual = {7:[(m,2 if m==125234926 else 1) for m in RESIDUAL7],
                             8:[(108979950,1)],9:[(108455662,1)]}
        assert sorted(residual) == expected_residual.get(z,[]), (z,residual)
        periodic = [m for m,_ in sorted(residual) if periodic_support(m)]
        assert periodic == {7:[125234926],8:[108979950],9:[108455662]}.get(z,[])
        print('PASS',z,'counts',sorted(counts.items()),'residual',sorted(residual),'periodic',periodic,flush=True)
    # The old and new unique residual supports differ by one extra coefficient.
    assert set(support(108455662)) == {t for t in TRIPLES if (t[0]-t[1]+t[2])%3}
    assert set(support(108979950))-set(support(108455662)) == {(2,0,1)}
    assert set(support(108455662)) < set(support(108979950))

if __name__ == '__main__':
    main()
