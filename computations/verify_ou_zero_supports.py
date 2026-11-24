#!/usr/bin/env python3
"""Enumerate support orbits and exclude the large-zero strata exactly.
Requires a C++17 compiler and SymPy. No probabilistic rank tests.
"""
from collections import Counter
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import subprocess
from sympy import Matrix, QQ, symbols
from sympy.polys.matrices import DomainMatrix

I = range(3)
TRIPLES = list(product(I, repeat=3))
PAIRS = list(product(I, repeat=2))


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


def obstruction(coeff, K):
    for length in (3,4,5):
        for rot in I:
            for i,j in PAIRS:
                expected = {3:3+int(i==j),4:5,5:7}[length]
                dim = path_dim(coeff,length,rot,i,j,K)
                if dim > expected:
                    return (length,rot,i,j,dim,expected)
    return None


def main():
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'support_orbits'
        src = Path(__file__).with_name('ou_zero_supports.cpp')
        subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
        data = [tuple(map(int,line.split())) for line in subprocess.check_output([str(exe)],text=True).splitlines()]
    assert Counter(n for n,_ in data) == {9:217,10:53,11:4,12:1}
    counts = Counter()
    residual = []
    for n,mask in data:
        sup = support(mask)
        _, pivots = incidence(sup).T.rref()
        free = [t for i,t in enumerate(sup) if i not in pivots]
        # A full-row-rank integer exponent matrix gives a surjective torus
        # map over C, so pivot coefficients can all be set equal to one.
        variables = symbols('u:'+str(len(free)))
        K = QQ.frac_field(*variables) if variables else QQ
        coeff = {t: K.one if t in sup else K.zero for t in TRIPLES}
        for t,v in zip(free,variables):
            coeff[t] = K.from_sympy(v)
        witness = obstruction(coeff,K)
        counts[n,len(free),witness[0] if witness else 0] += 1
        if witness is None:
            residual.append((n,mask,free))
    assert residual == [(9,108455662,[(2,2,2)])]
    assert counts == {(9,0,4):201,(9,0,5):10,(9,1,4):5,(9,1,0):1,
                      (10,0,4):52,(10,0,5):1,(11,0,4):4,(12,0,4):1}, counts
    # The sole remaining support is i-j+k != 0 mod 3; j -> -j gives i+j+k != 0.
    assert set(support(residual[0][1])) == {t for t in TRIPLES if (t[0]-t[1]+t[2])%3}
    print('PASS: exact support-orbit and path-rank classification:', sorted(counts.items()), flush=True)
    print('PASS: sole large-zero survivor has Latin-square zero support and one parameter')

if __name__ == '__main__':
    main()
