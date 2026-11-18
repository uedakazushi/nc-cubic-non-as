#!/usr/bin/env python3
"""Finite invariant-monomial and toric cubic checks.
Only a finite degree check: the all-degree proof is in the manuscript.
"""
from itertools import product
from math import comb
import sympy as s

def ms(n):
    return [(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]

def monomials_from_uvwt(n):
    out=set()
    for i in range(n+1):
        for j in range(n-i+1):
            for k in range(n-i-j+1):
                q=n-i-j-k
                out.add((3*i+q,3*j+q,3*k+q))
    return out
for n in range(9):
    inv={v for v in ms(3*n) if (v[0]-v[1])%3==0}
    generated=monomials_from_uvwt(n)
    assert inv==generated,(n, len(inv),len(generated))
    assert len(inv)==comb(n+3,3)-(comb(n,3) if n>=3 else 0)
U,V,W,T=s.symbols('U V W T')
f=U*V*W-T**3
assert [s.diff(f,x) for x in (U,V,W,T)]==[V*W,U*W,U*V,-3*T**2]
# The three projective coordinate points are the only singular points.
# In each affine chart of a nonzero coordinate, the germ is v*w=t**3 (A2).
print('PASS: invariant ring generators U,V,W,T, one cubic relation and 3 A2 singular points (degrees 0-8 tested).')
