#!/usr/bin/env python3
"""Numerical and determinant identities for the plane-cubic reduction.

The nonreduced classification and the full-support reconstruction on the
rank-two kernel locus are proved geometrically in the paper/notes. These
finite checks do not test AS regularity or the existence of a marking.
"""
from itertools import product
import sympy as sp


def nonreduced_degrees():
    # The degree of O_L(a) extended to 3L, and the weighted degree on 2L+M.
    a, b = sp.symbols('a b', integer=True)
    assert sp.expand(sum(a-j+1 for j in range(3))) == 3*a
    assert sp.expand((a+1)+a+(b+1)-2) == 2*a+b
    # Normalize the first block to (0,0); a_j is monotone and a_2 <= 1.
    candidates = []
    for a1, a2 in product(range(2), repeat=2):
        if a1 > a2:
            continue
        T = ((0,0), (a1,1-2*a1), (a2,2-2*a2))
        assert all(T[j][0]-T[i][0] >= 0 and T[j][1]-T[i][1] >= -1
                   for i in range(3) for j in range(i,3))
        # For a1 or a2 nonzero, one adjacent difference has simple-line degree -1.
        if a2:
            assert T[1][1]-T[0][1] == -1 or T[2][1]-T[1][1] == -1
            assert T[2] == (1,0)
        candidates.append(T)
    assert candidates == [((0,0),(0,1),(0,2)),
                          ((0,0),(0,1),(1,0)),
                          ((0,0),(1,-1),(1,0))]
    print('Nonreduced degrees: three candidates; two have a composition obstruction.')


def determinant_projection():
    a,b,c,d = sp.symbols('a b c d')
    t = sp.symbols('t0:6')
    C = sp.Matrix(3,3,sp.symbols('c0:9'))
    X = sp.Matrix([[1,1,1],[1,a,b],[1,c,d]])
    A,B,Cc,D,E,F = t
    coeff = (C[0,0]*C[1,1]*C[2,2], -C[0,0]*C[1,2]*C[2,1],
             -C[0,2]*C[1,1]*C[2,0], C[0,1]*C[1,2]*C[2,0],
             C[0,2]*C[1,0]*C[2,1], -C[0,1]*C[1,0]*C[2,2])
    z = (a*d,b*c,a,b,c,d)
    assert sp.expand(C.multiply_elementwise(X).det()-sum(x*y for x,y in zip(coeff,z))) == 0
    cubic = A*D*E-B*Cc*F
    assert sp.expand(cubic.subs(dict(zip(t,z)))) == 0
    # Inverse on the projective torus; all six t_i are invertible there.
    lam = A/(Cc*F)
    pre = [lam*Cc,lam*D,lam*E,lam*F]
    image = [pre[0]*pre[3],pre[1]*pre[2],*pre]
    diff = [sp.factor((p-lam*q)*Cc**2*F**2) for p,q in zip(image,t)]
    assert all(sp.rem(sp.Poly(x,B),sp.Poly(cubic,B)).is_zero for x in diff)
    print('Full-support determinant projection and torus inverse verified symbolically.')


if __name__ == '__main__':
    nonreduced_degrees()
    determinant_projection()
