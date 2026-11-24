#!/usr/bin/env python3
"""Exact finite certificate for the six-sheet formal lifting at 116309742, u=1.

The invertible Jacobian minor gives a solution over Q[[q]] by recursive
coefficient lifting. This script does not approximate the formal solution.
"""
from collections import Counter
from itertools import product
from sympy import Integer, Matrix, Symbol, expand

MASK = 116309742
P = [(0,0,0),(1,1,2),(2,2,0),(0,1,1),(2,0,2),(1,2,1)]
I = range(3)


def coeff(i,j,k):
    return (MASK >> (9*i+3*j+k)) & 1


def c(a,b,d):
    return coeff(P[a%6][0],P[b%6][1],P[d%6][2])


def main():
    for i,j,k in P:
        assert coeff(i,j,k)==0
        assert all(coeff(i,j,x) for x in I if x!=k)
        assert all(coeff(i,x,k) for x in I if x!=j)
        assert all(coeff(x,j,k) for x in I if x!=i)
    assert Counter((d,v[d]) for v in P for d in I)=={(d,i):2 for d,i in product(I,repeat=2)}
    assert all(all(x!=y for x,y in zip(P[a],P[(a+1)%6])) for a in range(6))
    assert all(c(a,b,d)==c(a+3,b+3,d+3) for a,b,d in product(range(6),repeat=3))
    q=Symbol('q')
    A={(a,d):Symbol(f'a{a}{d}') if d else Integer(1) for a,d in product(I,range(6))}
    B={(a,d):Symbol(f'b{a}{d}') if d else Integer(1) for a,d in product(I,range(6))}
    variables=[v for v in (*A.values(),*B.values()) if v!=1]
    equations=[]
    for a,d in product(I,range(6)):
        E=expand(sum(c(a,a+p,a+d)*A[a,p]*B[(a+p)%3,(d-p)%6]
                     *q**((p+(d-p)%6-d)//6) for p in range(6)))
        if d==0:
            assert E.subs(q,0)==0
            E=expand(E/q)
        equations.append(E)
    A0=[[1,1,1],[1,1,1],[1,0,0],[0,-1,0],[0,0,-1],[0,0,0]]
    B0=[[1,1,1],[-1,-1,-1],[0,0,1],[1,0,0],[0,-1,0],[0,0,0]]
    base={**{A[a,d]:A0[d][a] for a,d in product(I,range(1,6))},
          **{B[a,d]:B0[d][a] for a,d in product(I,range(1,6))},q:0}
    assert all(E.subs(base)==0 for E in equations)
    J=Matrix(equations).jacobian(variables).subs(base)
    pivots=[0,1,2,3,5,6,7,8,10,11,12,13,15,16,17,18,20,21]
    assert J.shape==(18,30)
    assert J[:,pivots].det()==-8
    # In the forward representation each row a has the order-one entry
    # A[a,a+1], and the other sheet b over i_a is different from a+1.
    # The two source sheets b,a+1 lie outside j_a; the selected minor has
    # one unit entry and leading term +/-t. For reversal B^t, use a-1.
    for axis in (0,2):
        step=1 if axis==0 else -1
        for a in range(6):
            b,= [b for b in range(6) if b!=a and P[b][axis]==P[a][axis]]
            h=(a+step)%6
            assert h!=b and P[b][1]!=P[a][1] and P[h][1]!=P[a][1]
            # Leading coefficient in A or B^t is respectively +1 or -1.
            assert (A0 if axis==0 else B0)[1][a%3 if axis==0 else h%3] in (1,-1)
    print('PASS: isolated zeros, double coverage, 18 formal equations, Jacobian minor -8, bounded heads in both orientations')

if __name__=='__main__':
    main()
