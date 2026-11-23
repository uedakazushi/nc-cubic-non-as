#!/usr/bin/env python3
"""Exact comparison of residual IV-type directed algebras.

A section of a line bundle of degree vector d is a tuple of polynomials
on the three normalizations, with equal constant terms and sum of linear
coefficients equal to lambda times the common constant term.
All assertions are polynomial identities over Q, with only the indicated
nonzero gluing differences required for the rank certificates.
"""
from collections import Counter
from itertools import combinations, permutations, product
import sympy as sp
from verify_ou_reduced_multidegrees import patterns, invariants, normalize, sub, forward
from verify_ou_reduced_algebraic_moves import terminal, edges, unit

u, v, w = (sp.symbols(f'{b}0:3') for b in 'uvw')
params = (u, v, w)
z = (0, 0, 0)
A = ((z,)*3, ((-1,1,1),(0,0,1),(0,1,0)), ((-1,1,2),(-1,2,1),(0,1,1)))
B = ((z,z,(0,1,-1)), ((-1,1,1),(0,1,0),(0,1,0)), ((-1,2,1),(-1,2,1),(0,1,1)))
C = ((z,z,(0,1,-1)), ((0,1,0),(1,1,-1),(1,1,-1)), ((1,2,-1),(1,2,-1),(2,1,-1)))
D = ((z,(1,-1,0),(1,0,-1)), ((0,0,1),(0,1,0),(1,0,0)), ((0,1,1),)*3)
E = ((z,(0,1,-1),(1,0,-1)), ((0,0,1),(0,1,0),(1,0,0)), ((0,0,2),(0,1,1),(1,0,1)))
F = ((z,(1,-1,0),(1,0,-1)), ((0,0,1),(1,-1,1),(1,0,0)), ((0,0,2),(1,-1,2),(1,0,1)))
A2 = ((z,)*3, ((0,0,1),(0,1,0),(1,0,0)), ((0,1,1),(1,0,1),(1,1,0)))
B2 = ((z,z,(0,1,-1)), ((1,0,0),(0,1,0),(0,1,0)), ((1,1,0),(1,1,0),(0,1,1)))
D2 = ((z,(0,1,-1),(1,0,-1)), ((0,1,0),(1,0,0),(1,1,-1)), ((1,1,0),)*3)


def residual():
    ts = [T for T in patterns(3) if not any(invariants(*T)[0])]
    known = {T for T in ts if terminal(T)}
    adj = {T: [S for S, _ in edges(T)] for T in ts}
    while True:
        new = known | {T for T in ts if any(S in known for S in adj[T])}
        if new == known:
            break
        known = new
    return [T for T in ts if T not in known]


def orbit(T):
    return min(normalize(*(tuple(tuple(x[i] for i in p) for x in block)
                           for block in T)) for p in permutations(range(3)))


def section(d, lam):
    """One nonzero degree-one section, using the conductor condition."""
    assert sum(d) == 1
    if unit(d):
        return tuple((1, lam) if a else (1, 0) for a in d)
    if sorted(d) == [-1, 1, 1]:
        i, j = [i for i, a in enumerate(d) if a == 1]
        return tuple((0, 1) if k == i else (0, -1) if k == j else (0, 0)
                     for k in range(3))
    if sorted(d) == [-1, 0, 2]:
        return tuple((0, 0, 1) if a == 2 else (0, 0) for a in d)
    raise ValueError(d)


def product_sections(a, b):
    out = []
    for f, g in zip(a, b):
        h = [0]*4
        for i, x in enumerate(f):
            for j, y in enumerate(g):
                h[i+j] += x*y
        out.extend(h)
    return out


def guaranteed_nonzero(f, T):
    if f in (1, -1):
        return True
    for block, lam in zip(T, params):
        for i, j in combinations(range(3), 2):
            if block[i] == block[j] and f in (lam[i]-lam[j], lam[j]-lam[i]):
                return True
    return False


def relation(T, i, k):
    U, V, W = T
    P = sp.Matrix.hstack(*(sp.Matrix(product_sections(
        section(sub(V[j], U[i]), v[j]-u[i]),
        section(sub(W[k], V[j]), w[k]-v[j]))) for j in range(3)))
    # A minor equal to +/-1 or a required nonzero gluing difference
    # ensures the two rows stay independent on the entire admissible locus.
    for a, b in combinations(range(P.rows), 2):
        r = tuple(sp.factor(x) for x in P[a, :].cross(P[b, :]))
        if any(guaranteed_nonzero(x, T) for x in r):
            assert all(sp.expand(x) == 0 for x in P*sp.Matrix(r))
            return r
    raise AssertionError(('no universal rank-two minor', T, i, k))


def proportional(a, b):
    return all(sp.expand(a[j]*b[l]-a[l]*b[j]) == 0
               for j, l in combinations(range(3), 2))


def check_pattern(T):
    assert all(sum(x) == b for b, block in enumerate(T) for x in block)
    assert all(forward(x, y) for b, block in enumerate(T) for x in block
               for later in T[b:] for y in later)
    assert not any(invariants(*T)[0])


def match(X, Y, uv=(), vw=()):
    check_pattern(X)
    check_pattern(Y)
    # All necessary gluing inequalities for Y follow from those for X.
    assert all((X[b][i] == X[b][j]) == (Y[b][i] == Y[b][j])
               for b in range(3) for i, j in combinations(range(3), 2))
    for i, k in product(range(3), repeat=2):
        a = relation(X, i, k)
        b = relation(Y, i, k)
        a = tuple(x*(-1)**(int((i,j) in uv)+int((j,k) in vw))
                  for j, x in enumerate(a))
        assert proportional(a, b), (i, k, a, b)


def main():
    assert Counter(orbit(T) for T in residual()) == Counter({A:3,B:6,C:6,D:3,E:3,F:3})
    for X, Y, uv, vw in ((A,A2,(),((2,0),)), (B,B2,(),()),
                         (C,B2,((2,0),),((0,2),)), (D,D2,((2,0),),())):
        assert invariants(*Y)[1] == (1,1,1)
        match(X, Y, uv, vw)
    match(E, F)
    for i, k in product(range(3), repeat=2):
        if i == k:
            target = tuple(0 if j == i else 1 if j == (i+1)%3 else -1 for j in range(3))
        else:
            target = tuple(1 if j == i else -1 if j == k else -u[i]+v[i]+v[k]-w[k]
                           for j in range(3))
        assert proportional(relation(E,i,k), target)
    print('IV: 24 residual candidates in six component-permutation orbits;')
    print('18 candidates match positive w=(1,1,1) algebras for all allowed parameters;')
    print('6 candidates give the same explicit five-parameter relation family.')


if __name__ == '__main__':
    main()
