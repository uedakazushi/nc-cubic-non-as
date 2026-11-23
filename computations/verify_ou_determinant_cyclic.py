#!/usr/bin/env python3
"""Finite identities for cyclic syzygies and determinant-plane recognition.

The smoothness and reconstruction theorems are geometric proofs, not
consequences of the illustrative tensor tests at the end of this script.
"""
from itertools import combinations, permutations, product
import sympy as sp


P = tuple(permutations(range(3)))
E = tuple(p for p in P if sum(p[i] > p[j] for i, j in combinations(range(3), 2)) % 2 == 0)
O = tuple(p for p in P if p not in E)


def matching_checks():
    assert len(E) == len(O) == 3
    expected = {frozenset((i, j) for j in range(3)) for i in range(3)} | {
                frozenset((i, j) for i in range(3)) for j in range(3)}
    got = set()
    for os in permutations(O):
        shared = []
        for p, q in zip(E, os):
            common = [(i, p[i]) for i in range(3) if p[i] == q[i]]
            assert len(common) == 1
            shared += common
        got.add(frozenset(shared))
    assert got == expected
    print('All six even/odd matchings are rows or columns.')


def euler_checks():
    I, J, Z = sp.eye(3), sp.ones(3), sp.zeros(3)
    D = I.row_join(J).row_join(2*J).col_join(Z.row_join(I).row_join(J)).col_join(Z.row_join(Z).row_join(I))
    S = D.inv()*D.T
    assert S.T*D == (I+3*J).row_join(5*J).row_join(7*J).col_join(
        (2*J).row_join(I+3*J).row_join(5*J)).col_join(J.row_join(2*J).row_join(I+3*J))
    print('Cyclic two-step and three-step Euler dimensions verified.')


def classical_checks():
    v = [sp.eye(7)[:, i] for i in range(7)]
    h, es = v[0], v[1:]
    zero = sp.zeros(7, 1)
    sumv = lambda ls: sum(ls, zero)
    K = 3*h-sumv(es)
    U = (zero, -2*h+sumv(es), -h+sumv(es[:3]))
    V = tuple(es[:3])
    W = tuple(h-e for e in es[3:])
    nextU = tuple(u+K for u in U)
    lines = {tuple(e) for e in es} | {tuple(h-a-b) for a,b in combinations(es,2)} | {
        tuple(2*h-sumv(es)+e) for e in es}
    assert len(lines) == 27
    for X,Y in ((U,V),(V,W),(W,nextU)):
        ds = {tuple(y-x) for x in X for y in Y}
        assert len(ds) == 9 and ds <= lines
        assert sumv(Y)-sumv(X) == K
    print('All three cyclic classical grids consist of nine distinct lines.')


def dual_checks():
    z = sp.symbols('x1:4 y1:4')
    # Independent torus coordinates subject to the binomial identity.
    X = z[:3]
    Y = z[3:5]+(sp.prod(X)/(z[3]*z[4]),)
    P0 = sp.prod(X)
    lam = [P0/x for x in X]+[-P0/y for y in Y]
    assert sp.cancel(sp.prod(lam[:3])+sp.prod(lam[3:])) == 0
    for start, block in ((0, X), (3, Y)):
        for j in range(3):
            assert sp.cancel(sp.prod(lam[start+i] for i in range(3) if i != j)-P0*block[j]) == 0
    print('Torus gradient duality verified.')


def permutation_products(C):
    return sp.Matrix([sp.prod(C[i,p[i]] for i in range(3)) for p in E] +
                     [-sp.prod(C[i,p[i]] for i in range(3)) for p in O])


def inverse_products(v):
    # Undo signs and use the normalized representative X from the paper.
    terms = {tuple(p): v[j] for j,p in enumerate(E)} | {
             tuple(p): -v[3+j] for j,p in enumerate(O)}
    z0,z1,z2,z3,z4,z5 = (terms[p] for p in ((0,1,2),(0,2,1),(2,1,0),(1,2,0),(2,0,1),(1,0,2)))
    a,b,c,d = [sp.cancel(z0*t/(z2*z5)) for t in (z2,z3,z4,z5)]
    C = sp.Matrix([[1,1,1],[1,a,b],[1,c,d]])
    w = permutation_products(C)
    assert all(sp.cancel(w[i]*v[0]-w[0]*v[i]) == 0 for i in range(6))
    return C


def syzygy_checks(Cs):
    # Tensor entries c_ijk; each C_k has rows i and columns j.
    for axis in range(3):
        rest = [i for i in range(3) if i != axis]
        slices = []
        for a in range(3):
            def entry(b,c):
                ijk = [0]*3
                ijk[axis], ijk[rest[0]], ijk[rest[1]] = a,b,c
                return Cs[ijk[2]][ijk[0],ijk[1]]
            slices.append(sp.Matrix(3,3,entry))
        for C,D in combinations(slices,2):
            R = C.multiply_elementwise(D.applyfunc(lambda t: 1/t))
            assert R.rank() >= 2


def coefficient_examples():
    # A secant of the dual binomial cubic, with its third intersection.
    C0 = sp.Matrix([[1,1,1],[1,2,3],[1,5,7]])
    C1 = sp.Matrix([[1,1,1],[1,11,13],[1,17,19]])
    v0,v1 = map(permutation_products,(C0,C1))
    t,s = sp.symbols('t s')
    L = t*v0+s*v1
    f = sp.factor(sp.prod(L[:3])+sp.prod(L[3:]))
    roots = sp.solve(f.subs(t,1),s)
    q = next(q for q in roots if q != 0)
    v2 = L.subs({t:1,s:q})
    assert all(v2) and sp.prod(v2[:3])+sp.prod(v2[3:]) == 0
    assert sp.Poly(f.subs(t,1),s).degree() == 2
    assert sp.gcd(sp.Poly(f.subs(t,1),s),sp.Poly(sp.diff(f.subs(t,1),s),s)).degree() == 0
    C2 = inverse_products(v2)
    Cs = [C0,C1,C2]
    assert sp.Matrix.hstack(*(permutation_products(C) for C in Cs)).rank() == 2
    syzygy_checks(Cs)
    print('A rank-two full-support example passes the cubic syzygy tests (not a helix test).')


def rotation_example():
    Cs = [sp.ones(3), sp.Matrix([[1,1,1],[1,2,3],[2,3,4]]),
          sp.Matrix([[1,2,3],[2,3,1],[3,5,4]])]
    syzygy_checks(Cs)
    u,v,w = coords = sp.symbols('u v w')
    for step in range(3):
        L = sp.Matrix.hstack(*(permutation_products(C) for C in Cs)).T
        assert L.rank() == 3
        K = sp.Matrix.hstack(*L.nullspace())
        for e,o in product(range(3),range(3,6)):
            assert K[[e,o],:].rank() == 2
            assert L[:,[e,o]].rank() == 2
        bad = [i for i,C in enumerate(Cs)
               if all(D.multiply_elementwise(C.applyfunc(lambda t:1/t)).det()==0 for D in Cs)]
        assert bad == ([0] if step == 0 else [])
        z = K*sp.Matrix(coords)
        f = sp.prod(z[:3])-sp.prod(z[3:])
        assert f != 0
        if step == 0:
            ds = [sp.diff(f,t) for t in coords]
            for t in coords:
                others = [q for q in coords if q != t]
                gb = sp.groebner([d.subs(t,1) for d in ds],*others)
                assert list(gb) == [q-1 for q in others]
            Hess = sp.hessian(f.subs(w,1),(u,v)).subs({u:1,v:1})
            assert Hess.det() != 0
        else:
            # Factorization into a line and a nondegenerate conic.
            fac = sp.factor_list(f, *coords)[1]
            assert sorted((sp.Poly(g,*coords).total_degree(),e) for g,e in fac) == [(1,1),(2,1)]
            conic = next(g for g,e in fac if sp.Poly(g,*coords).total_degree()==2)
            assert sp.hessian(conic,coords).det() != 0
        Cs = [sp.Matrix(3,3,lambda i,j: Cs[j][k,i]) for k in range(3)]
    print('Rotation example: integral bad kernel, two reducible good grids; all cubic syzygy tests pass.')


if __name__ == '__main__':
    matching_checks()
    euler_checks()
    classical_checks()
    dual_checks()
    coefficient_examples()
    rotation_example()
