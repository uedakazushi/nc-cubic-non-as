#!/usr/bin/env python3
"""Universal identities for the residual IV relation family.

The proof of the heart/rank obstruction is categorical; this script checks
only its finite path-space, formal-disc relation and Euler-matrix identities.
"""
from itertools import product
import sympy as sp


def relation(i, k, x):
    if i == k:
        return sp.Matrix([0 if j == i else 1 if j == (i+1)%3 else -1
                          for j in range(3)])
    return sp.Matrix([1 if j == i else -1 if j == k else x[i, k]
                      for j in range(3)])


def path_space_checks():
    x = {(i, k): sp.Symbol(f'x{i}{k}') for i, k in product(range(3), repeat=2) if i != k}
    eye = sp.eye(3)
    # In k^3/<relation>, the two paths via j != i span dimension
    # two off the diagonal (the omitted path has nonzero coefficient),
    # and dimension one on the diagonal. The same checks give the dual
    # injective complex. The diagonal kernel has coefficients +/-1.
    for i, k in product(range(3), repeat=2):
        r = relation(i, k, x)
        for omitted in (i, k):
            basis = [eye[:, j] for j in range(3) if j != omitted]
            det = sp.Matrix.hstack(*basis, r).det()
            if i == k:
                assert det == 0 and r[omitted] == 0
                assert all(abs(r[j]) == 1 for j in range(3) if j != omitted)
            else:
                assert abs(det) == 1
    print('Universal projective/injective path-space identities verified.')


def disc_checks():
    t = sp.Symbol('t')
    a = sp.symbols('a0:3')
    b = sp.symbols('b0:3')
    x = {(i, k): -a[i]-b[k] for i, k in product(range(3), repeat=2) if i != k}
    assert sp.expand(x[0, 1]+x[1, 2]+x[2, 0]-x[0, 2]-x[2, 1]-x[1, 0]) == 0
    # Right modules: arrows act backwards, by multiplication by these scalars.
    def sign(i, j):
        return -1 if {i, j} == {1, 2} else 1
    A = sp.Matrix(3, 3, lambda i, j: 1-a[i]*t if i == j else sign(i, j)*t)
    B = sp.Matrix(3, 3, lambda j, k: 1+b[j]*t if j == k else sign(j, k)*t)
    for i, k in product(range(3), repeat=2):
        r = relation(i, k, x)
        assert sp.expand(sum(r[j]*A[i, j]*B[j, k] for j in range(3))) == 0
        if i == k:
            assert A[i, i].subs(t, 0) == B[i, i].subs(t, 0) == 1
        else:
            assert abs(A[i, k]/t) == abs(B[i, k]/t) == 1
    print('All nine formal-disc relation identities verified over Q[a,b,t].')


def euler_checks():
    I, J, Z = sp.eye(3), sp.ones(3), sp.zeros(3)
    D = I.row_join(J).row_join(2*J).col_join(Z.row_join(I).row_join(J)).col_join(Z.row_join(Z).row_join(I))
    S = D.inv()*D.T
    N = S-sp.eye(9)
    ones = sp.ones(9, 1)
    assert N**3 == sp.zeros(9)
    assert ones.T*S == ones.T
    assert (N**2).T*D == 3*sp.ones(9)
    M = N.T*D
    assert M+M.T == 6*sp.ones(9)
    basis = sp.eye(9)
    qs = [basis[:, i]+basis[:, 6+i]-sum((basis[:, 3+j] for j in range(3) if j != i), sp.zeros(9, 1)) for i in range(3)]
    assert sum(qs, sp.zeros(9, 1)) == D.inv()*ones
    for i, q in enumerate(qs):
        assert S*q == q
        assert D*q == basis[:, i]+basis[:, 3+i]+basis[:, 6+i]
    assert sp.Matrix(3, 3, lambda i, j: (qs[i].T*D*qs[j])[0]) == 3*sp.eye(3)-sp.ones(3)
    # Fixed Euler-dimension vectors have equal sums in all three blocks.
    for q in N.nullspace():
        d = D*q
        assert sum(d[:3]) == sum(d[3:6]) == sum(d[6:])
    print('Euler rank, affine A2 pairing and fixed block-sum identities verified.')


if __name__ == '__main__':
    path_space_checks()
    disc_checks()
    euler_checks()
