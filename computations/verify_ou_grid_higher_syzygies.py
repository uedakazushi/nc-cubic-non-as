#!/usr/bin/env python3
"""Exact fourth-path formula checks and symbolic exclusion of triple grid failure.

The fraction-field ranks give upper bounds at EVERY specialization because
all larger minors of the original polynomial matrices vanish identically.
No pivot is required to survive specialization.
"""
from itertools import product
from sympy import QQ, symbols, eye, ones, zeros, BlockMatrix
from ou_path_rank import I, CUBE, TRIPLES, rotate, path_matrix, exact_rank


def bad_rectangle(c, i, ii, j, jj):
    ratios = [c[i,j,k]*c[ii,jj,k]/(c[i,jj,k]*c[ii,j,k]) for k in I]
    return ratios[0] == ratios[1] == ratios[2]


def check_four(c):
    for n in I:
        r = rotate(c, n)
        for i, j in product(I, repeat=2):
            b = sum(bad_rectangle(r,i,ii,j,jj)
                    for ii,jj in product(I,repeat=2) if ii != i and jj != j)
            cols, rows = path_matrix(r,4,i,j,QQ)
            assert cols-exact_rank(rows,QQ) == max(5,4+b)


def check_cube_patterns():
    # Three bad rectangles, after choosing 0 in each intersection of pairs.
    # a,b,d = 1 iff the two pairs on the corresponding axis agree.
    for a,b,d in product((1,2),repeat=3):
        Z = ({(1,1,k) for k in (1,2)} | {(a,j,1) for j in (1,2)}
             | {(i,b,d) for i in (1,2)})
        counts = []
        for axis in I:
            lines = [{p for p in CUBE if tuple(p[e] for e in I if e != axis)==t}
                     for t in product((1,2),repeat=2)]
            counts.append(sum(L <= Z for L in lines))
        assert (max(counts)>1) == (sum(x==2 for x in (a,b,d))>=2)


def main():
    J, E, Z = ones(3), eye(3), zeros(3)
    D = BlockMatrix([[E,J,2*J],[Z,E,J],[Z,Z,E]]).as_explicit()
    S = D.inv()*D.T
    assert (S.T*D)[0,3] == 5
    assert ((S**2).T*D)[0,0] == 10
    assert ((S**2).T*D)[0,3] == 12
    check_cube_patterns()
    for free in ([], [(1,2,2),(2,1,2),(2,2,1),(2,2,2)],
                 [(1,2,2),(2,1,1),(2,2,1),(2,2,2)], list(CUBE)):
        c = dict.fromkeys(TRIPLES,QQ.one)
        for t,v in zip(sorted(free),range(2,2+len(free))):
            c[t] = QQ(v)
        check_four(c)
    print('PASS: fourth-path formula and eight cube patterns', flush=True)

    variables = symbols('a b c d')
    K = QQ.frac_field(*variables)
    for free, rotation, length, expected_rank in (
        ([(1,2,2),(2,1,1),(2,2,1),(2,2,2)], 2, 6, 231),
        ([(1,2,2),(2,1,2),(2,2,1),(2,2,2)], 0, 7, 716),
    ):
        c = dict.fromkeys(TRIPLES,K.one)
        for t,v in zip(free,variables):
            c[t] = K.from_sympy(v)
        cols, rows = path_matrix(rotate(c,rotation),length,0,0,K)
        nr = len(rows)
        r = exact_rank(rows,K)
        assert r == expected_rank
        print(f'PASS: length {length}, matrix {nr} x {cols}, '
              f'generic rank {r}, quotient dimension {cols-r}', flush=True)


def check_mixed_obstructions():
    # Double torus failure with connected support: d_ijk=m_ij*z_k.
    # The analytic lower bound in the paper applies to all such tensors.
    variables = symbols('a b c d t')
    K = QQ.frac_field(*variables)
    a, b, c, d, t = map(K.from_sympy, variables)
    coeff = dict.fromkeys(TRIPLES, K.one)
    m = {(1,1):a, (1,2):b, (2,1):c, (2,2):d}
    for i,j,k in sorted(CUBE):
        coeff[i,j,k] += m[i,j]*(K.one if k == 1 else t)
    cols, rows = path_matrix(rotate(coeff,1),5,0,0,K)
    assert (len(rows), cols, exact_rank(rows,K)) == (108,81,73)
    print('PASS: double torus failure, length 5 generic quotient dimension 8', flush=True)

    # The only remaining mixed pattern after the fourth-path restrictions.
    variables = symbols('q u v w')
    K = QQ.frac_field(*variables)
    q,u,v,w = map(K.from_sympy, variables)
    coeff = dict.fromkeys(TRIPLES, K.one)
    for t, delta in {(2,2,1):q, (1,1,2):u, (1,2,2):u*w,
                     (2,1,2):u*v, (2,2,2):u*v*w}.items():
        coeff[t] += delta
    cols, rows = path_matrix(rotate(coeff,2),6,0,0,K)
    assert (len(rows), cols, exact_rank(rows,K)) == (405,243,231)
    print('PASS: mixed grid/torus failure, length 6 generic quotient dimension 12', flush=True)


if __name__ == '__main__':
    main()
    check_mixed_obstructions()
