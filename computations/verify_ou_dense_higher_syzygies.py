#!/usr/bin/env python3
"""Check cube supports and two rational higher-path examples (not a proof for all parameters)."""
from itertools import product, permutations
from sympy import Matrix, Rational

I = range(3)
P = list(product((1, 2), repeat=3))
supports = []
for bits in product((0, 1), repeat=8):
    S = {p for p, b in zip(P, bits) if b}
    if any({p[d] for p in S} != {1, 2} for d in I):
        continue
    faces = [{tuple(p[e] for e in I if e != d) for p in S if p[d] == v}
             for d in I for v in (1, 2)]
    if all(F == set(product({q[0] for q in F}, {q[1] for q in F})) for F in faces):
        supports.append(S)
assert sorted(map(len, supports)) == [2, 2, 2, 2, 8]
assert all(all(a != b for a, b in zip(*S)) for S in supports if len(S) == 2)


def rotate(c, n):
    return {t: c[t[-n:] + t[:-n]] for t in product(I, repeat=3)}


def slice_matrix(c, k):
    return Matrix(3, 3, lambda i, j: c[i, j, k])


def quotient_dim(c, length, start=0, end=0):
    words = [(start,) + w + (end,) for w in product(I, repeat=length-1)]
    idx = {w: i for i, w in enumerate(words)}
    rows = []
    for j in range(1, length):
        r = rotate(c, (j-1) % 3)
        for w in words:
            if w[j] != 0:
                continue
            row = [0] * len(words)
            for a in I:
                ww = w[:j] + (a,) + w[j+1:]
                row[idx[ww]] = r[ww[j-1:j+2]]
            rows.append(row)
    return len(words) - Matrix(rows).rank()


def check_example(c):
    assert all(c.values())
    for n in I:
        r = rotate(c, n)
        C = [slice_matrix(r, k) for k in I]
        for k, l in permutations(I, 2):
            ratio = Matrix(3, 3, lambda i, j: Rational(C[l][i, j], C[k][i, j]))
            assert ratio.rank() > 1
            if k == 0:
                assert ratio.det() == 0
        assert quotient_dim(r, 3) == 4
        assert quotient_dim(r, 3, 0, 1) == 3


x, y, z = (0, 1, 2), (0, 3, 5), (0, 7, 11)
dense = {(i, j, k): 1+x[i]*y[j]*z[k] for i, j, k in product(I, repeat=3)}
sparse = {t: 1 for t in product(I, repeat=3)}
sparse[1, 1, 1], sparse[2, 2, 2] = 2, 3
for c in (dense, sparse):
    check_example(c)
assert quotient_dim(dense, 4) == 5
assert quotient_dim(dense, 5) == 8
assert quotient_dim(sparse, 4) == 6
print('PASS: cube supports 4 antipodal pairs + full cube; path dimensions dense (5, 8), sparse (6)')
