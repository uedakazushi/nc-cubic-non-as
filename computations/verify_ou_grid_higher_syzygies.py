#!/usr/bin/env python3
"""Exact fourth-path formula checks and symbolic exclusion of triple grid failure.

The fraction-field ranks give upper bounds at EVERY specialization because
all larger minors of the original polynomial matrices vanish identically.
No pivot is required to survive specialization.
"""
from collections import Counter
from itertools import product
from sympy import QQ, symbols, eye, ones, zeros, BlockMatrix

I = range(3)
CUBE = set(product((1, 2), repeat=3))
TRIPLES = list(product(I, repeat=3))


def rotate(c, n):
    return {t: c[t[-n:] + t[:-n]] for t in TRIPLES}


def path_matrix(c, length, start, end, K):
    words = list(product(I, repeat=length-1))
    index = {w: i for i, w in enumerate(words)}
    rows = []
    for pos in range(length-1):
        r = rotate(c, pos % 3)
        for w in words:
            if w[pos] != 0:
                continue
            row = {}
            for a in I:
                v = w[:pos] + (a,) + w[pos+1:]
                full = (start,) + v + (end,)
                coef = r[full[pos:pos+3]]
                if coef:
                    row[index[v]] = coef
            rows.append(row)
    return len(words), rows


def exact_rank(rows, K):
    # Sparse Gaussian elimination over an EXACT field. Prefer unit pivots
    # and then a small fill-in score; all operations preserve rank.
    rank = 0
    while True:
        rows = [r for r in rows if r]
        if not rows:
            return rank
        count = Counter(j for r in rows for j in r)
        _, _, i, j = min(
            (0 if v in (K.one, -K.one) else 1,
             (len(r)-1)*(count[j]-1), i, j)
            for i, r in enumerate(rows) for j, v in r.items()
        )
        r = rows.pop(i)
        v = r.pop(j)
        if v != K.one:
            r = {k: w/v for k, w in r.items()}
        for m in rows:
            x = m.pop(j, K.zero)
            if not x:
                continue
            for k, w in r.items():
                val = m.get(k, K.zero)-x*w
                if val:
                    m[k] = val
                else:
                    m.pop(k, None)
        rank += 1


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


if __name__ == '__main__':
    main()
