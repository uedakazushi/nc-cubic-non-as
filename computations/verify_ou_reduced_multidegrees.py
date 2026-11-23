#!/usr/bin/env python3
"""Enumerate necessary cohomological multidegree patterns on a reduced cubic.
Components are labelled. Within each block the three vectors are sorted;
a common tensor twist makes the least vector of the first block zero.
This is NOT a classification of OU categories or an existence test for I_3.
"""
from itertools import product, combinations_with_replacement, permutations
from collections import Counter


def add(*xs):
    return tuple(map(sum, zip(*xs)))


def mul(n, x):
    return tuple(n*a for a in x)


def sub(x, y):
    return add(x, mul(-1, y))


def forward(x, y):
    """For total degree d>=0: every proper subcurve has deg(y-x)>=-1.
    For 2 or 3 components this is equivalent to the coordinate bounds.
    """
    d = sum(y)-sum(x)
    assert d >= 0
    return all(-1 <= t <= d+1 for t in sub(y, x))


def blocks(pool):
    return [T for T in combinations_with_replacement(sorted(pool), 3)
            if all(forward(x, y) for x in T for y in T)]


def normalize(U, V, W):
    a = min(U)
    return tuple(tuple(sorted(sub(x, a) for x in T)) for T in (U, V, W))


def invariants(U, V, W):
    A, B, C = add(*U), add(*V), add(*W)
    return add(A, C, mul(-2, B)), add(mul(5, B), mul(-3, A), mul(-2, C))


def classical(U, V, W):
    """One block ordering has plane polarization H and six smooth point degrees."""
    c = len(U[0])
    zero = (0,)*c
    delta, w = invariants(U, V, W)
    if delta != zero or not all(t > 0 for t in w):
        return False
    units = {tuple(int(i == j) for i in range(c)) for j in range(c)}
    Hs = {(1, 2), (2, 1)} if c == 2 else {(1, 1, 1)}
    B = add(*V)
    for a, b, v in permutations(U):
        H = add(B, mul(-2, a), mul(-1, v))
        if H not in Hs:
            continue
        P = [sub(x, a) for x in V] + [add(B, mul(-1, a), mul(-1, v), mul(-1, x)) for x in W]
        if any(p not in units for p in P):
            continue
        assert sub(b, a) == sub(add(*P), mul(2, H))
        assert sub(v, a) == sub(add(*P[:3]), H)
        assert w == sub(mul(3, H), add(*P))
        return True
    return False


def rotation(U, V, W, w):
    return normalize(V, W, tuple(add(x, w) for x in U))


def enumeration(c):
    z = (0,)*c
    p0 = [v for v in product(range(-1, 2), repeat=c) if sum(v) == 0 and v >= z]
    p1 = [v for v in product(range(-1, 3), repeat=c) if sum(v) == 1]
    p2 = [v for v in product(range(-1, 4), repeat=c) if sum(v) == 2]
    stat = Counter()
    for U in blocks(p0):
        if U[0] != z:
            continue
        for V in blocks(p1):
            if not all(forward(x, y) for x in U for y in V):
                continue
            for W in blocks(p2):
                if not all(forward(x, y) for x in U+V for y in W):
                    continue
                stat['necessary'] += 1
                delta, w = invariants(U, V, W)
                assert sum(delta) == 0 and sum(w) == 3
                if not all(t > 0 for t in w):
                    continue
                stat['positive_w'] += 1
                if delta != z:
                    continue
                stat['positive_w_delta_zero'] += 1
                T = (U, V, W)
                if classical(*T):
                    stat['classical_direct'] += 1
                found = False
                for _ in range(3):
                    assert invariants(*T) == (z, w)
                    found |= classical(*T)
                    T = rotation(*T, w)
                assert T == (U, V, W)
                assert found, (c, U, V, W)
                stat['classical_after_rotation'] += 1
    return dict(stat)


def main():
    expected = {2: (180, 14, 6, 6, 6), 3: (15012, 27, 9, 7, 9)}
    keys = ('necessary', 'positive_w', 'positive_w_delta_zero',
            'classical_direct', 'classical_after_rotation')
    for c in (2, 3):
        stat = enumeration(c)
        assert tuple(stat[k] for k in keys) == expected[c]
        print(c, stat)


if __name__ == '__main__':
    main()
