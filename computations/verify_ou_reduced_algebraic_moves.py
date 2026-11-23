#!/usr/bin/env python3
"""Sufficient algebraic line-bundle moves on reduced plane cubics.

No mapping-class model, gluing parameter, or existence of an intrinsic
marking is assumed by the integer computation. Every directed edge uses
only the numerical inequalities in the elementary-move lemma of the paper.
"""
from collections import Counter
from itertools import permutations, product
from verify_ou_reduced_multidegrees import (
    add, mul, sub, normalize, invariants, rotation, patterns,
)


def bounded(x, a, b):
    return all(a <= t <= b for t in x)


def unit(x):
    return sum(x) == 1 and bounded(x, 0, 1)


def reference(T):
    """The classical reference criterion, allowing nonample anticanonical degree."""
    U, V, W = T
    delta, w = invariants(*T)
    assert not any(delta)
    B = add(*V)
    for a, b, v in permutations(U):
        H = add(B, mul(-2, a), mul(-1, v))
        if not all(t > 0 for t in H):
            continue
        P = [sub(x, a) for x in V] + [add(B, mul(-1, a), mul(-1, v), mul(-1, x)) for x in W]
        if all(unit(p) for p in P):
            assert sum(H) == 3
            assert sub(b, a) == sub(add(*P), mul(2, H))
            assert sub(v, a) == sub(add(*P[:3]), H)
            assert w == sub(mul(3, H), add(*P))
            return True
    return False


def terminal(T):
    c = len(T[0][0])
    delta, w = invariants(*T)
    assert not any(delta)
    S = T
    for _ in range(3):
        assert invariants(*S) == (delta, w)
        if reference(S):
            return 'positive' if min(w) > 0 else 'boundary'
        S = rotation(*S, w)
    assert S == T
    zeros = sum(all(d[j] == 0 for B in T for d in B) for j in range(c))
    if zeros:
        assert zeros < c
        return f'contract_{zeros}'
    return None


def edges(T):
    """Each edge includes the explicit (v,q,k) witnessing the move."""
    U, V, W = T
    c = len(U[0])
    for v in product(range(-2, 4), repeat=c):
        if sum(v) != 1:
            continue
        if not (all(bounded(sub(v, x), -1, 2) for x in U)
                and all(bounded(sub(x, v), -1, 1) for x in V)
                and all(bounded(sub(x, v), -1, 2) for x in W)):
            continue
        for i, j in permutations(range(c), 2):
            q = tuple(int(a == i) - int(a == j) for a in range(c))
            if not all(bounded(add(sub(x, v), q), -1, 1) for x in V):
                continue
            for k in (1, -1):
                if k == 1:
                    if not all(unit(sub(v, x)) for x in U):
                        continue
                    if not all(bounded(add(sub(x, v), mul(2, q)), -1, 2) for x in W):
                        continue
                else:
                    if not all(unit(sub(x, v)) for x in W):
                        continue
                    if not all(bounded(add(sub(x, v), mul(2, q)), -2, 1) for x in U):
                        continue
                target = normalize(*(tuple(add(x, mul(b if k == 1 else 2-b, q))
                                           for x in B) for b, B in enumerate(T)))
                yield target, (v, q, k)


def classify(c):
    ts = [T for T in patterns(c) if not any(invariants(*T)[0])]
    universe = set(ts)
    labels = {T: terminal(T) for T in ts}
    paths = {T: () for T in ts if labels[T] is not None}
    adj = {T: list(edges(T)) for T in ts}
    assert all(S in universe for es in adj.values() for S, _ in es)
    changed = True
    while changed:
        changed = False
        for T in ts:
            if T in paths:
                continue
            for S, witness in adj[T]:
                if S in paths:
                    paths[T] = ((S, witness),) + paths[S]
                    labels[T] = labels[S]
                    changed = True
                    break
    counts = Counter()
    unresolved_w = Counter()
    for T in ts:
        if T not in paths:
            counts['remaining'] += 1
            unresolved_w[tuple(sorted(invariants(*T)[1]))] += 1
            continue
        path = paths[T]
        S = T
        for target, witness in path:
            assert (target, witness) in adj[S]
            S = target
        assert labels[T] == terminal(S)
        counts[('changed_' if path else 'direct_') + labels[T]] += 1
    expected = {
        2: {'direct_positive': 6, 'direct_boundary': 4, 'direct_contract_1': 2,
            'changed_positive': 4},
        3: {'direct_positive': 9, 'direct_boundary': 36, 'direct_contract_1': 42,
            'direct_contract_2': 3, 'changed_boundary': 36,
            'changed_contract_1': 18, 'remaining': 24},
    }
    assert dict(counts) == expected[c], counts
    print(c, dict(counts), 'edges', sum(map(len, adj.values())),
          'max_steps', max(map(len, paths.values())), 'unresolved_w', dict(unresolved_w))


if __name__ == '__main__':
    for c in (2, 3):
        classify(c)
