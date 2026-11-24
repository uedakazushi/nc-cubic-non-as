#!/usr/bin/env python3
"""Characteristic-zero path rank for c_021=c_022=0.

python-flint 0.9.0 accelerates the last polynomial elimination. The existing
SymPy backend computes the same fraction-field rank if FLINT is unavailable.
"""
from collections import Counter
from itertools import combinations
from math import lcm

from verify_ou_four_zero_supports import laurent_rank, rows_for
from verify_ou_zero_supports import TRIPLES, normalize, support


MASK = 134217343
SMALL_EXCLUSIONS = (134217343, 134217311, 134217335, 134200959, 134213247)


def flint_rank(rows, field):
    from flint import fmpz_mpoly_ctx

    ring = fmpz_mpoly_ctx.get(tuple(map(str, field.symbols)))
    zero = ring.constant(0)

    def to_integer_row(row):
        denominator = lcm(*(int(b.denominator) for f in row.values() for b in f.values()))
        return {j: ring.from_dict({e: int(b * denominator) for e, b in f.items()})
                for j, f in row.items()}

    def primitive(row):
        entries = sorted(row.values(), key=len)
        divisor = entries[0]
        for f in entries[1:]:
            if divisor.is_one():
                break
            divisor = divisor.gcd(f)
        if divisor.leading_coefficient() * row[min(row)].leading_coefficient() < 0:
            divisor = -divisor
        if divisor.is_one():
            return row
        return {j: f / divisor for j, f in row.items()}

    def distinct(rows):
        # Canonical polynomial strings allow hashing multivariate FLINT entries.
        return list({tuple(sorted((j, str(f)) for j, f in row.items())): row
                     for row in rows}.values())

    rows = distinct(primitive(to_integer_row(row)) for row in rows if row)
    rank = 0
    while rows:
        counts = Counter(j for row in rows for j in row)

        def cost(entry):
            i, j = entry
            row = rows[i]
            return (len(row[j]) != 1, (len(row) - 1) * (counts[j] - 1),
                    len(row[j]), sum(len(f) for f in row.values()))

        i, j = min(((i, j) for i, row in enumerate(rows) for j in row), key=cost)
        pivot = rows.pop(i)
        v = pivot[j]
        reduced = []
        for row in rows:
            if j not in row:
                reduced.append(row)
                continue
            x = row[j]
            divisor = v.gcd(x)
            a, b = v / divisor, x / divisor
            result = {k: a * row.get(k, zero) - b * pivot.get(k, zero)
                      for k in row.keys() | pivot.keys() if k != j}
            result = {k: value for k, value in result.items() if value}
            if result:
                reduced.append(primitive(result))
        rows = distinct(reduced)
        rank += 1
    return rank


def main():
    for mask in SMALL_EXCLUSIONS:
        zeros = set(TRIPLES) - set(support(mask))
        assert any(sum(x == y for x, y in zip(a, b)) == 2 for a, b in combinations(zeros, 2))
    coeffs, field, free = normalize(MASK)
    assert set(TRIPLES) - set(support(MASK)) == {(0, 2, 1), (0, 2, 2)}
    assert free == [(1, 1, 1), (1, 1, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2)]
    size, rows = rows_for(coeffs, 8, 1, 2, 0, field)
    assert size == 2187 and len(rows) == 5103
    try:
        import flint
    except ImportError:
        finish = None
    else:
        finish = flint_rank
    rank = laurent_rank(rows, field, finish=finish)
    assert rank == 2170 and size - rank == 17
    print('PASS length 8: 5103 rows, 2187 columns, generic rank 2170, quotient 17 > 15')


if __name__ == '__main__':
    main()
