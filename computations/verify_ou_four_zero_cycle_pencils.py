#!/usr/bin/env python3
"""Determine every divisorial failure of the four-zero complementary conic charts."""
from collections import Counter

from verify_ou_complement_cycles import verify_factor
from verify_ou_complement_line_conic import boundary_factors
from verify_ou_four_zero_boundary import tensor, u, v, w, q
from verify_ou_zero_supports import I, TRIPLES


CASES = (
    (134135535, (0,2,2), u, {v-1,w-1,q-1}),
    (134135535, (1,2,1), v, {u-1,w-1,q-1}),
    (134135535, (1,1,2), w, {u-1,v-1,q-1}),
    (134135535, (0,1,1), q, {u-1,v-1,w-1}),
    (134147823, (1,2,1), v, {u+w-1,q-1}),
    (134147823, (0,1,1), q, {v-1,u+v*w-v}),
    (134151407, (1,2,1), v, {u+w-1,w+q-1}),
    (134151918, (0,2,2), u, {q-1,w+1}),
    (134151918, (0,1,1), q, {u-1,w+1}),
)


def main():
    counts = Counter()
    for mask, zero, fixed, expected in CASES:
        coefficients = tensor(mask, {fixed: 1})
        labels = [[j] + [i for i in I if i != j] for j in zero]
        c = {t: coefficients[tuple(labels[i][t[i]] for i in I)] for t in TRIPLES}
        variables = [i for i in (u,v,w,q) if i != fixed]
        assert c[0,0,0] == 0 and all(c[t] == 1 for t in TRIPLES if all(t))
        factors = boundary_factors(c, variables)
        assert factors == expected
        for f in sorted(factors, key=str):
            count, chart = verify_factor(c, variables, factors, f)
            counts[count] += 1
            print("PASS", mask, zero, f, count, chart)
    assert counts == {3: 6, 4: 16}
    print("PASS all 22 four-zero divisor incidences", dict(counts))


if __name__ == "__main__":
    main()
