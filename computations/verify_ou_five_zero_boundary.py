#!/usr/bin/env python3
"""Characteristic-zero rank witness for the 117374702 boundary u=v=1."""
from sympy import QQ, symbols
from verify_ou_zero_supports import normalize, path_dim


def main():
    c, K, free = normalize(117374702)
    assert free == [(2, 1, 1), (2, 1, 2), (2, 2, 2)]
    u, v, w = symbols('u:3')
    F = QQ.frac_field(w)
    c = {q: F.from_sympy(K.to_sympy(z).subs({u: 1, v: 1}))
         for q, z in c.items()}
    # All coefficients are polynomials in w. The fraction-field rank is an
    # upper bound on the specialized relation rank, even at exceptional w.
    assert path_dim(c, 7, 0, 0, 2, F) == 13
    print('PASS 117374702 u=v=1: rank 716, length-7 quotient 13 > 12')


if __name__ == '__main__':
    main()
