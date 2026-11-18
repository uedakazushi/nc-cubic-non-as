#!/usr/bin/env python3
"""Exact smoothness check for the degree-three Clifford discriminant.

Run with Python 3 and SymPy. All Groebner computations take place over QQ.
The four affine charts cover the Hirzebruch surface F_1.
"""

import sympy as sp


def main():
    if not __debug__:
        raise SystemExit("Verification requires assertions; do not use optimized Python.")
    s, t, c, d = sp.symbols("s t c d")
    q = sp.Matrix(
        [[c * s**2 + d * t, c * t**2, d],
         [c * t**2, c * t**2 + d * s, 0],
         [d, 0, c]]
    )
    delta = sp.expand(q.det())
    expected = (c**3 * s**2 * t**2 - c**3 * t**4
                + c**2 * d * s**3 + c**2 * d * t**3
                + c * d**2 * s * t - c * d**2 * t**2 - d**3 * s)
    assert sp.expand(delta - expected) == 0
    print("Coefficient field: QQ")
    print("Discriminant:", delta)
    for base in (s, t):
        for fibre in (c, d):
            affine = sp.expand(delta.subs({base: 1, fibre: 1}))
            variables = tuple(v for v in (s, t, c, d) if v not in (base, fibre))
            generators = [affine] + [sp.diff(affine, v) for v in variables]
            basis = sp.groebner(generators, *variables, order="lex", domain=sp.QQ)
            assert list(basis) == [1], (base, fibre, list(basis))
            print(f"Chart {base}={fibre}=1; variables {variables}; singularity Groebner basis: [1]")
    print("PASS: D is smooth in all four covering affine charts.")


if __name__ == "__main__":
    main()
