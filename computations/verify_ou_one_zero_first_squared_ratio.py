#!/usr/bin/env python3
"""The first-family divisor ab(ce-h)^2=(e-1)^2(h-e)^2, on all six graphs."""
from sympy import Poly, cancel, symbols

from verify_ou_one_zero_closed_divisors import LABEL
from verify_ou_one_zero_first_orbits import verify_orbits
from verify_ou_one_zero_first_two_slices import graph_count
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_two_zero_closed_strata import local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import surface_chart


ROTATE = ((1, 2, 0), (2, 0, 1), (1, 2, 0))
EXTRA = {
    ('ratio', 0, 0): (ROTATE, (7, (4, 3))),
    ('ratio', 1, 0): (ROTATE, (7, (2, 3))),
}
EXPECTED = [(0, 1, 2, 23), (6, 1, 1, 24), (5, 0, 3, 50),
            (2, 0, 2, 28), (3, 0, 2, 13), (0, 0, 3, 34)]


def verify():
    co, ps, units = family(0)
    a, c, d, e = ps
    b = (a*d*e-e+1)/a
    h = (a*d*e-a-e+1)/(a*(a*d-1))
    central = a*c*d-c-d+1
    old = a*a*c*d*e-a*c*e-a*d*e-a*e+2*a+2*e-2
    junction = a*a*c*d-a*c-a*d*e+a+e-1
    linear = a*a*c*d-a*c-a*d-a*e+2*a+e-1
    quadratic, = factors(b*(a-1)**2/(a*d-1)-a*c*(b-1)+b-a)-units
    guards = units | {d-1, central}
    f, = factors(a*b*(c*e-h)**2-(e-1)**2*(h-e)**2)-guards
    orbits = verify_orbits(co, ps, units, {d-1, central, old})
    indices = [next(i for i, orbit in enumerate(orbits) if g in orbit)
               for g in (junction, linear, quadratic, f)]
    assert len(set(indices)) == 4

    q, y, v = params = symbols('x0:3')
    back = {q: cancel((e-1)*(h-e)/(c*e-h)), y: e, v: cancel(h/e)}
    assert all(factors(val) <= guards for val in back.values())
    assert factors(back[q]**2-a*b)-guards == {f}
    B = 1+(q*q-1)*v
    go = {a: q*q/B, c: v+(v-1)*(y-1)/q, d: (q*q+y-1)*B/(q*q*y), e: y}
    # On the divisor B=b; hence each forward denominator is a unit.
    assert cancel(1+(a*b-1)*back[v]-b) == 0 and factors(b) <= guards
    for p, val in go.items():
        numerator = cancel(val.subs(back, simultaneous=True)-p).as_numer_denom()[0]
        assert Poly(numerator, *ps).rem(Poly(f, *ps)).is_zero
    for p, val in back.items():
        assert cancel(val.subs(go, simultaneous=True)-p) == 0
    assert cancel(f.subs(go, simultaneous=True)) == 0
    coeff = apply_sub(co, go)
    local = local_units(guards, go) | {q}
    cuts = parametric_check(permute(coeff, LABEL['C']), params)-local-{-g for g in local}
    assert len(cuts) == 6
    results = []
    for i, g in enumerate(sorted(cuts, key=str)):
        sub, vs = surface_chart(g, params, local)
        results.append(graph_count('ratio', i, apply_sub(coeff, sub), vs,
                                   local_units(local, sub), extra=EXTRA))
    assert results == EXPECTED, results
    print('PASS complete first-family squared-ratio divisor', results, flush=True)


if __name__ == '__main__':
    verify()
