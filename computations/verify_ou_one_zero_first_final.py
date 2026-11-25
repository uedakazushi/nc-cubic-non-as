#!/usr/bin/env python3
"""The last two orbit divisors of the AOU one-zero first family."""
from sympy import Poly, cancel, symbols

from verify_ou_one_zero_closed_divisors import LABEL
from verify_ou_one_zero_first_two_slices import FLIP, graph_count
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_two_zero_closed_strata import local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import surface_chart


EXTRA = {
    ('product', 0, 1): (FLIP, (5, (4, 4))),
    ('product', 5, 0): (FLIP, (7, (3, 2))),
    ('product', 5, 1): (FLIP, (7, (3, 3))),
}
EXPECTED = {
    'square': [(0, 0, 3, 34), (6, 0, 4, 52), (4, 0, 2, 22),
               (3, 0, 1, 4), (0, 0, 2, 20), (0, 0, 3, 36)],
    'product': [(3, 7, 1, 13), (5, 0, 1, 20), (6, 0, 4, 52), (3, 0, 1, 18),
                (5, 0, 3, 50), (5, 23, 0, 0), (2, 0, 2, 14), (8, 0, 2, 31)],
}


def coordinates(name):
    coeff, ps, units = family(0)
    a, c, d, e = ps
    params = symbols('x0:3')
    if name == 'square':
        p, s, v = params
        new_e = 1+s*p*p
    else:
        v, b_ratio, new_e = params
        s = b_ratio-1
        p = (new_e-1/v)/s
    new_b = (1-v)*(1+s)
    sub = {a: cancel(s/(v*(1+s))), c: cancel(v*(1+p)),
           d: cancel((1+s)*(new_b+v*new_e-1)/(new_e*s)), e: new_e}
    guards = units | {d-1, a*c*d-c-d+1}
    b = (a*d*e-e+1)/a
    h = (b-1)/(a*d-1)
    v_old = cancel(h/e)
    s_old = cancel(b/(1-v_old)-1)
    p_old = cancel(c/v_old-1)
    old = (p_old, s_old, v_old) if name == 'square' else (v_old, 1+s_old, e)
    back = dict(zip(params, old))
    for x, value in back.items():
        assert factors(value) <= guards
        assert cancel(value.subs(sub, simultaneous=True)-x) == 0
    equation = (b*(c*e-h)**2-(e-h)*(h*h-2*h*c+c*c*e) if name == 'square'
                else b*e*(c*e-h)-(e-h)*(h*e-h+c*e-e))
    f, = factors(equation)-guards
    assert cancel(f.subs(sub, simultaneous=True)) == 0
    # The graph value new_e agrees with the original unit e on this divisor.
    e_diff = cancel(new_e.subs(back, simultaneous=True)-e).as_numer_denom()[0]
    assert Poly(e_diff, *ps).rem(Poly(f, *ps)).is_zero and e in guards
    for unit in (s, 1+s, v):
        assert factors(cancel(unit.subs(back, simultaneous=True))) <= guards
    allowed = factors(s) | factors(1+s) | factors(v) | factors(new_e)
    for x, value in sub.items():
        den = value.as_numer_denom()[1]
        assert factors(den) <= allowed
        diff = cancel(value.subs(back, simultaneous=True)-x).as_numer_denom()[0]
        assert Poly(diff, *ps).rem(Poly(f, *ps)).is_zero
    return params, apply_sub(coeff, sub), local_units(guards, sub)


def graph(name, f, params, local):
    r, b_ratio, e = params
    if name != 'product' or (f-(r*e-1)**2-r*r*(e-1)).expand() != 0:
        return surface_chart(f, params, local)
    w = symbols('w')
    back = e-1/r
    sub = {e: 1-w*w, r: 1/(1-w-w*w)}
    assert r in local and cancel(f.subs(sub, simultaneous=True)) == 0
    assert cancel(back.subs(sub, simultaneous=True)-w) == 0
    for x, value in sub.items():
        diff = cancel(value.subs(w, back)-x).as_numer_denom()[0]
        assert Poly(diff, r, e).rem(Poly(f, r, e)).is_zero
    # On the original divisor the putative denominator equals 1/r.
    diff = cancel(1-back-back**2-1/r).as_numer_denom()[0]
    assert Poly(diff, r, e).rem(Poly(f, r, e)).is_zero
    return sub, (b_ratio, w)


def verify_divisor(name):
    params, coeff, local = coordinates(name)
    cuts = parametric_check(permute(coeff, LABEL['C']), params)-local-{-f for f in local}
    assert len(cuts) == len(EXPECTED[name])
    counts = []
    for i, f in enumerate(sorted(cuts, key=str)):
        sub, vs = graph(name, f, params, local)
        graph_units = local_units(local, sub)
        assert set(vs) <= graph_units
        for value in sub.values():
            assert factors(value.as_numer_denom()[1]) <= graph_units
        counts.append(graph_count(name, i, apply_sub(coeff, sub), vs, graph_units, extra=EXTRA))
    assert counts == EXPECTED[name], counts
    print('PASS all points of first-family divisor', name, counts, flush=True)


if __name__ == '__main__':
    for divisor in EXPECTED:
        verify_divisor(divisor)
