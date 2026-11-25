#!/usr/bin/env python3
"""All points of the square-difference divisor in the first one-zero family."""
from sympy import cancel, symbols

from verify_ou_one_zero_closed_divisors import LABEL
from verify_ou_one_zero_first_orbits import verify_orbits
from verify_ou_one_zero_first_two_slices import graph_count
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_two_zero_closed_strata import local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import surface_chart


EXTRA = {('axis', 0, 0): (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (7, (5, 2)))}
EXPECTED = [(0, 0, 2, 29), (3, 0, 1, 18), (0, 0, 0, 0), (5, 0, 1, 18)]


def verify():
    coeff, ps, units = family(0)
    a, c, d, e = ps
    b, t, w, y = symbols('b t w y')
    central = a*c*d-c-d+1
    guards = units | {d-1, central}
    back = {b: (a*d*e-e+1)/a, t: 1/(a*d-1)}
    change = {d: (t+1)/(t*a), e: t*(a*b-1)}
    for p, val in back.items():
        assert factors(val) <= guards and cancel(val.subs(change, simultaneous=True)-p) == 0
    for p, val in change.items():
        assert cancel(val.subs(back, simultaneous=True)-p) == 0
    polynomial = b*t*(a-1)**2+b*(a*c-1)**2-a*(c-1)**2
    f, = factors(polynomial.subs(back, simultaneous=True))-guards
    old = a*a*c*d*e-a*c*e-a*d*e-a*e+2*a+2*e-2
    junction = a*a*c*d-a*c-a*d*e+a+e-1
    linear = a*a*c*d-a*c-a*d-a*e+2*a+e-1
    quadratic, = factors((b*t*(a-1)**2-a*c*(b-1)+b-a).subs(back, simultaneous=True))-units
    h = (back[b]-1)*back[t]
    ratio, = factors(a*back[b]*(c*e-h)**2-(e-1)**2*(h-e)**2)-guards
    orbits = verify_orbits(coeff, ps, units, {d-1, central, old})
    indices = [next(i for i, orbit in enumerate(orbits) if g in orbit)
               for g in (junction, linear, quadratic, ratio, f)]
    assert len(set(indices)) == 5

    # The two other orbits admit common coordinates (p, s, v), with e a graph.
    h0 = (back[b]-1)*back[t]
    residual = (back[b]*(c*e-h0)**2-(e-h0)*(h0*h0-2*h0*c+c*c*e),
                back[b]*e*(c*e-h0)-(e-h0)*(h0*e-h0+c*e-e))
    p, s0, v0 = symbols('p s0 v0')
    new_b = (1-v0)*(1+s0)
    inverse = {p: cancel((c*e-h0)/h0), v0: cancel(h0/e)}
    inverse[s0] = cancel((back[b]+inverse[v0]-1)/(1-inverse[v0]))
    assert all(factors(z) <= guards for z in inverse.values())
    for equation, new_e in zip(residual, (1+s0*p*p, 1/v0+s0*p)):
        remaining, = factors(equation)-guards
        indices.append(next(i for i, orbit in enumerate(orbits) if remaining in orbit))
        # Keep e free to check its unique value and every coordinate denominator.
        general = {a: s0/(v0*(1+s0)), c: v0*(1+p),
                   d: (1+s0)*(new_b+v0*e-1)/(e*s0)}
        for z, val in general.items():
            assert cancel(val.subs(inverse, simultaneous=True)-z) == 0
            assert factors(cancel(val.as_numer_denom()[1].subs(inverse, simultaneous=True))) <= guards
        for z, val in inverse.items():
            assert cancel(val.subs(general, simultaneous=True)-z) == 0
        unit = cancel(equation.subs(general, simultaneous=True)/(e-new_e))
        assert unit and factors(cancel(unit.subs(inverse, simultaneous=True))) <= guards
    assert len(set(indices)) == len(orbits) == 7

    base = apply_sub(coeff, change)
    local = local_units(guards, change) | {b, t}
    assert factors(b*(a-1)**2) <= local
    t_value = cancel((a*(c-1)**2-b*(a*c-1)**2)/(b*(a-1)**2))
    assert cancel(polynomial.subs(t, t_value)) == 0
    base = apply_sub(base, {t: t_value})
    local = local_units(local, {t: t_value})
    to_w = {c: (w+1)/(a+w), b: a/y}
    from_w = {w: (a*c-1)/(1-c), y: a/b}
    assert factors(1-c) | factors(b) <= local
    for p, val in to_w.items():
        assert cancel(val.subs(from_w, simultaneous=True)-p) == 0
        assert factors(cancel(val.as_numer_denom()[1].subs(from_w, simultaneous=True))) <= local
    for p, val in from_w.items():
        assert cancel(val.subs(to_w, simultaneous=True)-p) == 0
    assert cancel(t_value.subs(to_w, simultaneous=True)-(y-w*w)/(a+w)**2) == 0
    base = apply_sub(base, to_w)
    local = local_units(local, to_w) | {y, a+w}

    v, u, s = params = symbols('x0:3')
    axis = {a: v, w: 0, y: s}
    result = graph_count('axis', 0, apply_sub(base, axis), (v, s), local_units(local, axis), extra=EXTRA)
    assert result == (2, 12, 0, 0)
    open_chart = {a: (v-1)*u, w: u, y: (1+s*v)*u*u}
    inverse = {v: (a+w)/w, u: w, s: (y-w*w)/(w*(a+w))}
    for p, val in inverse.items():
        assert factors(val.as_numer_denom()[1]) <= local | {w}
        assert cancel(val.subs(open_chart, simultaneous=True)-p) == 0
    for p, val in open_chart.items():
        assert cancel(val.subs(inverse, simultaneous=True)-p) == 0
    coeff = apply_sub(base, open_chart)
    guards = local_units(local | {w}, open_chart)
    cuts = parametric_check(permute(coeff, LABEL['C']), params)-guards-{-g for g in guards}
    assert len(cuts) == 4
    results = []
    for i, g in enumerate(sorted(cuts, key=str)):
        sub, vs = surface_chart(g, params, guards)
        results.append(graph_count('difference', i, apply_sub(coeff, sub), vs,
                                   local_units(guards, sub), extra=EXTRA))
    assert results == EXPECTED, results
    print('PASS complete first-family square-difference divisor', result, results, flush=True)


if __name__ == '__main__':
    verify()
