#!/usr/bin/env python3
"""Two first-family divisors using the nontrivial stabilizer and its C chart."""

from sympy import Poly, cancel, prod, symbols

from verify_ou_one_zero_closed_divisors import LABEL, points, three_charts
from verify_ou_one_zero_first_junction import ROTATE
from verify_ou_one_zero_first_orbits import verify_orbits
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart, reindexings
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import curve_parameters
from verify_ou_two_zero_closed_strata import entire_pencil, local_cuts, local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, surface_chart


FLIP = ((1, 0, 2), (0, 1, 2), (2, 0, 1))
EXTRA = {
    ('linear', 1, 2): (FLIP, (5, (3, 4))),
    ('linear', 4, 3): (FLIP, (5, (2, 4))),
    ('linear', 5, 3): (ROTATE, (5, (4, 3))),
    ('quadratic', 1, 0): (FLIP, (5, (2, 2))),
    ('quadratic', 2, 0): (FLIP, (7, (2, 3))),
    ('quadratic', 2, 1): (FLIP, (7, (2, 2))),
    ('quadratic', 3, 1): (FLIP, (7, (5, 6))),
    ('quadratic', 3, 2): (FLIP, (5, (2, 2))),
    ('quadratic', 5, 0): (FLIP, (5, (2, 3))),
}


def restriction(fs, sub, units):
    vals = [cancel(f.subs(sub, simultaneous=True)) for f in fs]
    assert all(vals)
    assert all(factors(v.as_numer_denom()[1]) <= units | {-u for u in units} for v in vals)
    return set().union(*(factors(v) for v in vals))-units-{-u for u in units}


def normalized_parent(name):
    co, ps, units = family(0)
    a, c, d, e = ps
    b, t = symbols('b t')
    parent = parametric_check(permute(co, LABEL['C']), ps)-units
    zero = (0, 2, 2)
    _, _, vs, mons = chart(134217471, zero)
    expected = {a: 1/a, c: a*c, d: 1/d, e: a*d*e/(a*d*e-e+1)}
    for phi in reindexings(134217471, zero, zero):
        sub = {v: cancel(co[phi[w]]*prod(co[phi[j]]**n for j, n in powers))
               for v, (w, powers) in zip(vs, mons)}
        if all(cancel(sub[p]-expected[p]) == 0 for p in ps):
            assert all(cancel(co[w].subs(expected, simultaneous=True)-sub[v]) == 0
                       for v, (w, _) in zip(vs, mons))
            break
    else:
        raise AssertionError('First-family stabilizer not recovered')
    transported = restriction(parent, expected, units)
    guards = units | {d-1, a*c*d-c-d+1}
    if name == 'linear':
        f = a*a*c*d-a*c-a*d-a*e+2*a+e-1
        change = {}
        coords = (a, d, e)
    else:
        f = (a*a*c*d*d*e-a*a*c*d-a*a*d*e+a*a*d-2*a*c*d*e+a*c*d+a*c
             -a*d*d*e+2*a*d*e+a*e-2*a+c*e-c+d*e-d-2*e+2)
        change = {d: (t+1)/(a*t), e: (a*b-1)*t}
        back = {b: (a*d*e-e+1)/a, t: 1/(a*d-1)}
        for p, val in change.items():
            assert cancel(val.subs(back, simultaneous=True)-p) == 0
        for p, val in back.items():
            assert cancel(val.subs(change, simultaneous=True)-p) == 0
            assert factors(val) <= guards
        guards = local_units(guards, change) | {b, t}
        coords = (a, b, t)
    assert f in parent
    graph, = restriction({f}, change, guards)
    polynomial = Poly(graph, c)
    assert polynomial.degree() == 1
    if factors(polynomial.LC())-guards:
        assert not (factors(polynomial.TC())-guards)
    sub = {c: cancel(-polynomial.TC()/polynomial.LC())}
    coeff = apply_sub(apply_sub(co, change), sub)
    locations = {p: cancel(change.get(p, p).subs(sub, simultaneous=True)) for p in ps}
    local = local_units(guards, sub) | factors(polynomial.LC())
    cuts = restriction(transported, locations, local)
    assert len(cuts) == (8 if name == 'linear' else 7)
    return coords, coeff, local, sorted(cuts, key=str)


def parametrize(g, params, guards):
    a, b, t = params
    if (g-((b*(a-1))**2+b*(a-1)-b+1)).expand():
        return surface_chart(g, params, guards)
    w = symbols('w')
    inverse = b*(a-1)
    sub = {b: w*w+w+1, a: (w+1)**2/(w*w+w+1)}
    assert cancel(g.subs(sub, simultaneous=True)) == 0
    assert cancel(inverse.subs(sub, simultaneous=True)-w) == 0
    assert cancel(sub[b].subs(w, inverse)-b) == g.expand()
    numerator = cancel(sub[a].subs(w, inverse)-a).as_numer_denom()[0]
    assert Poly(numerator, a, b).rem(Poly(g, a, b)).is_zero and b in guards
    return sub, (w, t)


def graph_count(name, index, coeff, vs, local):
    charts = [fs-{-u for u in local} for fs in three_charts(coeff, vs, local)]
    curves = set.intersection(*charts)
    rest = [fs-curves for fs in charts]
    base = component_bases(vs, rest, local | curves) if all(rest) else []
    off = points(coeff, vs, local, base)
    exceptional = pencils = fibres = 0
    for j, curve in enumerate(sorted(curves, key=str)):
        if (name, index, j) in EXTRA:
            label, sample = EXTRA[name, index, j]
            fs = local_cuts(coeff, vs, label, sample)-local-{-u for u in local}
            assert curve not in fs and -curve not in fs
            base = component_bases(vs, [{curve}, fs, fs], local) if fs else []
            exceptional += points(coeff, vs, local, base)
        else:
            fibres += entire_pencil(coeff, local, curve_parameters(curve, vs))
            pencils += 1
    result = off, exceptional, pencils, fibres
    print('PASS graph', name, index, result, flush=True)
    return result

EXPECTED = {
    'quadratic': [(3, 0, 2, 13), (5, 10, 0, 0), (3, 18, 0, 0), (0, 16, 1, 3),
                  (1, 0, 3, 55), (3, 4, 1, 22), (0, 0, 1, 2)],
    'linear': [(8, 0, 2, 31), (5, 7, 2, 9), (5, 0, 1, 18), (4, 0, 2, 22),
               (2, 11, 3, 38), (5, 2, 3, 50), (1, 0, 0, 0)],
}


def verify_graphs(name):
    params, coeff, guards, cuts = normalized_parent(name)
    if name == 'linear':
        e = params[-1]
        assert cuts[-1] == e**2-3*e+3
        cuts = cuts[:-1]
    results = []
    for i, g in enumerate(cuts):
        sub, vs = parametrize(g, params, guards)
        section = apply_sub(coeff, sub)
        local = local_units(guards, sub)
        results.append(graph_count(name, i, section, vs, local))
    assert results == EXPECTED[name], results
    return params, coeff, guards


def verify():
    for name in ('quadratic', 'linear'):
        verify_graphs(name)
    co, ps, units = family(0)
    a, c, d, e = ps
    central = a*c*d-c-d+1
    old = a*a*c*d*e-a*c*e-a*d*e-a*e+2*a+2*e-2
    junction = a*a*c*d-a*c-a*d*e+a+e-1
    linear = a*a*c*d-a*c-a*d-a*e+2*a+e-1
    b = (a*d*e-e+1)/a
    t = 1/(a*d-1)
    quadratic, = factors(b*t*(a-1)**2-a*c*(b-1)+b-a)-units
    orbits = verify_orbits(co, ps, units, {central, d-1, old})
    representatives = (junction, linear, quadratic)
    indices = [next(i for i, orbit in enumerate(orbits) if f in orbit) for f in representatives]
    assert len(set(indices)) == 3
    print('PASS seven graphs on each divisor and distinct orbits', flush=True)


if __name__ == '__main__':
    verify()
