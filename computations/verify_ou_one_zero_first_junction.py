#!/usr/bin/env python3
"""The c=h orbit on the first one-zero family, including every closed fibre."""

from sympy import Poly, cancel, symbols

from verify_ou_one_zero_closed_divisors import LABEL, points, three_charts
from verify_ou_one_zero_first_orbits import verify_orbits
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import curve_parameters
from verify_ou_two_zero_closed_strata import entire_pencil, local_cuts, local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, surface_chart


ROTATE = ((1, 2, 0), (2, 0, 1), (1, 2, 0))
EXTRA = {
    (0, 1): (ROTATE, (5, (2, 2))),
    (3, 0): (ROTATE, (5, (2, 3))),
    (3, 1): (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (5, (2, 3))),
    (5, 0): (ROTATE, (7, (3, 2))),
}


def junction_chart(g, params, guards):
    a, b, h = params
    special = (a*(b-1))**2+a*(b-1)-a+1
    if (g-special).expand():
        return surface_chart(g, params, guards)
    w = symbols('w')
    inverse = a*(b-1)
    sub = {a: w*w+w+1, b: (w+1)**2/(w*w+w+1)}
    assert cancel(g.subs(sub, simultaneous=True)) == 0
    assert cancel(inverse.subs(sub, simultaneous=True)-w) == 0
    assert cancel(sub[a].subs(w, inverse)-a) == g.expand()
    diff = cancel(sub[b].subs(w, inverse)-b).as_numer_denom()[0]
    assert Poly(diff, a, b).rem(Poly(g, a, b)).is_zero and a in guards
    return sub, (w, h)


def verify():
    co, ps, units = family(0)
    a, c, d, e = ps
    b, h = symbols('b h')
    central = a*c*d-c-d+1
    previous = a*a*c*d*e-a*c*e-a*d*e-a*e+2*a+2*e-2
    f = a*a*c*d-a*c-a*d*e+a+e-1
    orbits = verify_orbits(co, ps, units, {central, d-1, previous})
    orbit, = [o for o in orbits if f in o]
    assert len(orbit) == 2
    original_guards = units | {central, d-1}
    change = {d: (b+h-1)/(a*h), e: (a*b-1)*h/(b-1)}
    back = {b: (a*d*e-e+1)/a, h: (a*d*e-a-e+1)/(a*(a*d-1))}
    for v, val in change.items():
        assert cancel(val.subs(back, simultaneous=True)-v) == 0
    for v, val in back.items():
        assert cancel(val.subs(change, simultaneous=True)-v) == 0
    for invertible in (a, back[h], back[b], back[b]-1):
        assert factors(invertible) <= original_guards
    guards = local_units(original_guards, change) | {b, h}
    assert factors(cancel(f.subs(change, simultaneous=True)))-guards == {c-h}
    sub = {c: h}
    section = apply_sub(apply_sub(co, change), sub)
    guards = local_units(guards, sub)
    params = (a, b, h)
    cuts = parametric_check(permute(section, LABEL['B']), params)-guards-{-g for g in guards}
    assert len(cuts) == 12
    counts = []
    for i, g in enumerate(sorted(cuts, key=str)):
        sub, vs = junction_chart(g, params, guards)
        coeff = apply_sub(section, sub)
        local = local_units(guards, sub)
        charts = [fs-{-k for k in local} for fs in three_charts(coeff, vs, local)]
        curves = set.intersection(*charts)
        rest = [fs-curves for fs in charts]
        base = component_bases(vs, rest, local|curves) if all(rest) else []
        off = points(coeff, vs, local, base)
        exceptional = pencils = fibres = 0
        for j, curve in enumerate(sorted(curves, key=str)):
            if (i, j) in EXTRA:
                label, sample = EXTRA[i, j]
                fs = local_cuts(coeff, vs, label, sample)-local-{-k for k in local}
                assert curve not in fs and -curve not in fs
                base = component_bases(vs, [{curve}, fs, fs], local) if fs else []
                exceptional += points(coeff, vs, local, base)
            else:
                fibres += entire_pencil(coeff, local, curve_parameters(curve, vs))
                pencils += 1
        counts.append((off, exceptional, pencils, fibres))
        print('PASS section', i, counts[-1], flush=True)
    assert counts == [(6, 0, 1, 10), (10, 0, 0, 0), (2, 0, 2, 14), (7, 8, 1, 2),
                      (4, 0, 2, 18), (8, 1, 0, 0), (8, 0, 3, 25), (11, 0, 1, 14),
                      (0, 0, 0, 0), (4, 0, 2, 17), (0, 0, 2, 9), (3, 0, 1, 4)], counts
    print('PASS first-family c=h orbit', counts, flush=True)


if __name__ == '__main__':
    verify()
