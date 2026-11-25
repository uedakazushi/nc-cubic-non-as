#!/usr/bin/env python3
"""One more whole divisor and a seven-representative cover of the first family."""
from sympy import Poly, cancel, prod, symbols

from verify_ou_one_zero_closed_divisors import LABEL, points, three_charts
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart, reindexings
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import curve_parameters
from verify_ou_two_zero_closed_strata import entire_pencil, local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, surface_chart


def verify_orbits(co,ps,units,settled):
    zero = (0,2,2)
    mask = 134217471
    _,_,vs,mons = chart(mask,zero)
    a,b,c,d,e,h = vs
    maps = []
    for phi in reindexings(mask,zero,zero):
        sub = {v:cancel(co[phi[w]]*prod(co[phi[t]]**i for t,i in pw))
               for v,(w,pw) in zip(vs,mons)}
        if all(cancel(f.subs(sub,simultaneous=True))==0 for f in
               (a*b-a*d*e+e-1,a*a*d*h-a*d*e-a*h+a+e-1)):
            maps.append(sub)
    assert len(maps)==2
    parent = parametric_check(permute(co,LABEL['C']),ps)-units
    assert len(parent)==12 and settled<=parent
    left = parent-settled
    sizes = []
    groups = []
    while left:
        f = min(left,key=str)
        orbit = set()
        for sub in maps:
            value = cancel(f.subs(sub,simultaneous=True))
            assert factors(value.as_numer_denom()[1])<=units
            image = factors(value)-units
            assert len(image)==1
            orbit.update(image)
        assert f in orbit and not orbit&settled
        groups.append(orbit&left)
        sizes.append(len(orbit&left))
        print('REP',f,flush=True)
        left -= orbit
    assert sorted(sizes)==[1,1,1,1,1,2,2]
    return groups


def verify():
    co,ps,units = family(0)
    a,c,d,e = ps
    b = symbols('b')
    central = a*c*d-c-d+1
    f = a*a*c*d*e-a*c*e-a*d*e-a*e+2*a+2*e-2
    verify_orbits(co,ps,units,{central,d-1,f})
    change = {d:(a*b+e-1)/(a*e)}
    inverse = {b:(a*d*e-e+1)/a}
    assert cancel(change[d].subs(inverse)-d)==0
    assert cancel(inverse[b].subs(change)-b)==0
    local = local_units(units|{central,d-1},change)|{b}
    graph, = factors(cancel(f.subs(change)))-local
    poly = Poly(graph,c)
    assert poly.degree()==1 and factors(poly.LC())<=local
    sub = {c:cancel(-poly.TC()/poly.LC())}
    coeff = apply_sub(apply_sub(co,change),sub)
    local = local_units(local,sub)|factors(poly.LC())
    params = (a,b,e)
    cuts = parametric_check(permute(coeff,LABEL['B']),params)-local-{-v for v in local}
    assert len(cuts)==6
    counts = []
    for z in sorted(cuts,key=str):
        chart_sub,variables = surface_chart(z,params,local)
        section = apply_sub(coeff,chart_sub)
        guards = local_units(local,chart_sub)
        charts = [fs-{-g for g in guards} for fs in three_charts(section,variables,guards)]
        curves = set.intersection(*charts)
        rest = [fs-curves for fs in charts]
        pts = component_bases(variables,rest,guards|curves) if all(rest) else []
        counts.append((points(section,variables,guards,pts),len(curves),
                       sum(entire_pencil(section,guards,curve_parameters(g,variables))
                           for g in sorted(curves,key=str))))
    assert counts==[(0,1,5),(0,0,0),(0,0,0),(0,1,2),(0,2,9),(1,0,0)],counts
    print('PASS first-family divisor and remaining seven representatives',counts,flush=True)


if __name__=='__main__':
    verify()
