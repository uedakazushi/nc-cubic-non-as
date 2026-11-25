#!/usr/bin/env python3
"""All eight common divisors of the two parent charts in the second one-zero family."""
from sympy import Poly, QQ, cancel, groebner, prod

from verify_ou_one_zero_closed_divisors import LABEL, full_graph, points, three_charts
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart, reindexings
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import curve_parameters, t
from verify_ou_two_zero_closed_strata import entire_pencil, local_cuts, local_units, primary_intersection
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, surface_chart


CROSSING_CHARTS = (
    (((1,0,2),(0,1,2),(2,0,1)), (7,(2,5,6))),
    (((1,2,0),(2,0,1),(1,2,0)), (7,(6,2,3))),
    (((1,0,2),(0,1,2),(2,1,0)), (7,(2,5,6))),
    (LABEL['D'], (7,(6,4,4))),
    (LABEL['E'], (7,(3,2,6))),
)


def divisor(a,c,d,h):
    return (a-1)*c**2*(a*d-1) + a*(d-1)*(a*c*d-c-d+1)*h


def parent_orbits(co, ps, units):
    a,c,d,h = ps
    zero = (0,2,2)
    mask = 134217471
    _,_,variables,mons = chart(mask,zero)
    A,B,C,D,E,H = variables
    checks = (A*B-A*D*E+E-1, A*C*D*E-A*C-A*D*H+A*H-C*E+C)
    formulas = (divisor(A,C,D,H), A*C-1, A*C-C+D-1, A*C*D-A*D+A-C,
                A*C**2*D-C**2-D+1, A*C*D+A*D**2*H-A*D*H-C,
                A*C**2+A*C*D-A*C-C**2-D+1)
    orbits = [set() for _ in formulas]
    for phi in reindexings(mask,zero,zero):
        vals = {v:cancel(co[phi[w]]*prod(co[phi[k]]**i for k,i in pw))
                for v,(w,pw) in zip(variables,mons)}
        if any(cancel(eq.subs(vals,simultaneous=True)) for eq in checks):
            continue
        for formula,orbit in zip(formulas,orbits):
            image = cancel(formula.subs(vals,simultaneous=True))
            fac = factors(image)-units
            assert len(fac)==1 and factors(image.as_numer_denom()[1])<=units
            orbit.update(fac)
    parent = [parametric_check(permute(co,LABEL[s]),ps)-units for s in 'AC']
    common = parent[0]&parent[1]
    assert len(common)==8 and list(map(len,orbits))==[4,4,8,2,4,1,8]
    settled = orbits[0]|orbits[1]|{a*c*d-c-d+1,a*c*d-c+d-1}
    assert common<=settled
    remaining = parent[1]-settled
    assert remaining<=set.union(*orbits[2:])
    assert [len(remaining&orb) for orb in orbits[2:]]==[2,1,4,1,2]
    print('PASS both four-element orbits and all eight parent factors',flush=True)


def surfaces(co, ps, units, shared):
    totals = [0,0,0]
    for g in sorted(shared,key=str):
        sub,vs = surface_chart(g,ps,units)
        coeff = apply_sub(co,sub)
        local = local_units(units,sub)
        charts = three_charts(coeff,vs,local)
        common = set.intersection(*charts)
        rest = [F-common for F in charts]
        if all(rest):
            totals[0] += points(coeff,vs,local,component_bases(vs,rest,local|common))
        for f in sorted(common,key=str):
            totals[1] += 1
            totals[2] += entire_pencil(coeff,local,curve_parameters(f,vs))
    assert totals==[8,3,29],totals
    print('PASS common surfaces',totals,flush=True)


def crossings(co, ps, units, shared, charts):
    a,c,d = ps
    guards = units|shared|factors(a*c*d-c-d+1)
    comps = primary_intersection(ps,guards,sorted([F-guards for F in charts],key=len))
    assert not comps[0] and len(comps[1])==25
    remaining = {b:groebner(b,*ps,order='grevlex',domain=QQ) for b in comps[1]}
    counts = []
    for label,sample in CROSSING_CHARTS:
        cuts = local_cuts(co,ps,label,sample)-guards
        done = [b for b,gb in remaining.items() if all(gb.reduce(f)[1] for f in cuts)]
        npoints = 0
        for b in done:
            residual = primary_intersection(ps,guards,[cuts],initial=b) if cuts else {0:[],1:[]}
            assert not residual[1]
            npoints += points(co,ps,guards,residual[0])
            del remaining[b]
        counts.append((len(done),npoints))
    assert [n for n,_ in counts]==[4,5,9,5,1],counts
    b, = remaining
    full_graph(co,ps,guards,b,{a:t**2,c:t**2,d:1/(t**2-t+1)},c+1-1/d)
    print('PASS 25 crossing curves',counts,flush=True)


def verify():
    co,ps,units = family(1)
    parent_orbits(co,ps,units)
    a,c,d,h = ps
    poly = Poly(divisor(*ps),h)
    assert poly.degree()==1 and factors(poly.TC())<=units
    # The constant coefficient is invertible, hence the leading coefficient
    # never vanishes on the divisor and no vertical fibre is lost.
    sub = {h:cancel(-poly.TC()/poly.LC())}
    coeff = apply_sub(co,sub)
    local = local_units(units,sub)|factors(poly.LC())
    vs = (a,c,d)
    charts = three_charts(coeff,vs,local)
    shared = set.intersection(*charts)
    assert len(shared)==3
    surfaces(coeff,vs,local,shared)
    crossings(coeff,vs,local,shared,charts)


if __name__=='__main__':
    verify()
