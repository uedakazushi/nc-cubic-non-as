#!/usr/bin/env python3
"""Exhaustive tests on further orbit representatives in the second one-zero family."""
from sympy import Poly, QQ, cancel, groebner

from verify_ou_one_zero_closed_divisors import points, three_charts
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_one_zero_second_common import parent_orbits
from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import apply_sub
from verify_ou_two_zero_closed_strata import local_cuts, local_units, primary_intersection


def verify_graph(co, ps, guards, equation, pivot, samples, expected):
    graph = Poly(equation,pivot)
    assert graph.degree()==1 and factors(graph.LC())<=guards
    sub = {pivot:cancel(-graph.TC()/graph.LC())}
    params = tuple(v for v in ps if v!=pivot)
    local = local_units(guards,sub)
    coeff = apply_sub(co,sub)
    charts = three_charts(coeff,params,local)
    assert not set.intersection(*charts)
    base = primary_intersection(params,local,sorted(charts,key=len))
    assert not base[0]
    remaining = {b:groebner(b,*params,order='grevlex',domain=QQ) for b in base[1]}
    counts = []
    for label,sample in samples:
        cuts = local_cuts(coeff,params,label,sample)-local
        covered = [b for b,G in remaining.items() if all(G.reduce(f)[1] for f in cuts)]
        count = 0
        for b in covered:
            if cuts:
                residual = primary_intersection(params,local,[cuts],initial=b)
                assert not residual[1]
                count += points(coeff,params,local,residual[0])
            del remaining[b]
        counts.append((len(covered),count))
    assert not remaining and counts==expected, (remaining,counts)
    print('PASS divisor',equation,counts,flush=True)


def verify():
    co,ps,units = family(1)
    guards = units|parent_orbits(co,ps,units)
    a,c,d,h = ps
    samples = [
        (((1,2,0),(2,0,1),(1,2,0)),(7,(3,6,2))),
        (((1,0,2),(0,1,2),(2,1,0)),(7,(5,2,4))),
        (((2,1,0),(2,0,1),(0,1,2)),(7,(4,5,2))),
    ]
    verify_graph(co,ps,guards,a*c*d-a*d+a-c,d,samples,[(9,17),(5,0),(3,6)])
    samples = [
        (((1,0,2),(0,1,2),(2,0,1)),(7,(3,3,1))),
        (((1,2,0),(2,0,1),(1,2,0)),(11,(2,7,3))),
        (((1,0,2),(0,1,2),(2,1,0)),(7,(3,3,1))),
    ]
    verify_graph(co,ps,guards,a*c**2*d-c**2-d+1,a,samples,[(2,8),(3,6),(3,0)])


if __name__=='__main__':
    verify()
