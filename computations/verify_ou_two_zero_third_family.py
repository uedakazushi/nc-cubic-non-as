#!/usr/bin/env python3
"""Full recognition on the third two-zero family, away from the shared surfaces."""
from sympy import Poly, cancel, groebner, QQ

from verify_ou_two_zero_closed_strata import (
    LABELS, apply_sub, component_bases, entire_pencil, family, local_cuts, local_units,
    open_witness, parametric_check, permute, primary_intersection, t,
)


def five_chart_intersection(coeff,params,units):
    labels = [LABELS[k] for k in 'ACB']+[
        ((2,1,0),(2,0,1),(0,1,2)),((2,1,0),(0,1,2),(2,0,1))]
    charts = [parametric_check(permute(coeff,p),params)-units for p in labels]
    common = set.intersection(*charts[:3])
    assert len(common)==5
    charts = sorted((s-common for s in charts),key=len)
    return primary_intersection(params,units|common,charts)


def verify():
    coeff,params,units = family(2)
    c,d,e = params
    components = five_chart_intersection(coeff,params,units)
    assert (len(components[0]),len(components[1]))==(10,6)
    for basis in components[0]:
        print('PASS isolated point',basis,open_witness(coeff,params,basis,units),flush=True)

    rational = {
        2*d*e-d-e:({c:t,d:1/(1-t),e:1/(1+t)},c,(1-c)*(1+c)),
        c+d-1:({c:1-t,d:t,e:t**3/(2*t-1)},d,2*d-1),
        d*d-2*d-e+2:({c:-2*(t-1)**3/t**2,d:t,e:t*t-2*t+2},d,d)}
    quadratic = {
        d*d+d-1:(1-d*e,(5,(2,2)),4),
        d*d-2*d+2:((d*e-3*e+2)/2,(5,(4,2)),8)}
    found = set()
    for basis in components[1]:
        eq, = set(basis)&(set(rational)|set(quadratic)|{d+1})
        assert eq not in found
        found.add(eq)
        if eq==d+1:
            # The full d=-1 slice is covered by verify_ou_two_zero_residual_curves.py.
            continue
        if eq in rational:
            sub,inverse,denom = rational[eq]
            assert all(cancel(h.subs(sub,simultaneous=True))==0 for h in basis)
            gb = groebner(basis,*params,order='grevlex',domain=QQ)
            assert all(gb.reduce(cancel(v.subs(t,inverse)-z).as_numer_denom()[0])[1]==0
                       for z,v in sub.items())
            exceptional = primary_intersection(params,units,[{denom}],initial=basis)
            assert not exceptional[1]
            for pt in exceptional[0]:
                print('PASS vertical fibre',pt,open_witness(coeff,params,pt,units),flush=True)
            print('PASS rational curve',eq,entire_pencil(coeff,units,sub),flush=True)
            continue
        value,sample,count = quadratic[eq]
        gb = groebner(basis,*params,order='grevlex',domain=QQ)
        assert gb.reduce(cancel(c-value).as_numer_denom()[0])[1]==0
        assert all(Poly(cancel(h.subs(c,value)).as_numer_denom()[0],d,e).rem(Poly(eq,d,e)).is_zero
                   for h in basis)
        local = local_units(units,{c:value})
        co = apply_sub(coeff,{c:value})
        cuts = local_cuts(co,(d,e),LABELS['A'],sample)-local
        assert cuts and eq not in cuts
        points = component_bases((d,e),[{eq},cuts,cuts],local)
        assert len(points)==count
        for pt in points:
            print('PASS quadratic pencil point',eq,pt,open_witness(co,(d,e),pt,local),flush=True)
    assert len(found)==6
    print('PASS entire third two-zero family',flush=True)


if __name__=='__main__':
    verify()
