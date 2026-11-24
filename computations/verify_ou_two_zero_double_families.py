#!/usr/bin/env python3
"""All four d=1 three-parameter families on the two-zero support 134217455."""
from argparse import ArgumentParser

from sympy import Poly, cancel

from verify_ou_complement_cycles import face_dets
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check, verify_point
from verify_ou_three_zero_remaining import LABELS, central_generation, curve_parameters, t
from verify_ou_zero_supports import TRIPLES


def permute(c, labels):
    return {w:c[tuple(p[i] for p,i in zip(labels,w))] for w in TRIPLES}


def curve_chart(f, params):
    x,y = params
    cusp = 3*x**2*y-x**2+3*x*y**2-6*x*y+2*x+(y-1)**3
    conic = x**2+x*y-2*x+y**2-2*y+1
    if f == cusp.expand():
        chart = {x:1/(1+t**3),y:(t**3+t**2)/(1+t**3)}
    elif f == conic:
        chart = {x:1/(t**2+t+1),y:(t+1)**2/(t**2+t+1)}
    else:
        return curve_parameters(f,params)
    assert cancel(f.subs(chart,simultaneous=True))==0
    return chart


def verify_family(name):
    _,template,(a,b,c,d,e),_ = chart(134217455,(0,1,1))
    # Interchanging the zeros exchanges the two product equations on d=1.
    swap = ((0,1,2),(2,1,0),(2,1,0))
    assert permute(apply_sub(template,{d:1}),swap)==apply_sub(template,{a:b,b:a,c:e,d:1,e:c})
    if name=='product':
        coeff = apply_sub(template,{d:1,a:(b+c*e-e)/b})
        params = (b,c,e)
        surfaces = [(b+c-1,b), (b-1,b), (b-(c-1)**2,b),
                    (b**2+(c-1)*e,e), (b**2-(c-1)**2*e,e)]
    elif name=='sum':
        coeff = apply_sub(template,{d:1,e:1-a})
        params = (a,b,c)
        surfaces = [(2*a+b-c-1,c),(2*a*b+a*c-a+b**2-2*b-c+1,c),
                    (a*b+2*a*c-a-c**2,a),(a*b+a*c-c,a),
                    (b+c-1,b),(b-1,b),(c-1,c)]
    else:
        assert name=='linear'
        coeff = apply_sub(template,{d:1,a:b+c-e})
        params = (b,c,e)
        surfaces = [(2*b+2*c-e-1,e),(2*b+c-1,c),(2*b+c-e-1,e),
                    (2*b+c-e-2,e),(b+2*c-1,b),(b+c,b),(b+c-1,b),
                    (b+c-e-1,e),(b-1,b),(c+e,c),(c-1,c)]
    units = set().union(*(factors(v) for v in coeff.values() if v))
    units.update(set().union(*(factors(v) for v in face_dets(coeff))))
    cuts = parametric_check(permute(coeff,LABELS['A']),params)-units
    assert cuts == {f.expand() for f,_ in surfaces}
    counts = []
    for f,pivot in surfaces:
        poly = Poly(f,pivot)
        # An invertible constant coefficient also forces a nonzero
        # leading coefficient at every point of this graph.
        assert poly.degree()==1 and (factors(poly.LC())<=units or factors(poly.TC())<=units)
        value = cancel(-poly.TC()/poly.LC())
        co = apply_sub(coeff,{pivot:value})
        ps = tuple(z for z in params if z!=pivot)
        local = set().union(*(factors(cancel(v.subs(pivot,value))) for v in units))
        local.update(factors(value))
        reindexed = permute(co,LABELS['A'])
        if f==b+c-1:
            curves = central_generation(reindexed,ps)-local
        else:
            curves = parametric_check(reindexed,ps)-local
        points = 0
        for g in sorted(curves,key=str):
            change = curve_chart(g,ps)
            at = apply_sub(co,change)
            required = [cancel(v.subs(change,simultaneous=True)) for v in local]
            assert all(required)
            try:
                failed = parametric_check(at,(t,))
            except AssertionError:
                failed = central_generation(at)
            exceptional = failed-set().union(*(factors(v) for v in required))
            for h in sorted(exceptional,key=str):
                witness = verify_point(at,required,h)
                print('PASS closed factor',g,h,witness,flush=True)
                points += 1
        counts.append((len(curves),points))
        print('PASS entire divisor',name,f,counts[-1],flush=True)
    expected = {'product':(5,6,38),'sum':(7,18,24),'linear':(11,57,278)}
    assert (len(counts),sum(n for n,_ in counts),sum(n for _,n in counts))==expected[name]
    print('PASS entire family',name,counts,flush=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    families = ('product','sum','linear')
    parser.add_argument('--family',choices=families)
    selected = parser.parse_args().family
    for name in ((selected,) if selected else families):
        verify_family(name)
