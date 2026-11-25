#!/usr/bin/env python3
"""All points on three divisors of the one-zero families, and their relabellings."""
from argparse import ArgumentParser
from sympy import Poly, QQ, S, cancel, groebner, prod

from verify_ou_one_zero_residual_surfaces import family
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart, reindexings
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import LABELS, curve_parameters, t
from verify_ou_two_zero_closed_strata import entire_pencil, local_cuts, local_units, primary_intersection
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, open_witness, surface_chart


NAMES = ('d=1','ac=1','acd-c+d-1=0')
LABEL = dict(LABELS,D=((2,1,0),(2,0,1),(0,1,2)),E=((2,1,0),(0,1,2),(2,0,1)))
SAMPLES = (
    [('A',7,(3,6,4)),('C',7,(2,2,3)),('B',7,(2,6,2)),('D',7,(4,5,6)),('E',5,(4,4,4))],
    [('A',7,(2,6,5)),('C',11,(4,6,6)),('B',7,(2,4,6))],
    [])


def three_charts(co, ps, units):
    return [parametric_check(permute(co,LABEL[name]),ps)-units for name in 'ACB']


def points(co, ps, units, bases):
    for b in bases:
        print('PASS point',b,open_witness(co,ps,b,units),flush=True)
    return len(bases)


def common_surfaces(index, co, ps, units, common):
    counts = [0,0,0]
    previous = None
    for g in sorted(common,key=str):
        sub,vs = surface_chart(g,ps,units)
        coeff = apply_sub(co,sub)
        local = local_units(units,sub)
        charts = three_charts(coeff,vs,local)
        shared = set.intersection(*charts)
        rest = [F-shared for F in charts]
        if all(rest):
            counts[0] += points(coeff,vs,local,component_bases(vs,rest,local|shared))
        x,y = vs
        doubled = x*x*y*y-x*y**3+x*y+2*y**3-5*y*y+4*y-1
        for h in sorted(shared,key=str):
            if index==0 and h==doubled:
                if previous is None:
                    cuts = local_cuts(coeff,vs,LABEL['A'],(7,(2,2)))-local
                    assert h not in cuts
                    counts[0] += points(coeff,vs,local,component_bases(vs,[{h},cuts,cuts],local))
                    previous = sub
                else:
                    # The two graphs agree on this curve. The first graph's
                    # only denominator is invertible throughout the source torus.
                    assert set(previous)==set(sub)
                    for p,val in previous.items():
                        assert factors(val.as_numer_denom()[1])<=local
                        num = cancel(val-sub[p]).as_numer_denom()[0]
                        assert Poly(num,*vs).rem(Poly(h,*vs)).is_zero
                continue
            counts[1] += 1
            counts[2] += entire_pencil(coeff,local,curve_parameters(h,vs))
    print('PASS common surfaces',NAMES[index],counts,flush=True)
    assert counts==([27,12,54],[8,2,10],[2,4,22])[index],counts


def full_graph(co, ps, units, basis, sub, inverse):
    gb = groebner(basis,*ps,order='grevlex',domain=QQ)
    assert all(cancel(v.subs(sub,simultaneous=True))==0 for v in basis)
    assert cancel(inverse.subs(sub,simultaneous=True)-t)==0
    for p,v in sub.items():
        back = cancel(v.subs(t,inverse))
        assert gb.reduce(cancel(back-p).as_numer_denom()[0])[1]==0
        for f in factors(back.as_numer_denom()[1])-units:
            empty = primary_intersection(ps,units,[{f}],initial=basis)
            assert not any(empty.values())
    print('PASS crossing pencil',basis,entire_pencil(co,units,sub),flush=True)


def crossings(index, co, ps, units, charts, common):
    a,c,d = (co[w] for w in ((0,2,2),(1,2,0),(2,0,2)))
    guards = units|common|factors(cancel(a*c*d-c-d+1))
    cuts = sorted([F-guards for F in charts],key=len)
    components = primary_intersection(ps,guards,cuts)
    assert (len(components[0]),len(components[1]))==((63,17),(0,8),(0,3))[index]
    points(co,ps,guards,components[0])
    remaining = {b:groebner(b,*ps,order='grevlex',domain=QQ) for b in components[1]}
    for name,prime,point in SAMPLES[index]:
        fs = local_cuts(co,ps,LABEL[name],(prime,point))-guards
        covered = [b for b,G in remaining.items() if all(G.reduce(z)[1] for z in fs)]
        for b in covered:
            pts = primary_intersection(ps,guards,[fs],initial=b) if fs else {0:[],1:[]}
            assert not pts[1]
            points(co,ps,guards,pts[0])
            del remaining[b]
    if index==0:
        a,c,e = ps
        b, = remaining
        full_graph(co,ps,guards,b,{a:t,c:1/t,e:(t-2)/(t-1)},a)
    elif index==1:
        assert not remaining
    else:
        c,d,h = ps
        images = [((t+1)**2/(2*t),-t**2),
                  (t*(t+1)/(t-1),t*(1-t)/(t+1)),((1-t**2)/2,S.One)]
        for b,(de,hh) in zip(components[1],images):
            full_graph(co,ps,guards,b,{c:t,d:de,h:hh},c)


def symmetries():
    co,ps,u = family(1)
    a,c,d,h = ps
    zero = (0,2,2)
    mask = 134217471
    _,_,variables,mons = chart(mask,zero)
    A,B,C,D,E,H = variables
    orbit = set()
    for phi in reindexings(mask,zero,zero):
        vals = {v:cancel(co[phi[w]]*prod(co[phi[k]]**i for k,i in pw))
                for v,(w,pw) in zip(variables,mons)}
        eqs = (A*B-A*D*E+E-1,A*C*D*E-A*C-A*D*H+A*H-C*E+C)
        if any(cancel(z.subs(vals,simultaneous=True)) for z in eqs):
            continue
        orbit.update(factors(cancel(vals[A]*vals[C]-1))-u)
    assert orbit=={a*c-1,c-1,a*a*c*d*h-a*c*h-a*c-a*d*h+a*h+c,
                  a*c*d*h-a*c-a*d*h+a*h-c*h+c}
    print('PASS four divisors in the relabelling orbit',flush=True)


def verify(index):
    co,ps,units = family(int(index>0))
    a,c,d,_ = ps
    sub = ({d:S.One},{a:1/c},{a:(c-d+1)/(c*d)})[index]
    vs = tuple(p for p in ps if p not in sub)
    co = apply_sub(co,sub)
    units = local_units(units,sub)
    charts = three_charts(co,vs,units)
    common = set.intersection(*charts)
    assert len(common)==(5,3,5)[index]
    common_surfaces(index,co,vs,units,common)
    crossings(index,co,vs,units,charts,common)


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--divisor',type=int,choices=range(3))
    args = parser.parse_args()
    for index in (range(3) if args.divisor is None else (args.divisor,)):
        verify(index)
    symmetries()
