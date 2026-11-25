#!/usr/bin/env python3
"""Test the generic points of the remaining one-zero surfaces over Q."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import product

from sympy import Poly, QQ, S, cancel, count_ops, groebner, symbols

from verify_ou_complement_cycles import face_dets
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import LABELS
from verify_ou_two_zero_closed_strata import local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, open_witness


def family(index):
    _,template,(a,b,c,d,e,h),_ = chart(134217471,(0,2,2))
    if index==0:
        sub = {b:(a*d*e-e+1)/a,h:(a*d*e-a-e+1)/(a*(a*d-1))}
        params = a,c,d,e
    else:
        sub = {b:(c+d*h-h)/c,e:(a*c+a*d*h-a*h-c)/(a*c*d-c)}
        params = a,c,d,h
    co = apply_sub(template,sub)
    units = set().union(*(factors(v) for v in co.values() if v),
                        *(factors(v) for v in face_dets(co)))
    return co,params,units


def graph_options(f, params):
    for p in params:
        poly = Poly(f,p)
        if poly.degree()==1:
            yield p,cancel(-poly.TC()/poly.LC()),poly.LC(),poly.TC()


def surfaces(co, params, units, cuts):
    for f in sorted(cuts-units,key=str):
        print('PASS surface',f,open_witness(co,params,(f,),units),flush=True)
    return len(cuts-units)


def verticals(poly, pivot, params, units, accepts):
    # Off the graph the pivot is free. A two-dimensional component must be
    # the whole cylinder above a codimension-two component of (LC,TC).
    pol = Poly(poly,pivot)
    left,right = factors(pol.LC()),factors(pol.TC())
    assert not left&right
    if not left or not right:
        return
    rest = tuple(z for z in params if z!=pivot)
    independent = {z for z in units if not z.has(pivot)}|{S.One}
    for basis in component_bases(rest,[left,right,right],independent):
        gb = groebner(basis,*params,order='grevlex',domain=QQ)
        if all(gb.reduce(z)[1] for z in units):
            assert not accepts(lambda z: gb.reduce(z)[1]==0),(poly,basis)


def common_graph(co, params, units, f):
    choices = []
    for p,value,lead,_ in graph_options(f,params):
        coeff = apply_sub(co,{p:value})
        cost = (bool(factors(lead)-units),sum(count_ops(z) for z in coeff.values()),
                len(str(value)),str(p))
        choices.append((cost,p,value,lead,coeff))
    _,p,value,lead,coeff = min(choices)
    rest = tuple(z for z in params if z!=p)
    local = local_units(units,{p:value})|factors(lead)
    cuts = parametric_check(permute(coeff,LABELS['C']),rest)
    count = surfaces(coeff,rest,local,cuts)
    verticals(f,p,params,units,lambda contains: True)
    return count


def first_changed_graph(co, params, units, f):
    # x=h and y=b. The torus units ad-1, h and e make this change invertible.
    a,c,d,e = params
    x,y = symbols('x y')
    forward = {d:(x+y-1)/(a*x),e:(a*y-1)*x/(y-1)}
    backward = {x:(a*d*e-a-e+1)/(a*(a*d-1)),y:(a*d*e-e+1)/a}
    for z,v in forward.items():
        assert cancel(v.subs(backward,simultaneous=True)-z)==0
    for z,v in backward.items():
        assert cancel(v.subs(forward,simultaneous=True)-z)==0
    assert all(factors(z)<=units for z in (a,backward[x],backward[y]-1,a*d-1))
    local = local_units(units,forward)|{x,a}
    q, = factors(cancel(f.subs(forward,simultaneous=True)))-local
    _,value,lead,_ = next(graph_options(q,(c,)))
    assert factors(lead)<=local
    new = apply_sub(apply_sub(co,forward),{c:value})
    local = local_units(local,{c:value})|factors(lead)
    # The previously resolved central divisor is a component of one cut here.
    special = cancel((a*c*d-c-d+1).subs(forward,simultaneous=True).subs(c,value))
    assert factors(special)-local=={2*a*x*y-a*y-2*x-y+2}
    local |= factors(special)
    rest = (a,y,x)
    return surfaces(new,rest,local,parametric_check(permute(new,LABELS['C']),rest))


def nonlinear(co, params, units, f, g):
    choices = []
    for x in (params[1],params[-1]):
        F,G = Poly(f,x),Poly(g,x)
        assert F.degree()==G.degree()==2
        remainder = (G.LC()*f-F.LC()*g).expand()
        assert remainder
        cuts = factors(remainder)-units
        if cuts and all(Poly(q,x).degree()==1 for q in cuts):
            choices.append((sum(len(str(q)) for q in cuts),str(x),x,cuts))
    _,_,x,graphs = min(choices)
    count = 0
    rest = tuple(z for z in params if z!=x)
    for q in graphs:
        _,value,lead,_ = next(graph_options(q,(x,)))
        sub = {x:value}
        local = local_units(units,sub)|factors(lead)
        polys = [cancel(z.subs(sub)).as_numer_denom()[0] for z in (f,g)]
        nonzero = [factors(z) for z in polys if z]
        assert nonzero
        count += surfaces(apply_sub(co,sub),rest,local,set.intersection(*nonzero))
        # A unit coefficient of either quadratic already excludes a whole
        # vertical cylinder. Otherwise test its generic point by reduction.
        if not any(cf and factors(cf)<=units for z in (f,g) for cf in Poly(z,x).all_coeffs()):
            verticals(q,x,params,units,lambda contains: contains(f) and contains(g))
    return count


def crossings(co, params, units, common, left, right, part):
    @lru_cache(None)
    def preferred(f):
        choices = [(bool(factors(lead)-units and factors(const)-units),
                    len(str(value)),str(p),p,value,lead)
                   for p,value,lead,const in graph_options(f,params)]
        return min(choices) if choices else None

    guards = units|common
    groups,quadratics = {},[]
    for f,g in product(sorted(left-common,key=str),sorted(right-common,key=str)):
        a,b = preferred(f),preferred(g)
        if not a and not b:
            quadratics.append((f,g))
            continue
        q,other = (f,g) if a and (not b or a[:3]<=b[:3]) else (g,f)
        groups.setdefault(q,set()).add(other)
    count = 0
    if part!='quadratic':
        for q,others in sorted(groups.items(),key=lambda pair:str(pair[0])):
            _,_,_,p,value,lead = preferred(q)
            sub = {p:value}
            local = local_units(guards,sub)|factors(lead)
            rest = tuple(z for z in params if z!=p)
            cuts = set()
            for polynomial in others:
                at = cancel(polynomial.subs(sub))
                assert at
                cuts.update(factors(at))
            count += surfaces(apply_sub(co,sub),rest,local,cuts)
            verticals(q,p,params,guards,lambda contains: any(map(contains,others)))
    if part!='linear':
        for f,g in quadratics:
            count += nonlinear(co,params,guards,f,g)
    return len(groups),len(quadratics),count


def verify(index, part):
    co,params,units = family(index)
    a,c,d,_ = params
    labels = ('C','A') if index==0 else ('A','C')
    left,right = [parametric_check(permute(co,LABELS[label]),params)-units for label in labels]
    common = left&right
    assert (len(left),len(right),len(common))==((12,29,7),(18,19,8))[index]
    solved = a*c*d-c-d+1
    assert solved in common
    if part in ('all','common'):
        counts = []
        for i,f in enumerate(sorted(common,key=str)):
            if f==solved:
                continue
            check = first_changed_graph if index==0 and i==4 else common_graph
            counts.append(check(co,params,units,f))
        assert sum(counts)==(45,77)[index]
        print('PASS common divisors',index,counts,flush=True)
    if part!='common':
        groups,pairs,count = crossings(co,params,units,common,left,right,part)
        assert (groups,pairs)==((18,8),(20,0))[index]
        linear,quadratic = ((88,8),(98,0))[index]
        expected = {'linear':linear,'quadratic':quadratic,'all':linear+quadratic}
        assert count==expected[part]
        print('PASS crossing graphs',index,groups,pairs,count,flush=True)


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--family',type=int,choices=(0,1))
    parser.add_argument('--part',choices=('all','common','linear','quadratic'),default='all')
    args = parser.parse_args()
    for index in ((0,1) if args.family is None else (args.family,)):
        verify(index,args.part)
