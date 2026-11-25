#!/usr/bin/env python3
"""All nineteen common surfaces of the two-zero families and the entire second family."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import product
from math import prod
import os
import subprocess

from sympy import Poly, QQ, cancel, groebner, sympify

from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import PencilCertificate, apply_sub, parametric_check, verify_point
from verify_ou_three_zero_remaining import LABELS, central_generation, t
from verify_ou_two_zero_double_families import curve_chart, permute
from verify_ou_two_zero_residual_curves import component_bases, family, open_witness, surface_chart


def local_units(units, sub):
    values = [cancel(h.subs(sub,simultaneous=True)) for h in units]
    assert all(values)
    return set().union(*(factors(h) for h in values))


def entire_pencil(coeff, units, sub):
    at = apply_sub(coeff,sub)
    local = local_units(units,sub)
    try:
        failed = parametric_check(at,(t,))
    except AssertionError:
        failed = central_generation(at)
    for h in sorted(failed-local,key=str):
        print('PASS pencil fibre',h,verify_point(at,local,h),flush=True)
    return len(failed-local)


class LocalPencil(PencilCertificate):
    """Prefer exact pivots that remain units at the supplied prime and point."""
    def __init__(self, mask, coeff, field, sample=(7,(2,2))):
        self.prime,self.point = sample
        super().__init__(mask,coeff,field)

    @lru_cache(None)
    def residue(self, value):
        p = self.prime
        def evaluate(poly):
            return sum(int(k.numerator)*pow(int(k.denominator),-1,p)
                       * prod(pow(x,i,p) for x,i in zip(self.point,exponent))
                       for exponent,k in poly.items())%p
        a,b = evaluate(value.numer),evaluate(value.denom)
        return a*pow(b,-1,p)%p if b else None

    def eliminate(self, rows, width):
        todo = [dict(row) for row in rows]
        pivots = {}

        def clear(row, col, normalized):
            if col in row:
                multiple = row[col]
                for k,entry in normalized.items():
                    value = row.get(k,self.field.zero)-multiple*entry
                    if value:
                        row[k] = value
                    else:
                        row.pop(k,None)

        while any(todo):
            pos = next(((i,k) for i,row in enumerate(todo) for k in sorted(row)
                        if self.residue(row[k])),None)
            if pos is None:
                pos = next((i,min(row)) for i,row in enumerate(todo) if row)
            i,col = pos
            row = todo.pop(i)
            pivot = row[col]
            self.unit(pivot)
            normalized = {k:v/pivot for k,v in row.items()}
            pivots[col] = normalized
            for other in todo:
                clear(other,col,normalized)

        def reduce_row(row):
            result = dict(row)
            for col,val in pivots.items():
                clear(result,col,val)
            return result
        return [k for k in range(width) if k not in pivots],reduce_row


def local_cuts(coeff, params, label, sample):
    field = QQ.frac_field(*params)
    permuted = permute(coeff,label)
    cert = LocalPencil(0,{k:field.from_sympy(v) for k,v in permuted.items()},field,sample)
    cert.validate_paths()
    cert.check_generation()
    for r,i in product(range(3),repeat=2):
        cert.check_corner(r,i)
    return cert.exceptional


def quadratic_line(coeff, params, units):
    a,c,d = params
    sub = {a:(2*c-1)*d+1-c}
    co = apply_sub(coeff,sub)
    local = local_units(units,sub)
    q = 3*c*c-3*c+1
    cuts = local_cuts(co,(c,d),LABELS['A'],(7,(2,2)))-local
    assert cuts and q not in cuts
    points = component_bases((c,d),[{q},cuts,cuts],local)
    assert len(points)==5
    for basis in points:
        print('PASS quadratic pencil point',basis,open_witness(co,(c,d),basis,local),flush=True)


def surfaces(index, coeff, params, units, common):
    npoints = ncurves = nfibres = skipped = 0
    for h in sorted(common,key=str):
        sub,ps = surface_chart(h,params,units)
        co = apply_sub(coeff,sub)
        local = local_units(units,sub)
        charts = [parametric_check(permute(co,LABELS[label]),ps)-local for label in 'ACB']
        curves = set.intersection(*charts)
        rest = [fs-curves for fs in charts]
        points = component_bases(ps,rest,local|curves) if all(rest) else []
        npoints += len(points)
        for basis in points:
            print('PASS common surface point',index,h,basis,open_witness(co,ps,basis,local),flush=True)
        for g in sorted(curves,key=str):
            x,y = ps
            if index==1 and (g==3*x*x-3*x+1 or g==x*x-x+1):
                # Both equations cut out the same quadratic-field pencil.
                a,c,d = params
                curve = dict(zip(params,[cancel(sub.get(z,z)) for z in params]))
                for e in (3*c*c-3*c+1,a-(2*c-1)*d+c-1):
                    num = cancel(e.subs(curve,simultaneous=True)).as_numer_denom()[0]
                    assert Poly(num,x,y).rem(Poly(g,x,y)).is_zero
                skipped += 1
                continue
            nfibres += entire_pencil(co,local,curve_chart(g,ps))
            ncurves += 1
    assert (npoints,ncurves,nfibres,skipped)==(
        (26,10,32,0),(8,2,16,2),(12,10,26,0),(55,13,44,0))[index]
    if skipped:
        quadratic_line(coeff,params,units)
    print('PASS all common surfaces',index,len(common),npoints,ncurves,nfibres,flush=True)


def primary_intersection(params, units, charts, initial=()):
    """Walk through intersections of unions of hypersurfaces, retaining every dimension."""
    def polys(seq):
        return ','.join(str(v).replace('**','^') for v in sorted(seq,key=str))
    assert all(charts)
    source = ['LIB "primdec.lib";',f'ring r=0,({",".join(map(str,params))}),dp;',
              'short=0;option(redSB);',f'ideal U={polys(units)};']
    for i,fs in enumerate(charts):
        source.append(f'ideal C{i}={polys(fs)};')
    source.append('list cuts='+','.join('C'+str(i) for i in range(len(charts)))+';')
    source.append("""
proc descend(ideal P,int level){
 ideal H=std(P); int i,j; list L; ideal F,J;
 if(size(H)==1 && H[1]==1){return();}
 for(i=1;i<=size(U);i++){if(reduce(U[i],H)==0){return();}}
 if(level>size(cuts)){print("COMP "+string(dim(H))+" "+string(H));return();}
 F=cuts[level];
 for(i=1;i<=size(F);i++){if(reduce(F[i],H)==0){descend(H,level+1);return();}}
 if(dim(H)==0){return();}
 for(i=1;i<=size(F);i++){
   J=std(H+F[i]);
   if(size(J)!=1 || J[1]!=1){
     L=minAssGTZ(J);
     for(j=1;j<=size(L);j++){descend(L[j],level+1);}
   }
 }
}
""")
    source.append('descend(ideal('+(polys(initial) or '0')+'),1);print("DONE");quit;')
    out = subprocess.run([os.environ.get('SINGULAR','Singular'),'-q'],input='\n'.join(source),
                         capture_output=True,text=True,check=True).stdout
    assert 'DONE' in out.splitlines() and '? ' not in out,out[-1500:]
    result = {0:[],1:[]}
    for line in sorted(set(z for z in out.splitlines() if z.startswith('COMP '))):
        _,dim,row = line.split(' ',2)
        result[int(dim)].append(tuple(sympify(f.replace('^','**')) for f in row.split(',')))
    return result


def second_crossings(coeff, params, units, charts, common):
    fourth = ((2,1,0),(2,0,1),(0,1,2))
    cuts = charts+[parametric_check(permute(coeff,fourth),params)]
    components = primary_intersection(params,units|common,[fs-units-common for fs in cuts])
    print("Intersection counts",len(components[0]),len(components[1]),flush=True)
    assert (len(components[0]),len(components[1]))==(94,3)
    for basis in components[0]:
        print('PASS crossing point',basis,open_witness(coeff,params,basis,units),flush=True)
    a,c,d = params
    A,B = 3*a*a-3*a+1,-8*a*a+8*a-2
    v = 1/(1-t*t)
    changes = [
        ({a:1/t,c:t,d:(1+t-t*t)/t**3},c,1),
        ({a:t,c:t*t/(t*t+t-1),d:t**3+t*t-t},a,a*a+a-1),
        ({a:v,c:cancel(v*v/(2*v-1)),d:cancel(((-B+2*a*t*(2*a-1))/(2*A)).subs(a,v))},
         (2*A*d+B)/(2*a*(2*a-1)),(2*a-1)*A)]
    counts = []
    for basis,(sub,inverse,excluded) in zip(components[1],changes):
        assert all(cancel(h.subs(sub,simultaneous=True))==0 for h in basis)
        gb = groebner(basis,*params,order='grevlex',domain=QQ)
        poles = factors(cancel(inverse).as_numer_denom()[1])
        for z,val in sub.items():
            num = cancel(val.subs(t,inverse)-z).as_numer_denom()[0]
            assert gb.reduce(num)[1]==0
            denom = cancel(val).as_numer_denom()[1]
            poles.update(factors(cancel(denom.subs(t,inverse))))
        local = units|factors(excluded)
        for f in poles-local:
            bad = primary_intersection(params,local,[{f}],initial=basis)
            assert not any(bad.values()),(basis,f,bad)
        if excluded!=1:
            pts = primary_intersection(params,units,[{excluded}],initial=basis)
            assert not pts[1]
            for pt in pts[0]:
                print('PASS inverse-chart exceptional point',pt,open_witness(coeff,params,pt,units))
        counts.append(entire_pencil(coeff,units,sub))
    assert counts==[6,7,9]
    print('PASS entire second family away from common surfaces',counts,flush=True)


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--family',type=int,choices=range(4))
    parser.add_argument('--part',choices=('surfaces','crossings','all'),default='all')
    args = parser.parse_args()
    for index in (range(4) if args.family is None else (args.family,)):
        if args.part=='crossings' and index!=1:continue
        co,ps,u = family(index)
        charts = [parametric_check(permute(co,LABELS[label]),ps)-u for label in 'ACB']
        common = set.intersection(*charts)
        assert len(common)==(5,4,5,5)[index]
        if args.part!='crossings':surfaces(index,co,ps,u,common)
        if args.part!='surfaces' and index==1:second_crossings(co,ps,u,charts,common)
