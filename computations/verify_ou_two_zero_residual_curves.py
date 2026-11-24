#!/usr/bin/env python3
"""Test the generic point of every remaining two-zero curve over Q.

Requires SymPy and Singular with primdec.lib.  Finite-field points certify
nonempty characteristic-zero opens via a monic basis over Z localized at p.
"""
from argparse import ArgumentParser
from itertools import product
import os
import subprocess

from sympy import GF, Poly, cancel, groebner, lambdify, primerange, prod, symbols, sympify

from verify_ou_small_zero_corners import factors, verify_point
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check, verify_point as closed_fibre
from verify_ou_three_zero_remaining import LABELS, central_generation, t
from verify_ou_two_zero_codim_five import family
from verify_ou_two_zero_double_families import curve_chart, permute


def component_bases(params, charts, units):
    """List the codimension-two components of the triple intersection off common surfaces."""
    def entries(seq):
        return ','.join(str(f).replace('**','^') for f in seq)

    source = ['LIB "primdec.lib";',f'ring r=0,({entries(params)}),dp;',
              'short=0; option(redSB);',f'ideal U={entries(sorted(units,key=str))};']
    for i,fs in enumerate(charts):
        source.append(f'ideal C{i}={entries(sorted(fs,key=str))};')
    source.append("""
proc keep(ideal P,int needsThird)
{
  ideal H=std(P);
  if(size(H)==1 && H[1]==1){return();}
  if(dim(H)!=nvars(basering)-2){ERROR("unexpected dimension");}
  int w;
  for(w=1;w<=size(U);w++){if(reduce(U[w],H)==0){return();}}
  if(needsThird){
    for(w=1;w<=size(C2);w++){if(reduce(C2[w],H)==0){break;}}
    if(w>size(C2)){return();}
  }
  print("COMP "+string(H));
}
list L; int i,j,k,v;
for(i=1;i<=size(C0);i++){for(j=1;j<=size(C1);j++){
  if(C0[i]==C1[j]){
    for(k=1;k<=size(C2);k++){
      L=minAssGTZ(ideal(C0[i],C2[k]));
      for(v=1;v<=size(L);v++){keep(L[v],0);}
    }
  }else{
    L=minAssGTZ(ideal(C0[i],C1[j]));
    for(v=1;v<=size(L);v++){keep(L[v],1);}
  }
}}
print("DONE");quit;
""")
    out = subprocess.run([os.environ.get('SINGULAR','Singular'),'-q'],
                         input='\n'.join(source),text=True,capture_output=True,check=True).stdout
    assert 'DONE' in out.splitlines() and '? ' not in out,out[-1200:]
    bases = sorted({line[5:] for line in out.splitlines() if line.startswith('COMP ')})
    return [tuple(sympify(f.replace('^','**')) for f in row.split(',')) for row in bases]


def open_witness(coeff, params, basis, units):
    # Singular supplies a dp basis over Q.  Invert the denominators and leading
    # coefficients so that polynomial division by its monic normalization works
    # integrally.  Standard monomials then give a flat Z_(p)-model.
    bad = []
    for f in basis:
        poly = Poly(f,*params)
        bad.extend([poly.LC(order='grevlex').p,*(c.q for c in poly.coeffs())])
    equations = [lambdify(params,Poly(f,*params).clear_denoms()[1].as_expr(),modules='math')
                 for f in basis]
    nonzero = [lambdify(params,f,modules='math') for f in units]
    for prime in primerange(5,100):
        if any(n%prime==0 for n in bad):
            continue
        field = GF(prime)
        for point in product(range(prime),repeat=len(params)):
            if any(f(*point)%prime for f in equations) or any(f(*point)%prime==0 for f in nonzero):
                continue
            sub = dict(zip(params,point))
            values = {w:v.subs(sub) for w,v in coeff.items()}
            if any(not v.is_Rational or v.q%prime==0 for v in values.values()):
                continue
            try:
                verify_point({w:field(int(v.p))/field(int(v.q)) for w,v in values.items()},field)
            except (AssertionError,ZeroDivisionError):
                continue
            return prime,point
    raise AssertionError(('no flat regular-corner witness',basis))


def surface_chart(f, params, units):
    # A linear equation up+v=0 is a graph unless u=v=0.  The last
    # test checks that vertical locus in the localized parameter plane.
    choices = []
    for p in params:
        poly = Poly(f,p)
        if poly.degree()==1:
            leading,constant = poly.LC(),poly.TC()
            value = cancel(-constant/leading)
            vertical = bool(factors(leading)-units and factors(constant)-units)
            choices.append((vertical,len(str(value)),str(p),p,leading,constant,value))
    for vertical,_,_,p,leading,constant,value in sorted(choices):
        others = tuple(z for z in params if z!=p)
        if vertical:
            inv = symbols('inv')
            keep = prod(h for h in units if not h.has(p))
            if tuple(groebner([leading,constant,inv*keep-1],inv,*others,order='grevlex'))!=(1,):
                continue
        return {p:value},others
    a,c,d = params
    assert (f-(a-c)**2+c*(c-1)*(d-1)**2).expand()==0
    assert {c,c-1,d-1,a-c}<=units
    w = symbols('w')
    # Inverse w=(a-c)/(c*(d-1)), so 1-w^2=1/c is a unit.
    return {a:(1+(d-1)*w)/(1-w*w),c:1/(1-w*w)},(w,d)


def verify(index):
    coeff,params,units = family(index)
    charts = [parametric_check(permute(coeff,LABELS[label]),params)-units for label in 'ACB']
    common = set.intersection(*charts)
    rest = [fs-common for fs in charts]
    bases = component_bases(params,rest,units|common)
    assert len(bases)==(21,9,6,46)[index] and len(common)==(5,4,5,5)[index]
    for basis in bases:
        print('PASS crossing curve',index,basis,open_witness(coeff,params,basis,units),flush=True)
    counts = []
    for h in sorted(common,key=str):
        sub,ps = surface_chart(h,params,units)
        assert cancel(h.subs(sub,simultaneous=True))==0
        c = apply_sub(coeff,sub)
        values = [cancel(f.subs(sub,simultaneous=True)) for f in units]
        assert all(values)
        local = set().union(*(factors(g) for g in values))
        cuts = parametric_check(permute(c,LABELS['A']),ps)-local
        for f in sorted(cuts,key=str):
            print('PASS surface curve',index,h,f,open_witness(c,ps,(f,),local),flush=True)
        counts.append(len(cuts))
    assert sum(counts)==(47,28,38,47)[index]
    print('PASS two-zero curve components',index,len(bases),counts,flush=True)


def verify_minus_one(index):
    coeff,params,units = family(index)
    d = params[1] if index in (0,2) else params[2]
    assert str(d)=='v3'
    coeff = apply_sub(coeff,{d:-1})
    ps = tuple(v for v in params if v!=d)
    values = [cancel(v.subs(d,-1)) for v in units]
    assert all(values)
    local = set().union(*(factors(v) for v in values))
    charts = [parametric_check(permute(coeff,LABELS[label]),ps)-local for label in 'ACB']
    common = set.intersection(*charts)
    bases = component_bases(ps,[fs-common for fs in charts],local|common)
    assert len(bases)==(20,9,2,5)[index]
    for basis in bases:
        print('PASS isolated intersection',index,basis,open_witness(coeff,ps,basis,local),flush=True)
    counts = []
    for f in sorted(common,key=str):
        sub = curve_chart(f,ps)
        co = apply_sub(coeff,sub)
        values = [cancel(v.subs(sub,simultaneous=True)) for v in local]
        assert all(values)
        try:
            failed = parametric_check(co,(t,))
        except AssertionError:
            failed = central_generation(co)
        cuts = failed-set().union(*(factors(v) for v in values))
        for h in sorted(cuts,key=str):
            print('PASS closed fibre',index,f,h,closed_fibre(co,values,h),flush=True)
        counts.append(len(cuts))
    assert (len(common),sum(counts))==((4,38),(4,38),(4,33),(2,4))[index]
    print('PASS entire d=-1 slice',index,counts,flush=True)


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--family',type=int,choices=range(4))
    parser.add_argument('--part',choices=('curves','minus-one','all'),default='all')
    args = parser.parse_args()
    for index in (range(4) if args.family is None else (args.family,)):
        if args.part!='minus-one':
            verify(index)
        if args.part!='curves':
            verify_minus_one(index)
