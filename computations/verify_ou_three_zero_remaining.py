#!/usr/bin/env python3
"""All curve and closed-point specializations on the twelve remaining three-zero surfaces."""
from argparse import ArgumentParser

from sympy import Poly, QQ, cancel, factor_list

from verify_ou_complement_cycles import face_dets
from verify_ou_corner_comparison import t
from verify_ou_five_zero_supports import expected
from verify_ou_seven_zero_center import CyclicPaths, row_matrix
from verify_ou_small_zero_corners import factors, graph
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check, verify_point
from verify_ou_small_zero_surface_orbits import chart, representatives
from verify_ou_zero_supports import I, TRIPLES


# Coordinate choices only affect the pivot divisors that must be examined.
LABELS = {
    'A': ((1,0,2),(0,1,2),(2,0,1)),
    'B': ((1,0,2),(0,1,2),(2,1,0)),
    'C': ((1,2,0),(2,0,1),(1,2,0)),
}
PATTERNS = {134151919:'BCAAABAA',134217454:'CAAA'}


def curve_parameters(f, ps):
    x,y = ps
    charts = ({x:x,y:y},{x:x/y,y:y},{x:x,y:y/x},
              {x:x,y:x*y},{x:x*y,y:y})
    for change in charts:
        h = cancel(f.subs(change,simultaneous=True)).as_numer_denom()[0]
        # Remove only monomials invertible on the torus.
        divisors = [g for g,_ in factor_list(h)[1] if g not in ps]
        if not divisors:
            # A coordinate line x=0 or y=0 is already outside the torus.
            raise AssertionError(('coordinate divisor',f))
        assert len(divisors)==1
        for pivot,param in ((x,y),(y,x)):
            poly = Poly(divisors[0],pivot)
            if poly.degree()!=1:
                continue
            assert Poly(poly.LC(),param).gcd(Poly(poly.TC(),param)).degree()==0
            val = cancel(-poly.TC()/poly.LC())
            sub = {v:cancel(w.subs(pivot,val).subs(param,t)) for v,w in change.items()}
            assert cancel(f.subs(sub,simultaneous=True))==0
            return sub
    # This affine change is invertible, with parameter x-y at each point.
    poly = Poly(f.subs(y,x-y),x)
    assert poly.degree()==1
    assert Poly(poly.LC(),y).gcd(Poly(poly.TC(),y)).degree()==0
    val = cancel((-poly.TC()/poly.LC()).subs(y,t))
    sub = {x:val,y:val-t}
    assert cancel(f.subs(sub,simultaneous=True))==0
    return sub


def central_generation(c):
    field = QQ.frac_field(t)
    paths = CyclicPaths(0,field,{w:field.from_sympy(v) for w,v in c.items()})
    centers = paths.center(3)
    assert len(centers)>=2
    cuts = set()

    def record_poles(v):
        if not v.denom.is_ground:
            cuts.update(f.as_expr() for f,_ in v.denom.factor_list()[1])

    for z in centers:
        for v in z.values():
            record_poles(v)
    for n in range(2,6):
        for r,i,j in TRIPLES:
            basis,reductions = paths.space(n,r,i,j)
            assert len(basis)==expected(n,i,j)
            for row in reductions.values():
                for v in row.values():
                    record_poles(v)

    for r,i,j in TRIPLES:
        width = len(paths.space(5,r,i,j)[0])
        rows = []
        for z in centers:
            for tail in paths.space(2,r,i,j)[0]:
                row = {}
                for (rotation,loop),coefficient in z.items():
                    if rotation==r and loop[0]==i:
                        for k,v in paths.word(r,loop+tail[1:]).items():
                            row[k] = row.get(k,field.zero)+coefficient*v
                rows.append(row)
        matrix = row_matrix(rows,width,field)
        basis = matrix.transpose().rref()[1]
        assert len(basis)==width
        determinant = matrix.extract(basis,list(range(width))).det()
        cuts.update(factors(field.to_sympy(determinant)))
    for fn in (paths.space,paths.word):
        fn.cache_clear()
    return cuts


def verify_surface(c, ps, label):
    perms = LABELS[label]
    permuted = {w:c[tuple(p[i] for p,i in zip(perms,w))] for w in TRIPLES}
    units = set().union(*(factors(v) for v in c.values() if v))
    units.update(set().union(*(factors(v) for v in face_dets(c))))
    cuts = sorted(parametric_check(permuted,ps)-units,key=str)
    points = 0
    central = 0
    for g in cuts:
        sub = curve_parameters(g,ps)
        cc = apply_sub(c,sub)
        local = [cancel(v.subs(sub,simultaneous=True)) for v in units]
        assert all(local)
        try:
            remaining = parametric_check(cc,(t,))
        except AssertionError:
            remaining = central_generation(cc)
            central += 1
        excluded = set().union(*(factors(v) for v in local))
        for h in sorted(remaining-excluded,key=str):
            result = verify_point(cc,local,h)
            print('PASS closed fibre',g,h,result,flush=True)
            points += 1
        print('PASS entire curve',g,flush=True)
    return len(cuts),central,points


def verify(mask, selected=None):
    counts = []
    rows = representatives(mask)
    assert len(rows)==len(PATTERNS[mask])
    for index,((zero,f,g),label) in enumerate(zip(rows,PATTERNS[mask])):
        if selected is not None and selected!=index:
            continue
        _,coefficients,vs,_ = chart(mask,zero)
        pivot,val = graph(f,vs,require_unit=True)
        params = [v for v in vs if v!=pivot]
        c = apply_sub(coefficients,{pivot:val})
        pivot,val = graph(g,params,require_unit=True)
        c = apply_sub(c,{pivot:val})
        stats = verify_surface(c,tuple(v for v in params if v!=pivot),label)
        counts.append(stats)
        print('PASS entire surface',mask,index,stats,flush=True)
    assert counts
    print('PASS requested surfaces',mask,counts,flush=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--mask',type=int,choices=tuple(PATTERNS))
    parser.add_argument('--surface',type=int)
    args = parser.parse_args()
    for m in ((args.mask,) if args.mask else PATTERNS):
        verify(m,args.surface)
