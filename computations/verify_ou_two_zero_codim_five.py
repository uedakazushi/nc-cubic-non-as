#!/usr/bin/env python3
"""A certified regular-corner point on every remaining two-zero surface divisor."""
from argparse import ArgumentParser

from sympy import QQ, cancel, sqrt, sympify

from verify_ou_complement_cycles import face_dets
from verify_ou_small_zero_corners import factors, verify_point
from verify_ou_small_zero_surface_orbits import chart
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import LABELS
from verify_ou_two_zero_double_families import permute


# Coordinates on each irreducible divisor, in the order of its polynomial's string.
WITNESSES = (
    (
        '3 -1 -2', '-3 -2 4', '-3 2 6', '4 6 -2',
        '3 -4 -1', '-1/2 2 1', '2 4 5', '2 2 -2/3',
        '2 4 -1', '-2/3 2/3 1', '4/3 4 2/3', '5 -2/3 2',
        '1/2 4 5', '2 2 -1', '5 2/3 -1', '3 1/2 -1/2',
        '5 -1 4/3', '2 6 5',
    ),
    (
        '-3 4 -1', '1/2 3 1/2', '-2 -1/2 -2', '2/3 4/3 2',
        '4 -1 6', '-1 -2/3 2', '1/2 -2 -2/3', '-3 -1/2 -2/3',
        '2/3 1/2 2', '-1/2 -1 3/2', '4 2 -3', '3 4/3 3/2',
    ),
    (
        '2/3 -2/3 -1/2', '-2 2 4', '-1 3/2 4', '-1 1/2 4/3',
        '-1/2 4 -3', '2 2/3 -1/2', '-1/2 4/3 4/3', '-1/2 -1 1/2',
        '4 -3 3', '-2 -2 -4', '-1 -3 -2', '-3 -2 -2',
        '-4 -1 3',
    ),
    (
        '4 -4 -1', '2/3 1/2 2/3', '-1 2 -3', '3/2 2/3 3/2',
        '-2 -1/2 2', '2 3 2', '-1 4 -4', '-1 3 -3',
        '-2 -3 2', '2-sqrt(2) 2 -1', '2/3 1/2 2', '-1 -3 2/3',
        '-4 -3 4', '-1/2 1/2 -1', '-2/3 -4 3', '-2 -2/3 4/3',
        '6 -4 -1',
    ),
)


def family(index):
    _,template,(a,b,c,d,e),_ = chart(134217455,(0,1,1))
    forms = (
        ((c,d,e),{a:d-e,b:1-c}),
        ((a,c,d),{b:(a-a*c+c**2*d-c*d)/(a-c),e:(-a**2+a*c*d+a-c*d)/(a-c)}),
        ((c,d,e),{a:1-e,b:1-c}),
        ((a,c,d),{b:1-c,e:d*(1-a)}),
    )
    params,sub = forms[index]
    co = apply_sub(template,sub)
    units = set().union(*(factors(v) for v in co.values() if v))
    units.update(set().union(*(factors(v) for v in face_dets(co))))
    if index==1:
        # The original exterior face with determinant c-a excludes the vertical fibre.
        assert c-a in face_dets(template)
        _,den = cancel(sub[e]).as_numer_denom()
        assert factors(den)<=units
    return co, params, units | {d-1}


def verify(index):
    co,params,units = family(index)
    cuts = parametric_check(permute(co,LABELS['A']),params)-units
    assert len(cuts)==len(WITNESSES[index])
    for h,coordinates in zip(sorted(cuts,key=str),WITNESSES[index]):
        point = tuple(sympify(x) for x in coordinates.split())
        values = dict(zip(params,point))
        assert h.subs(values).expand()==0 and all(v.subs(values).expand()!=0 for v in units)
        field = QQ if all(x.is_Rational for x in point) else QQ.algebraic_field(sqrt(2))
        coeff = {w:field.from_sympy(v.subs(values)) for w,v in co.items()}
        verify_point(coeff,field)
        print('PASS divisor',index,h,point,flush=True)
    print('PASS all divisors',index,len(cuts),flush=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--family',type=int,choices=range(4))
    index = parser.parse_args().family
    for n in (range(4) if index is None else (index,)):
        verify(n)
