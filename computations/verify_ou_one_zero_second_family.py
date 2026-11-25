#!/usr/bin/env python3
"""Dimension at most two for the undecided locus in the second four-parameter one-zero family."""
from sympy import cancel

from verify_ou_complement_cycles import face_dets
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import LABELS, central_generation
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import open_witness


def verify():
    _,template,(a,b,c,d,e,h),_ = chart(134217471,(0,2,2))
    assert d-b in face_dets(template)
    sub = {b:(c+d*h-h)/c,e:(a*c+a*d*h-a*h-c)/(a*c*d-c)}
    equations = (a*b-a*d*e+e-1,a*c*d*e-a*c-a*d*h+a*h-c*e+c)
    assert all(cancel(eq.subs(sub,simultaneous=True))==0 for eq in equations)
    ps = a,c,d,h
    co = apply_sub(template,sub)
    units = set().union(*(factors(v) for v in co.values() if v),
                        *(factors(v) for v in face_dets(co)))
    labelled = permute(co,LABELS['A'])
    cuts = parametric_check(labelled,ps)-units
    special = a*c*d-c-d+1
    assert len(cuts)==18 and special in cuts
    for f in sorted(cuts-{special},key=str):
        print('PASS flat divisor witness',f,open_witness(co,ps,(f,),units),flush=True)
    specialized = apply_sub(labelled,{a:(c+d-1)/(c*d)})
    residual = central_generation(specialized,(c,d,h))
    print('PASS central generation on the exceptional divisor',len(residual),flush=True)
    print('PASS generic points of all eighteen one-zero divisors',flush=True)


if __name__=='__main__':
    verify()
