#!/usr/bin/env python3
"""The first one-zero family has residual dimension at most two."""
from sympy import cancel

from verify_ou_complement_cycles import face_dets
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_three_zero_remaining import LABELS
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import open_witness


def verify():
    _,template,(a,b,c,d,e,h),_ = chart(134217471,(0,2,2))
    assert d-b in face_dets(template)
    sub = {b:(a*d*e-e+1)/a,h:(a*d*e-a-e+1)/(a*(a*d-1))}
    eqs = (a*b-a*d*e+e-1,a*a*d*h-a*d*e-a*h+a+e-1)
    assert all(cancel(eq.subs(sub,simultaneous=True))==0 for eq in eqs)
    ps = a,c,d,e
    co = apply_sub(template,sub)
    units = set().union(*(factors(v) for v in co.values() if v),
                        *(factors(v) for v in face_dets(co)))
    cuts = parametric_check(permute(co,LABELS['C']),ps)-units
    special = a*c*d-c-d+1
    assert len(cuts)==12 and special in cuts
    for f in sorted(cuts-{special},key=str):
        print('PASS flat divisor witness',f,open_witness(co,ps,(f,),units),flush=True)
    a0 = (c+d-1)/(c*d)
    h0 = c*(c+d*e-1)/(c+d-1)
    e0 = (a*c+a*d*h-a*h-c)/(a*c*d-c)
    assert cancel(sub[h].subs(a,a0)-h0)==0
    assert cancel(e0.subs({a:a0,h:h0},simultaneous=True)-e)==0
    assert cancel(((c+d*h-h)/c).subs(h,h0)-sub[b].subs(a,a0))==0
    assert cancel(h0.diff(e)-c*d/(c+d-1))==0
    print('PASS invertible parameter change to the second family on the special divisor',flush=True)


if __name__=='__main__':
    verify()
