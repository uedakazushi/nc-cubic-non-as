#!/usr/bin/env python3
"""Close the entire second one-zero family after the two previously settled orbits."""
from sympy import Poly, cancel

from verify_ou_one_zero_closed_divisors import LABEL, points, three_charts
from verify_ou_one_zero_residual_surfaces import family
from verify_ou_one_zero_second_common import parent_orbit_data
from verify_ou_small_zero_corners import factors
from verify_ou_three_zero_double_pencils import apply_sub, parametric_check
from verify_ou_two_zero_closed_strata import local_units
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_residual_curves import component_bases, surface_chart


def verify():
    co,ps,units = family(1)
    settled,orbits = parent_orbit_data(co,ps,units)
    guards = units|settled
    a,c,d,h = ps
    cases = [(a*c-c+d-1,d,a*c*h-c-h,6),
             (a*c*d+a*d**2*h-a*d*h-c,h,a*c*d-a*c-a*d+1,8),
             (a*c**2+a*c*d-a*c-c**2-d+1,a,h-1,0)]
    for equation,pivot,remaining,count in cases:
        poly = Poly(equation,pivot)
        assert poly.degree()==1
        if pivot==a:
            # At d=1-c the constant coefficient is c(1-c), a unit.
            assert factors(poly.LC())-guards=={c+d-1}
            assert factors(poly.TC().subs(d,1-c))<=guards
        else:
            assert factors(poly.LC())<=guards
        sub = {pivot:cancel(-poly.TC()/poly.LC())}
        params = tuple(v for v in ps if v!=pivot)
        local = local_units(guards,sub)|factors(poly.LC())
        coeff = apply_sub(co,sub)
        cuts = parametric_check(permute(coeff,LABEL['B']),params)-local
        known = local_units(orbits[1]|orbits[2],sub)
        assert len(cuts)==4 and cuts-known=={remaining}
        next_sub,variables = surface_chart(remaining,params,local)
        second = apply_sub(coeff,next_sub)
        if pivot==a:
            # The last graph needs only the B test on the complement of the
            # two prior orbits. Here their pullbacks are nowhere identically zero.
            local = local_units(local|known,next_sub)
            assert not (parametric_check(permute(second,LABEL['B']),variables)-local)
        else:
            local = local_units(local,next_sub)
            charts = three_charts(second,variables,local)
            assert not set.intersection(*charts)
            pts = component_bases(variables,charts,local)
            assert points(second,variables,local,pts)==count
        print('PASS entire representative',equation,flush=True)


if __name__=='__main__':
    verify()
