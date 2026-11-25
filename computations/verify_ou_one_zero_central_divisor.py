#!/usr/bin/env python3
"""Certify all specializations on the common central-rank divisor of both one-zero families."""
import json
from itertools import product
from pathlib import Path

from sympy import QQ, cancel, sympify

from verify_ou_complement_cycles import face_dets
from verify_ou_seven_zero_center import CyclicPaths, I
from verify_ou_small_zero_corners import factors
from verify_ou_small_zero_surface_orbits import chart
from verify_ou_three_zero_double_pencils import apply_sub
from verify_ou_three_zero_remaining import LABELS, generation_factors
from verify_ou_two_zero_double_families import permute
from verify_ou_two_zero_closed_strata import local_units


def certified_cuts(coeff, variables, sections):
    field = QQ.frac_field(*variables)
    labelled = permute(coeff,LABELS['A'])
    paths = CyclicPaths(0,field,{w:field.from_sympy(v) for w,v in labelled.items()})
    centers = []
    for section in sections:
        z = {(r,tuple(word)):field.from_sympy(sympify(v)) for r,word,v in section}
        assert len(z)==len(section)
        for r,word in z:
            assert r in I and len(word)==4 and word[0]==word[-1]
            assert word in paths.space(3,r,word[0],word[-1])[0]
        for r,i,j in product(I,repeat=3):
            commutator = {}
            for (turn,word),v in z.items():
                if turn==r and word[0]==i:
                    for k,x in paths.word(r,word+(j,)).items():
                        commutator[k] = commutator.get(k,field.zero)+v*x
                if turn==(r+1)%3 and word[0]==j:
                    for k,x in paths.word(r,(i,)+word).items():
                        commutator[k] = commutator.get(k,field.zero)-v*x
            assert not any(commutator.values())
        centers.append(z)
    return generation_factors(paths,centers)


def verify():
    file = Path(__file__).resolve().parents[1]/'certificates/ou_one_zero_central_sections.json'
    first,second = json.loads(file.read_text())
    _,template,(a,b,c,d,e,h),_ = chart(134217471,(0,2,2))
    coeff = apply_sub(template,{b:(c+d*h-h)/c,e:(a*c+a*d*h-a*h-c)/(a*c*d-c)})
    units = set().union(*(factors(v) for v in coeff.values() if v),
                        *(factors(v) for v in face_dets(coeff)))
    sub = {a:(c+d-1)/(c*d)}
    cs = apply_sub(coeff,sub)
    local = local_units(units,sub)
    q = (c*(c-1)*(d-1)+(c+d-1)*h).expand()
    assert certified_cuts(cs,(c,d,h),first)-local=={q}
    sub = {h:-c*(c-1)*(d-1)/(c+d-1)}
    assert c+d-1 in local and cancel(q.subs(sub))==0
    assert not certified_cuts(apply_sub(cs,sub),(c,d),second)-local_units(local,sub)
    print('PASS common one-zero central divisor including its entire last surface',flush=True)


if __name__=='__main__':
    verify()
