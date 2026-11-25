#!/usr/bin/env python3
"""Finish the first and fourth two-zero families by exact cuts on each crossing curve."""
from argparse import ArgumentParser
from sympy import QQ, groebner

from verify_ou_two_zero_closed_strata import (
    LABELS, family, local_cuts, open_witness, primary_intersection,
)
from verify_ou_two_zero_third_family import five_chart_intersection


# All computations take place over Q. These small residues specify the order of pivots.
SAMPLES = {
    0: [(LABELS['A'],7,(2,6,3)),(LABELS['C'],7,(2,2,3)),
        (LABELS['B'],5,(3,2,3)),(((2,1,0),(2,0,1),(0,1,2)),7,(4,2,6))],
    3: [(LABELS['A'],7,(2,4,4)),(LABELS['B'],7,(4,5,4)),
        (((1,0,2),(1,2,0),(2,1,0)),5,(3,4,3))],
}


def verify(index):
    coeff,params,units = family(index)
    components = five_chart_intersection(coeff,params,units)
    assert (len(components[0]),len(components[1]))=={0:(114,5),3:(163,19)}[index]
    for basis in components[0]:
        print('PASS crossing point',index,basis,open_witness(coeff,params,basis,units),flush=True)
    remaining = {basis:groebner(basis,*params,order='grevlex',domain=QQ) for basis in components[1]}
    newly_covered = []
    points = 0
    for label,prime,point in SAMPLES[index]:
        cuts = local_cuts(coeff,params,label,(prime,point))-units
        covered = [basis for basis,gb in remaining.items() if all(gb.reduce(f)[1]!=0 for f in cuts)]
        for basis in covered:
            local = primary_intersection(params,units,[cuts],initial=basis) if cuts else {0:[],1:[]}
            assert not local[1]
            for pt in local[0]:
                print('PASS point on crossing curve',index,pt,open_witness(coeff,params,pt,units),flush=True)
            points += len(local[0])
            del remaining[basis]
        newly_covered.append(len(covered))
    assert not remaining
    assert (newly_covered,points)=={0:([1,2,1,1],21),3:([11,5,3],56)}[index]
    print('PASS entire two-zero family',index,newly_covered,points,flush=True)


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--family',type=int,choices=(0,3))
    args = parser.parse_args()
    for index in ((0,3) if args.family is None else (args.family,)):
        verify(index)
