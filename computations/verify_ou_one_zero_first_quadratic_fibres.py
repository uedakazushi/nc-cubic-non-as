#!/usr/bin/env python3
"""Conjugate exceptional fibres in the first one-zero family."""
from verify_ou_one_zero_first_two_slices import normalized_parent
from verify_ou_one_zero_closed_divisors import LABEL, points
from verify_ou_three_zero_double_pencils import parametric_check
from verify_ou_two_zero_closed_strata import local_cuts, primary_intersection
from verify_ou_two_zero_double_families import permute


SAMPLES = (
    (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (7, (2, 3, 6))),
    (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (7, (3, 4, 6))),
    (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (7, (2, 2, 4))),
    (((1, 2, 0), (2, 0, 1), (1, 2, 0)), (7, (2, 6, 6))),
    (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (7, (6, 2, 6))),
    (((1, 0, 2), (0, 1, 2), (2, 0, 1)), (7, (2, 5, 6))),
)
COUNTS = [9, 12, 14, 1, 10, 8]


def intersections(params, guards, cuts, basis):
    if not cuts:
        return {0: [], 1: []}
    # Singular's list-of-ideals interface requires at least two chart ideals.
    return primary_intersection(params, guards, [cuts, cuts], initial=basis)


def verify():
    params, coeff, guards, cuts = normalized_parent('linear')
    a, d, e = params
    special = e**2-3*e+3
    assert cuts[-1] == special
    regular = parametric_check(permute(coeff, LABEL['B']), params)-guards-{-v for v in guards}
    base = intersections(params, guards, regular, [special])
    assert not base[0]
    expected = [3*a*d+e-3, 3*a*a*d-6*a*d-2*a*e+3*a+2*e,
                3*a*a*d*d+a*a*d*e-3*a*a*d-a*d*e-a*d-a*e+2*a+e-1,
                a*d-e+1, a*a*d+2*a*d*e-3*a*d-a*e-e+2,
                a*a*d-a*d*e+e-1]
    assert base[1] == [(special, f) for f in expected]
    counts = []
    for i, basis in enumerate(base[1]):
        label, sample = SAMPLES[i]
        extra = local_cuts(coeff, params, label, sample)-guards-{-v for v in guards}
        residual = intersections(params, guards, extra, basis)
        assert not residual[1]
        counts.append(points(coeff, params, guards, residual[0]))
        print('PASS conjugate curve', i, counts[-1], flush=True)
    assert counts == COUNTS
    print('PASS both conjugate fibres and entire first-family representative', counts, flush=True)


if __name__ == '__main__':
    verify()
