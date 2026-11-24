#!/usr/bin/env python3
"""Exact smooth plane sections on the seven remaining one-to-three-zero supports."""
from sympy import Matrix, QQ, groebner, symbols

from verify_ou_four_zero_supports import laurent_rank, rows_for
from verify_ou_seven_zero_curves import determinant_matrix
from verify_ou_zero_supports import TRIPLES, incidence, normalize


SMOOTH_MASKS = {
    1: (134217471,),
    2: (134209279, 134217455),
    3: (134151919, 134184687, 134216943, 134217454),
}


def check_cut(mask, zero_count):
    c, field, free = normalize(mask)
    assert 27 - mask.bit_count() == zero_count
    assert len(free) == 8 - zero_count
    values = dict(zip(field.symbols, (2, 3, 5, 7, 11, 13, 17)))
    sample = {t: field.to_sympy(v).subs(values) for t, v in c.items()}
    matrix = determinant_matrix(sample)
    assert matrix.rank() == 3
    coords = symbols('x y z')
    plane = Matrix.hstack(*matrix.nullspace()) * Matrix(coords)
    cubic = plane[0]*plane[3]*plane[4] - plane[1]*plane[2]*plane[5]
    gradient = [cubic.diff(v) for v in coords]
    for v in coords:
        other = [w for w in coords if w != v]
        assert list(groebner([g.subs(v, 1) for g in gradient], *other, domain=QQ)) == [1]
    print(f'PASS {zero_count} zeros, mask {mask}, three unit gradient ideals')


def check_singular_exterior_face():
    core = [t for t in TRIPLES if t[1] != 0 and t[2] != 0]
    ordered = core + [t for t in TRIPLES if t != (0, 0, 0) and t not in core]
    pivots = incidence(ordered).T.rref()[1]
    free = [t for i, t in enumerate(ordered) if t not in core and i not in pivots]
    assert free == [(1, 0, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 0)]
    field = QQ.frac_field(*symbols('u:5'))
    coeff = {t: field.one if t in ordered else field.zero for t in TRIPLES}
    coeff.update(zip(free, field.gens))
    size, rows = rows_for(coeff, 6, 0, 0, 0, field)
    assert (len(rows), size) == (405, 243)
    assert laurent_rank(rows, field) == 231
    print('PASS singular exterior face, length 6, generic quotient 12 > 10')


if __name__ == '__main__':
    for zeros, masks in SMOOTH_MASKS.items():
        for support_mask in masks:
            check_cut(support_mask, zeros)
    check_singular_exterior_face()
