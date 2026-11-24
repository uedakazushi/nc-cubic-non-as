#!/usr/bin/env python3
"""Line and conic on the generic complementary binomial of every few-zero support."""
from itertools import product
from random import Random

from sympy import Matrix, expand, factor_list, symbols

from verify_ou_few_zero_cuts import SMOOTH_MASKS
from verify_ou_zero_supports import I, TRIPLES, incidence, support


x, y = symbols("x y")


def pencil_matrix(c):
    l, p, q, r, s = variables = symbols("l p q r s")
    equations = [
        -c[0,1,0]*(c[1,0,0]*c[2,2,0]*s-c[2,0,0]*c[1,2,0]*q)
        +c[0,2,0]*(c[1,0,0]*c[2,1,0]*r-c[2,0,0]*c[1,1,0]*p)
    ]
    for k in (1, 2):
        equations.append(c[0,0,k]*l + c[0,1,k]*(c[2,0,k]*q-c[1,0,k]*s)
                         + c[0,2,k]*(c[1,0,k]*r-c[2,0,k]*p))
    return Matrix(equations).jacobian(variables)


def pencil(c):
    basis = pencil_matrix(c).nullspace()
    assert len(basis) == 2
    return Matrix.hstack(*basis)*Matrix([x, y])


def independent(f, g):
    return f.diff(x)*g.diff(y) != f.diff(y)*g.diff(x)


def verify(c):
    l, p, q, r, s = pencil(c)
    delta = p*s-q*r
    for f in (l, p, q, r, s):
        assert f != 0 and delta.subs({x: f.diff(y), y: -f.diff(x)}) != 0
    assert all(independent(f, g) for f, g in ((p, r), (p, q), (r, s)))
    for axis in I:
        face = Matrix(2, 2, lambda i,j: c[(0,i+1,j+1) if axis == 0 else
                                          (i+1,0,j+1) if axis == 1 else (i+1,j+1,0)])
        assert face.det() != 0
    conic = Matrix([[l,delta,delta], [1,p,q], [1,r,s]])
    line = Matrix([[1,0,0], [r,1,1], [p,1,1]])
    for k in I:
        ck = Matrix(3, 3, lambda i,j: c[i,j,k])
        if k == 0:
            u, v = c[0,2,0], -c[0,1,0]
            vn = Matrix([-(u*c[1,1,0]*p+v*c[1,2,0]*q)/c[1,0,0], u, v])
            tail = -ck[1:,1:].inv()*Matrix([c[1,0,0]*r, c[2,0,0]*p])
            vl = Matrix([1, *tail])
        else:
            u = c[2,0,k]*q-c[1,0,k]*s
            v = c[1,0,k]*r-c[2,0,k]*p
            assert independent(u, v)
            vn, vl = Matrix([delta,u,v]), Matrix([0,1,-1])
        for matrix, vector in ((conic,vn), (line,vl)):
            assert all(expand(f) == 0 for f in ck.multiply_elementwise(matrix)*vector)


def boundary_factors(c, variables):
    # A bracket is the determinant of the pencil matrix bordered by two forms.
    # For coordinate forms these are signed 3-by-3 complementary minors.
    mat = pencil_matrix(c)
    brackets = Matrix.zeros(5)
    for i in range(5):
        for j in range(i+1,5):
            cols = [k for k in range(5) if k not in (i,j)]
            brackets[i,j] = (-1)**(i+j+1)*mat[:,cols].det()
            brackets[j,i] = -brackets[i,j]
    polys = [brackets[i,j] for i,j in ((1,2),(1,3),(2,4),(3,4))]
    polys.append(brackets[0,1]*brackets[0,4]-brackets[0,2]*brackets[0,3])
    for k in (1,2):
        b = Matrix([0,0,c[2,0,k],0,-c[1,0,k]])
        e = Matrix([0,-c[2,0,k],0,c[1,0,k],0])
        polys.append((b.T*brackets*e)[0])

    def factors(entries):
        out = set()
        for f in entries:
            assert expand(f) != 0
            out.update(p for p, _ in factor_list(f)[1] if p not in variables)
        return out

    faces = [Matrix(2, 2, lambda i,j: c[(0,i+1,j+1) if axis == 0 else
                                         (i+1,0,j+1) if axis == 1 else (i+1,j+1,0)]).det()
             for axis in I]
    return factors(polys) - factors(faces)


def normalize_binomial(sup):
    cube = sorted(t for t in sup if all(t))
    axes = sorted(t for t in sup if t.count(0) == 2)
    rows = cube + axes + sorted(sup - set(cube) - set(axes))
    pivots = incidence(rows).T.rref()[1]
    assert len(pivots) == 19 and len(cube) == 8
    free = [t for i,t in enumerate(rows) if i not in pivots and t not in cube]
    assert len(free) == len(sup)-20
    variables = symbols('v:' + str(len(free)))
    return {t: variables[free.index(t)] if t in free else int(t in sup) for t in TRIPLES}, variables


def main():
    count = 0
    for masks in SMOOTH_MASKS.values():
        for mask in masks:
            zeros = set(TRIPLES) - set(support(mask))
            for zero in sorted(zeros):
                if any(all(i != j for i, j in zip(zero, t)) for t in zeros):
                    continue
                labels = [[i]+[j for j in I if i != j] for i in zero]
                # One reproducible integer witness for openness on each binomial.
                rng = Random(1111)
                c = {}
                for t in product(I, repeat=3):
                    old = tuple(labels[i][j] for i, j in enumerate(t))
                    c[t] = 0 if old in zeros else 1 if all(t) else rng.randint(2,70)
                verify(c)
                normalized, variables = normalize_binomial({t for t in TRIPLES if c[t]})
                cuts = []
                for rotation in I:
                    rotated = {t: normalized[t[rotation:]+t[:rotation]] for t in TRIPLES}
                    cuts.append(boundary_factors(rotated, variables))
                assert cuts[0] == cuts[1] == cuts[2]
                expected = {134217471: 6, 134217455: 5, 134151919: 4,
                            134184687: 3, 134217454: 4}
                assert len(cuts[0]) == expected[mask]
                print("Boundary factors", mask, zero, sorted(map(str,cuts[0])))
                count += 1
                print("PASS complementary binomial", mask, zero)
    assert count == 10


if __name__ == "__main__":
    main()
