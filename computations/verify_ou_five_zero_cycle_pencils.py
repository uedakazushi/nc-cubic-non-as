#!/usr/bin/env python3
"""Verify the relative three-line and four-cycle models of the five-zero pencils."""
from itertools import combinations
from functools import reduce
from sympy import Matrix, Poly, QQ, S, cancel, fraction, gcd, groebner, symbols
from verify_ou_complement_cycles import V, W
from verify_ou_complement_line_conic import x, y, pencil_matrix
from verify_ou_zero_supports import I, TRIPLES

a = symbols("a")
def families():
    raw = [
        ([[0,1,1],[1,1,0],[1,1,-1]], [[1,0,1],[1,1,1],[0,1,1]], [[1,1,0],[1,1,1],[a,1,1]]),
        ([[0,1,1],[1,1,0],[1,1,a]], [[1,0,1],[1,1,1],[1,1,1]], [[1,1,0],[0,1,1],[1,1,1]]),
    ]
    for n, slices in enumerate(raw,3):
        c = {(i,j,k):S(slices[k][i][j]) for i,j,k in TRIPLES}
        if n == 3:
            lam, tau = 1/a, -1
            forms = Matrix([-y,x,lam*x,y,x-y])
            e = x-y-lam*y
            matrices = [
                Matrix([[-y,e,e],[1,1,lam],[x,y,x-y]]),
                Matrix([[1,0,0],[y,1,1],[x,1,1]]),
                Matrix([[-1,x,x],[1,y,tau*y-x],[0,1,tau]]),
            ]
            units = a
        else:
            c = {(i,j,k):c[k,i,j] for i,j,k in TRIPLES}
            lam = a/(a-1)
            forms = Matrix([y,x,lam*x,-lam*x,y])
            e = y+lam**2*x
            matrices = [
                Matrix([[y,e,e],[1,1,lam],[x,-lam*x,y]]),
                Matrix([[1,1,x],[1,1,y+lam*x],[0,0,1]]),
                Matrix([[1,x,0],[1,x+y,1],[0,y,1]]),
                Matrix([[1,0,0],[y,1,1],[x,1,1]]),
            ]
            units = a*(a-1)
        assert all(cancel(g) == 0 for g in pencil_matrix(c)*forms)
        yield c,matrices,units

def unit(f, allowed):
    variables = sorted(allowed.free_symbols, key=str)
    p = Poly(f,*variables,domain=QQ)
    if p.is_zero:
        return False
    return Poly(allowed, *variables, domain=QQ).rem(p.sqf_part()).is_zero

def numerator(v, allowed):
    n,d = fraction(cancel(v))
    assert unit(d,allowed), d
    return n

def empty_locus(entries, allowed, chart_variables=()):
    parameters = sorted(allowed.free_symbols, key=str)
    entries = [numerator(v, allowed) for v in entries if v != 0]
    if any(f.free_symbols.issubset(parameters) and unit(f, allowed) for f in entries):
        return
    inverse = symbols('inverse')
    ideal = groebner(entries + [inverse*allowed-1], inverse, *chart_variables,
                     *parameters, order="grevlex", domain=QQ)
    assert list(ideal) == [1]


def nowhere_common(entries, allowed):
    for at, variable in ((x, y), (y, x)):
        empty_locus([f.subs(at, 1) for f in entries], allowed, (variable,))


def no_parameter_zero(entries, allowed):
    empty_locus(entries, allowed)


def kernel_vector(matrix, allowed):
    crosses = [matrix.row(i).cross(matrix.row(j)) for i, j in combinations(I, 2)]
    nowhere_common([v for cross in crosses for v in cross], allowed)
    cross = next(cross for cross in crosses if any(cancel(v) != 0 for v in cross))
    denominator = reduce(lambda l, r: l*r, (fraction(cancel(v))[1] for v in cross), S.One)
    entries = [cancel(denominator*v) for v in cross]
    common = reduce(gcd, entries, 0)
    vector = [cancel(v/common) for v in entries]
    nowhere_common(vector, allowed)
    assert all(cancel(f) == 0 for f in matrix*Matrix(vector))
    return vector


def verify_family(c, matrices, allowed):
    count = len(matrices)
    for k in I:
        coefficient_slice = Matrix(3, 3, lambda i, j: c[i, j, k])
        kernels = []
        for h, matrix in enumerate(matrices):
            relation = coefficient_slice.multiply_elementwise(matrix)
            assert cancel(relation.det()) == 0
            vector = kernel_vector(relation, allowed)
            for j, v in enumerate(vector):
                if v != 0:
                    degrees = {sum(mon)+V[count][j][h] for mon in Poly(v, x, y).monoms()}
                    assert degrees == {W[count][k][h]}
            kernels.append(vector)

        for i in I:
            rows = []
            for matrix, vector in zip(matrices, kernels):
                polys = [Poly(matrix[i, j]*vector[j], x, y) for j in I]
                for mon in set().union(*(p.monoms() for p in polys)):
                    rows.append([p.coeff_monomial(mon) for p in polys])
            compositions = Matrix(rows)
            assert all(cancel(v) == 0 for v in compositions*Matrix([c[i, j, k] for j in I]))
            minors = [compositions.extract(rs, js).det() for rs in combinations(range(len(rows)), 2)
                      for js in combinations(I, 2)]
            no_parameter_zero(minors, allowed)
            for j in I:
                no_parameter_zero(compositions[:, j], allowed)


if __name__ == "__main__":
    for coefficients, matrices, allowed in families():
        verify_family(coefficients, matrices, allowed)
        print(f"PASS {len(matrices)} components on the full open {allowed} != 0")
