#!/usr/bin/env python3
"""Exact path ranks and simultaneous three-line models on three-zero strata."""
from itertools import combinations, product
from sympy import Matrix, QQ, Poly, cancel, expand, groebner, symbols

from verify_ou_five_zero_supports import expected
from verify_ou_four_zero_triangle import triangle_components
from verify_ou_seven_zero_curves import determinant_matrix
from verify_ou_zero_supports import I, TRIPLES, normalize, path_dim


a, b, x, y, z, h = symbols('a b x y z h')


def tensor(mask, values):
    c, field, free = normalize(mask)
    assert len(free) == len(values) == 5
    return {t: values[free.index(t)] if t in free else field.to_sympy(v) for t, v in c.items()}


def path_checks():
    field = QQ.frac_field(a, b)
    for mask, values, n, r, i, j, dim in [
        (134184687, (1-(b+1)/a, a, -b, b+1, b), 5, 1, 0, 2, 8),
        (134216943, (a, b, b, a*b, b), 6, 0, 1, 1, 11),
        (134151919, ((b-1+a)/b, a, a, 1-a, b), 5, 0, 2, 0, 8),
        (134151919, (a, b, a*b-b+1, 1-b, b), 5, 2, 0, 0, 8),
        (134151919, ((a+b-1)/a, a, b, 1-b, b), 5, 1, 0, 2, 8),
    ]:
        c = {t: field.from_sympy(f) for t, f in tensor(mask, values).items()}
        assert path_dim(c, n, r, i, j, field) == dim > expected(n, i, j)
        print('PASS path quotient', mask, n, r, i, j, dim)


def unit(polys, domain, *variables):
    nums = [cancel(f).as_numer_denom()[0] for f in polys]
    assert list(groebner(nums + [h*domain-1], *variables, h, a, b)) == [1]


def triangle_check(mask):
    if mask == 134216943:
        w, values = a, (a, b, 1-a*b, 1-b, b)
        domain = a*b*(1-a*b)*(1-b)*(1-(a+1)*b)
        kernels = [
            [[-x,-1,1], [-a*x-z,-1,1], [b*x-b*z-x,-1,1]],
            [[-1,a*x,a*y], [0,-1,1], [1-a*b-b,a*(a*b*(x+y)-x),a*(b*(x+y)-y)]],
            [[-1,1,1], [a*y-z,-y,z], [-1,1,1]],
        ]
        w_degree = [(1,1,0), (1,0,1), (1,1,0)]
    else:
        w, values = 1-a, (a, b, b, b-a*b, b)
        domain = a*b*(1-a)*(b-1)*(1-(2-a)*b)
        kernels = [
            [[-x,-1,1], [-z,-1,1], [b*((a-1)*x-z),-1,1]],
            [[1,(a-1)*x,(a-1)*y], [-1,y-a*(x+y),x], [0,-1,1]],
            [[-1,1,1], [-1,1,1], [b*((1-a)*y-z),b*(y+z)-y,(a-1)*b*(y+z)+z]],
        ]
        w_degree = [(1,1,0), (1,1,0), (1,0,1)]
    original = tensor(mask, values)
    coeff = {t: original[t[-1:] + t[:-1]] for t in TRIPLES}
    slices = [Matrix(3,3,lambda i,j: coeff[i,j,k]) for k in I]
    matrices = triangle_components(w,x,y,z)
    vs = [(x,z), (x,y), (y,z)]
    p = Matrix([w*(x+y),x+y,x+w*y-z,x,y,z])
    assert all(cancel(f)==0 for f in determinant_matrix(coeff)*p)
    assert expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]-(x+y)*(x-z)*(w*y-z))==0
    v_degrees = [(0,1,0), (1,0,0), (1,0,0)]
    for l, mat in enumerate(matrices):
        for k in I:
            v = Matrix(kernels[l][k])
            m = slices[k].multiply_elementwise(mat)
            assert all(cancel(f)==0 for f in m*v)
            for j in I:
                if v[j]:
                    assert {sum(power) for power in Poly(v[j],*vs[l]).monoms()} == {
                        w_degree[k][l]-v_degrees[j][l]}
            minors = [m[rows,cols].det() for rows,cols in product(combinations(I,2),repeat=2)]
            for t in vs[l]:
                other = next(u for u in vs[l] if u!=t)
                unit([f.subs(t,1) for f in v], domain, other)
                unit([f.subs(t,1) for f in minors], domain, other)
    for k in I:
        v = Matrix(kernels[0][k])
        for i in (1,2):
            prods = [expand(matrices[0][i,j]*v[j]) for j in (1,2)]
            unit([Matrix([[f.coeff(x),f.coeff(z)] for f in prods]).det()],domain)
            unit(Poly(matrices[0][i,0]*v[0],x,z).coeffs(),domain)
        assert matrices[0][0,0]==0 and all(matrices[0][0,j]*v[j]!=0 for j in (1,2))
        other_coeffs = [q for l in (1,2) for q in
                        Poly(matrices[l][0,0]*kernels[l][k][0],*vs[l]).coeffs()]
        unit(other_coeffs,domain)
    print('PASS three-line family',mask)


if __name__ == '__main__':
    path_checks()
    for support_mask in (134216943,134151919):
        triangle_check(support_mask)
