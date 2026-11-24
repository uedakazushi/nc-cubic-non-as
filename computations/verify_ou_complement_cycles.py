#!/usr/bin/env python3
"""Exact generic cycle constructions on all 43 complementary-pencil divisors."""
from collections import Counter
from itertools import combinations, product
from random import Random

from sympy import Matrix, Poly, QQ, S, cancel, expand

from verify_ou_complement_line_conic import (
    boundary_factors, normalize_binomial, pencil, pencil_matrix, x, y,
)
from verify_ou_few_zero_cuts import SMOOTH_MASKS
from verify_ou_zero_supports import I, TRIPLES, support


U = {
    3: [(0,1,-1),(1,0,-1),(0,0,0)],
    4: [(-1,0,0,1),(0,0,0,0),(-1,1,0,0)],
}
V = {
    3: [(1,1,-1),(1,0,0),(1,0,0)],
    4: [(0,0,0,1),(0,0,1,0),(0,1,0,0)],
}
W = {
    3: [(1,1,0),(2,0,0),(2,0,0)],
    4: [(0,1,0,1),(1,0,1,0),(1,0,1,0)],
}


def actions():
    for rot in I:
        for row,col in product(((0,1,2),(0,2,1)),repeat=2):
            perm = {}
            for u in TRIPLES:
                v = (row[u[0]],col[u[1]],u[2])
                perm[u] = v[-rot:]+v[:-rot]
            yield rot,row[1],col[1],perm


def face_dets(c):
    faces = []
    for k in I:
        faces.append(Matrix(2,2,lambda i,j: c[(0,i+1,j+1) if k == 0 else
                                              (i+1,0,j+1) if k == 1 else (i+1,j+1,0)]).det())
    return faces


def rational_point(f, variables, avoid):
    rng = Random(219)
    for _ in range(4):
        for u in variables:
            if Poly(f,u).degree() != 1:
                continue
            at = {z:rng.randint(2,70) for z in variables if z != u}
            equation = Poly(f.subs(at),u)
            if equation.degree() == 1:
                at[u] = -equation.TC()/equation.LC()
                if at[u] and all(g.subs(at) != 0 for g in avoid):
                    assert f.subs(at) == 0
                    return at
    raise AssertionError(('no rational witness', f))


def model(pencil_forms, count):
    ell,p,q,r,s = pencil_forms
    basis = Matrix([p,r if count == 3 else s]).jacobian((x,y))
    if basis.det() == 0:
        return None
    ell,p,q,r,s = pencil_forms.jacobian((x,y))*basis.inv()*Matrix([x,y])
    lam, m = q.diff(x), ell.diff(y)
    if q != lam*x or not lam or not m:
        return None
    if count == 3:
        sigma,tau = s.diff(x),s.diff(y)
        e = s-lam*y
        if not all((sigma,tau,tau-lam,ell.subs({x:tau-lam,y:-sigma}))):
            return None
        return [
            Matrix([[ell,e,e],[1,1,lam],[x,y,s]]),
            Matrix([[1,0,0],[y,1,1],[x,1,1]]),
            Matrix([[m,(tau-lam)*x,(tau-lam)*x],[1,x+y,lam*x+tau*y],[0,1,tau]]),
        ]
    mu = r.diff(x)
    if r != mu*x or not mu or ell.subs({x:1,y:lam*mu}) == 0:
        return None
    e = y-lam*mu*x
    return [
        Matrix([[ell,e,e],[1,1,lam],[x,mu*x,y]]),
        Matrix([[m,1,x],[1,1,y+lam*x],[0,0,1]]),
        Matrix([[m,x,0],[1,x+y,1],[0,y,1]]),
        Matrix([[1,0,0],[y,1,1],[x,1,1]]),
    ]


def verify_fibre(c, matrices):
    n = len(matrices)
    for k in I:
        ck = Matrix(3,3,lambda i,j: c[i,j,k])
        kernels = []
        for h,mat in enumerate(matrices):
            m = ck.multiply_elementwise(mat)
            assert expand(m.det()) == 0
            cross = [m.row(i).cross(m.row(j)) for i,j in combinations(I,2)]
            gcd = Poly(0,x,y)
            for v in cross:
                for g in v:
                    gcd = gcd.gcd(Poly(g,x,y))
            assert gcd.total_degree() == 0
            v = next(v for v in cross if any(expand(g) != 0 for g in v))
            common = Poly(0,x,y)
            for g in v:
                common = common.gcd(Poly(g,x,y))
            kernel = [cancel(g/common.as_expr()) for g in v]
            for j,g in enumerate(kernel):
                if g:
                    assert {sum(mon)+V[n][j][h] for mon in Poly(g,x,y).monoms()} == {W[n][k][h]}
            kernels.append(kernel)
        for i in I:
            rows = []
            for mat,v in zip(matrices,kernels):
                polys = [Poly(mat[i,j]*v[j],x,y) for j in I]
                for mon in set().union(*(f.monoms() for f in polys)):
                    rows.append([f.coeff_monomial(mon) for f in polys])
            combinations_matrix = Matrix(rows)
            assert combinations_matrix.rank() == 2 and all(any(combinations_matrix[:,j]) for j in I)


def verify_factor(c, variables, factors, f):
    at = rational_point(f,variables,(factors-{f}) | set(face_dets(c)))
    fibre = {t:S(v).subs(at) for t,v in c.items()}
    for rot,row,col,perm in actions():
        numeric = {t:fibre[v] for t,v in perm.items()}
        forms = pencil(numeric)
        for n in (3,4):
            matrices = model(forms,n)
            if matrices is None:
                continue
            # [P,Q]=0 (also [P,R]=0 in the four-cycle). A complementary minor
            # of the three pencil equations computes each bracket. Divisibility
            # proves the identity on the whole irreducible divisor.
            eqs = pencil_matrix({t:c[v] for t,v in perm.items()})
            divisor = Poly(f,*variables,domain=QQ)
            for j in ((2,) if n == 3 else (2,3)):
                minor = eqs[:,[i for i in range(5) if i not in (1,j)]].det()
                assert Poly(minor,*variables,domain=QQ).rem(divisor).is_zero
            verify_fibre(numeric,matrices)
            return n,(rot,row,col)
    raise AssertionError(('no cycle chart',f))


def numerical_checks():
    for n in (3,4):
        blocks = (U[n],V[n],W[n])
        for i in I:
            for j in range(i,3):
                for a,b in product(blocks[i],blocks[j]):
                    difference = [z-v for z,v in zip(a,b)]*2
                    assert max(sum(difference[l:l+d]) for l in range(n)
                               for d in range(1,n)) <= 1
    # Exceptional divisors are total transforms. The only proximity is the
    # third centre on the fourth exceptional curve, whose strict class is E4-E3.
    collection = [(-2,1,1,1,1,1,1),(-1,1,1,1,0,0,0),(0,)*7]
    collection.extend(tuple(int(j == i) for j in range(7)) for i in (1,2,3))
    collection.extend(tuple(int(j == 0)-int(j == i) for j in range(7)) for i in (4,5,6))
    for first,second in combinations(collection,2):
        d = [i-j for i,j in zip(first,second)]
        assert 2+d[0]*(d[0]+3)-sum(z*(z-1) for z in d[1:]) == 0
        for coeff in (d,[k-i for k,i in zip((-3,1,1,1,1,1,1),d)]):
            strict = coeff[1:]
            strict[2] += strict[3]
            assert coeff[0] < 0 or (coeff[0] == 0 and min(strict) < 0)


def main():
    numerical_checks()
    counts = Counter()
    for masks in SMOOTH_MASKS.values():
        for mask in masks:
            zeros = set(TRIPLES)-set(support(mask))
            for zero in sorted(zeros):
                if any(all(i != j for i,j in zip(zero,t)) for t in zeros):
                    continue
                labels = [[i]+[j for j in I if j!=i] for i in zero]
                sup = {t for t in TRIPLES if tuple(labels[i][j] for i,j in enumerate(t)) not in zeros}
                c,variables = normalize_binomial(sup)
                factors = boundary_factors(c,variables)
                for f in sorted(factors,key=str):
                    n,chart = verify_factor(c,variables,factors,f)
                    counts[n] += 1
                    print('PASS',mask,zero,f,'cycle',n,'chart',chart,flush=True)
    assert counts == {3:29,4:14}
    print('PASS all 43 divisorial cases',dict(counts))


if __name__ == '__main__':
    main()
