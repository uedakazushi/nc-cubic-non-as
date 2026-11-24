#!/usr/bin/env python3
"""Exact cycle-chart divisors and regular quadratic corner covers at their witnesses."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from random import Random

from sympy import Matrix, Poly, QQ, S, cancel, factor_list, groebner, symbols

from verify_ou_complement_cycles import actions, face_dets
from verify_ou_complement_line_conic import boundary_factors, normalize_binomial, pencil_matrix
from verify_ou_corner_comparison import CornerCertificate
from verify_ou_few_zero_cuts import SMOOTH_MASKS
from verify_ou_zero_supports import I, TRIPLES, support


def graph(polynomial, variables, require_unit=False):
    """Write an irreducible hypersurface as a rational graph; track vertical fibres."""
    choices = []
    for v in variables:
        p = Poly(polynomial,v)
        if p.degree() != 1:
            continue
        unit = len(Poly(p.LC(),*(w for w in variables if w!=v)).terms()) == 1
        if require_unit and not unit:
            continue
        value = cancel(-p.TC()/p.LC())
        choices.append((not unit, len(str(value)), str(v), v, value))
    assert choices, (polynomial,variables,require_unit)
    return min(choices)[-2:]


def chart_conditions(coefficients):
    matrix = pencil_matrix(coefficients)

    @lru_cache(None)
    def bracket(i,j):
        if i > j:
            return -bracket(j,i)
        cols = [k for k in range(5) if k not in (i,j)]
        return (-1)**(i+j+1) * matrix[:,cols].det()

    l,p,q,r,s = range(5)
    b = bracket
    if cancel(b(p,q)) != 0:
        return
    d = b(p,r)
    yield (d,b(q,r),b(p,l),b(s,r),b(p,s),b(l,s)*d-b(q,r)*b(l,r))
    if cancel(b(p,r)) == 0:
        d = b(p,s)
        yield (d,b(q,s),b(r,s),b(p,l),b(l,s)*d**2+b(q,s)*b(r,s)*b(p,l))


def factors(value):
    return {f for part in cancel(value).as_numer_denom() for f,_ in factor_list(part)[1]}


def residual_divisors(coefficients, variables, divisor, require_single_chart=False):
    pivot, value = graph(divisor,variables,require_unit=True)
    params = [p for p in variables if p != pivot]
    c = {t: cancel(S(z).subs(pivot,value)) for t,z in coefficients.items()}
    units = set(params) | factors(value)
    for determinant in face_dets(c):
        assert determinant != 0
        units.update(factors(determinant))
    failures = []
    for _,_,_,perm in actions():
        for entries in chart_conditions({t:c[h] for t,h in perm.items()}):
            if all(cancel(z) != 0 for z in entries):
                failures.append(set().union(*(factors(z) for z in entries))-units)
    assert failures
    common = set.intersection(*failures)
    if require_single_chart:
        assert common in failures
    return c,params,units,sorted(common,key=str)


def witness(c, params, units, divisor):
    pivot, value = graph(divisor,params)
    rng = Random(871)
    for _ in range(16):
        at = {v:S(rng.randint(2,17)) for v in params if v != pivot}
        at[pivot] = value.subs(at)
        if not at[pivot].is_Rational or not all(h.subs(at) != 0 for h in units):
            continue
        assert divisor.subs(at) == 0
        return {t:QQ.from_sympy(v.subs(at)) for t,v in c.items()}
    raise AssertionError(('no rational witness',divisor))


def verify_quadratic_cover(cert, rotation, vertex):
    loops = cert.loops(rotation,vertex)
    rows,free = cert.products(rotation,[a+b[1:] for a,b in product(loops,repeat=2)])
    table = cert.matrix([[row.get(k,QQ.zero) for k in free] for row in rows])
    relations = table.transpose().nullspace().to_list()
    assert len(relations) == 3 and len(free) == 6

    # H tensor R and R tensor H in the 27-dimensional space of cubic words.
    left = [[relations[j][3*v+w] if u == h else QQ.zero for u,v,w in TRIPLES]
            for h,j in product(I,repeat=2)]
    right = [[relations[j][3*u+v] if w == h else QQ.zero for u,v,w in TRIPLES]
             for j,h in product(I,repeat=2)]
    syzygies = cert.matrix(left+right).transpose().nullspace().to_list()
    assert len(syzygies) == 1
    for entries in (syzygies[0][:9],syzygies[0][9:]):
        assert cert.matrix([entries[k:k+3] for k in (0,3,6)]).det()

    coords = symbols('x:3')
    m = Matrix([[sum(coords[i]*QQ.to_sympy(row[3*i+j]) for i in I)
                 for j in I] for row in relations])
    minors = [m.extract(a,b).det() for a in combinations(I,2) for b in combinations(I,2)]
    for j in I:
        basis = groebner([z.subs(coords[j],1) for z in minors],
                         *(coords[:j]+coords[j+1:]),domain=QQ)
        assert tuple(basis) == (1,)


def verify_point(coefficients):
    # All 27 coefficients, including zero entries, are passed explicitly.
    cert = CornerCertificate(0,coefficients,QQ)
    cert.validate_paths()
    cert.check_generation()
    for rotation,vertex in product(I,repeat=2):
        verify_quadratic_cover(cert,rotation,vertex)
    for fn in (cert.quotient,cert.unit,cert.paths.space,cert.paths.word):
        fn.cache_clear()


def main():
    counts = Counter()
    divisors = 0
    for masks in SMOOTH_MASKS.values():
        for mask in masks:
            zeros = set(TRIPLES)-set(support(mask))
            for zero in sorted(zeros):
                if any(all(i!=j for i,j in zip(zero,t)) for t in zeros):
                    continue
                labels = [[i]+[j for j in I if j!=i] for i in zero]
                sup = {t for t in TRIPLES if tuple(labels[i][j] for i,j in enumerate(t)) not in zeros}
                c,variables = normalize_binomial(sup)
                for f in sorted(boundary_factors(c,variables),key=str):
                    reduced,params,units,remaining = residual_divisors(c,variables,f)
                    for g in remaining:
                        verify_point(witness(reduced,params,units,g))
                        print('PASS quadratic covers',mask,zero,f,g,flush=True)
                    counts[mask] += len(remaining)
                    divisors += 1
    assert divisors == 43
    assert counts == {134217471:30,134217455:42,134151919:42,134184687:6,134217454:42}
    print('PASS all 162 codimension-three incidences',dict(counts))


if __name__ == '__main__':
    main()
