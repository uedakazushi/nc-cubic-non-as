#!/usr/bin/env python3
"""All specializations of the complementary binomial on three-zero support 134184687."""
from itertools import product

from sympy import GF, Matrix, Poly, QQ, S, cancel, factor_list, primerange, symbols

from verify_ou_complement_cycles import actions, face_dets
from verify_ou_complement_line_conic import boundary_factors, normalize_binomial
from verify_ou_corner_comparison import CornerCertificate, t
from verify_ou_seven_zero_center import CyclicPaths
from verify_ou_small_zero_corners import chart_conditions, factors, graph
from verify_ou_zero_supports import I, TRIPLES, support


MASK = 134184687


class PencilCertificate(CornerCertificate):
    def relation_normals(self, table):
        # Clear each row over the coefficient polynomial ring before the determinant.
        relations = table.transpose().nullspace(divide_last=True).to_list()
        assert len(relations) == 3
        coords = symbols('X:3')
        rows = []
        for row in relations:
            denominator = row[0].denom.ring.one
            for v in row:
                denominator = denominator.lcm(v.denom)
            integral = [v.numer*denominator.exquo(v.denom) for v in row]
            common = integral[0]
            for v in integral[1:]:
                common = common.gcd(v)
            entries = [v.exquo(common).as_expr() for v in integral]
            rows.append([sum(coords[i]*entries[3*i+j] for i in I) for j in I])
        determinant = Matrix(rows).det(method='domain-ge')
        if determinant != 0:
            for f,_ in factor_list(determinant)[1]:
                form = Poly(f,*coords)
                if form.total_degree() == 1:
                    yield [self.field.from_sympy(form.coeff_monomial(x)) for x in coords]


def check_corners(coefficients, field):
    if field.is_FiniteField:
        paths = CyclicPaths(0,field,coefficients)
        central = [(3,z) for z in paths.center(3)]
        if len(central) > 1:
            for n in (5,6,7):
                if all(a == b for a,b in (paths.product_rank(n,r,i,j,central) for r,i,j in TRIPLES)):
                    return set()
    cls = PencilCertificate if field.is_FractionField else CornerCertificate
    cert = cls(0,coefficients,field)
    cert.validate_paths()
    cert.check_generation()
    for r,i in product(I,repeat=2):
        cert.check_corner(r,i)
    divisors = cert.exceptional.copy()
    for func in (cert.quotient,cert.unit,cert.paths.space,cert.paths.word):
        func.cache_clear()
    return divisors


def verify_point(c, units, h):
    p = Poly(h,t,domain=QQ)
    assert p.degree() > 0 and p.is_irreducible
    for value in units:
        if any(Poly(w,t,domain=QQ).rem(p).is_zero for w in cancel(value).as_numer_denom()):
            return 'outside'
    if p.degree() == 1:
        root = -p.TC()/p.LC()
        try:
            check_corners({w:QQ.from_sympy(v.subs(t,root)) for w,v in c.items()},QQ)
            return 'rational'
        except AssertionError:
            pass
    for prime in primerange(2,150):
        field = GF(prime)
        if any(int(k.q)%prime == 0 for k in p.all_coeffs()):
            continue
        for root in range(prime):
            if p.eval(root)%prime or not p.diff().eval(root)%prime:
                continue
            values = [v.subs(t,root) for v in units]
            if any(not v.is_Rational or int(v.q)*int(v.p)%prime == 0 for v in values):
                continue
            at = {w:v.subs(t,root) for w,v in c.items()}
            if any(not v.is_Rational or int(v.q)%prime == 0 for v in at.values()):
                continue
            try:
                check_corners({w:field(int(v.p))/field(int(v.q)) for w,v in at.items()},field)
                return prime,root
            except AssertionError:
                continue
    raise AssertionError(('no verified fibre',h))


def rational_curve(f, ps):
    x,y = ps
    # The coordinate changes are isomorphisms of the two-dimensional torus.
    for change in ({x:x,y:y},{x:x/y,y:y},{x:x,y:y/x}):
        h = cancel(f.subs(change,simultaneous=True)).as_numer_denom()[0]
        for pivot,param in ((x,y),(y,x)):
            poly = Poly(h,pivot)
            if poly.degree() == 1:
                leading = Poly(poly.LC(),param)
                constant = Poly(poly.TC(),param)
                assert leading.gcd(constant).degree() == 0
                v = cancel(-poly.TC()/poly.LC())
                sub = {z:cancel(a.subs(pivot,v).subs(param,t)) for z,a in change.items()}
                assert cancel(f.subs(sub,simultaneous=True)) == 0
                return sub
    raise AssertionError(('no torus graph',f))


def apply_sub(c, sub):
    return {w:cancel(S(v).subs(sub,simultaneous=True)) for w,v in c.items()}


def parametric_check(c, params):
    field = QQ.frac_field(*params)
    return check_corners({w:field.from_sympy(v) for w,v in c.items()},field)


def check_surface(c, ps, units):
    labels = ((1,0,2),(0,1,2),(2,0,1))
    changed = {w:c[tuple(p[j] for p,j in zip(labels,w))] for w in TRIPLES}
    cuts = parametric_check(changed,ps)-units
    assert len(cuts) in (4,5,6)
    point_count = 0
    for f in sorted(cuts,key=str):
        sub = rational_curve(f,ps)
        cc = apply_sub(c,sub)
        local_units = [cancel(v.subs(sub,simultaneous=True)) for v in units]
        assert all(local_units)
        excluded = set().union(*(factors(v) for v in local_units))
        remaining = parametric_check(cc,(t,))-excluded
        for h in sorted(remaining,key=str):
            where = verify_point(cc,local_units,h)
            print('PASS closed fibre',f,h,where,flush=True)
            point_count += 1
    return len(cuts),point_count


def main():
    zero = (0,2,2)
    order = [[i]+[j for j in I if j!=i] for i in zero]
    sup = set(support(MASK))
    sup = {w for w in TRIPLES if tuple(order[i][w[i]] for i in I) in sup}
    c,vs = normalize_binomial(sup)
    a,b,d,e = vs
    boundary = (a+b-1, a+b*d-b*e-1, a-b*e+e)
    assert boundary_factors(c,vs) == set(boundary)
    # A single expanded-cycle chart on each divisor fails only on intersections
    # with the other two. This checks set-theoretic coverage, including intersections.
    for f in boundary:
        p,v = graph(f,vs,require_unit=True)
        cc = apply_sub(c,{p:v})
        us = (set(vs)-{p}) | factors(v) | set().union(*(factors(w) for w in face_dets(cc)))
        expected = set().union(*(factors(h.subs(p,v)) for h in boundary if h!=f))-us
        assert len(expected)==2
        all_failures = []
        for _,_,_,perm in actions():
            for entries in chart_conditions({w:cc[h] for w,h in perm.items()}):
                if all(cancel(z)!=0 for z in entries):
                    all_failures.append(set().union(*(factors(z) for z in entries))-us)
        assert expected in all_failures

    surfaces = [
        ({a:1-b,d:e+1},(b,e)),
        ({a:1-b,e:-1},(b,d)),
        ({e:b*d-1,a:1-b*d+b*(b*d-1)},(b,d)),
    ]
    counts = []
    for (i,j),(sub,ps) in zip(((0,1),(0,2),(1,2)),surfaces):
        for k in (i,j):
            assert cancel(boundary[k].subs(sub,simultaneous=True)) == 0
        cc = apply_sub(c,sub)
        units = set().union(*(factors(v) for v in cc.values() if v))
        units.update(set().union(*(factors(v) for v in face_dets(cc))))
        counts.append(check_surface(cc,ps,units))
        print('PASS entire double-pencil surface',sub,counts[-1],flush=True)
    assert [n for n,_ in counts] == [6,5,4]
    print('PASS support 134184687',counts,flush=True)


if __name__ == '__main__':
    main()
