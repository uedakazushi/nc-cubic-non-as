#!/usr/bin/env python3
"""Central cubic quotients and three-generator corner covers of four-zero pencils."""
from functools import lru_cache, reduce
from itertools import product

from sympy import (GF, Matrix, Poly, QQ, S, cancel, factor, factor_list, fraction,
                   gcd, lcm, symbols)
from sympy.polys.matrices import DomainMatrix

from verify_ou_five_zero_supports import expected
from verify_ou_four_zero_boundary import params, tensor
from verify_ou_seven_zero_center import CyclicPaths
from verify_ou_zero_supports import I, TRIPLES


t = symbols('t')
CURVES = (
    (134135535, (1,1,1,t)), (134135535, (1,1,t,t)),
    (134147823, (1-t,1,t,(t-1)/t)), (134147823, (-t,1,t,1)),
    (134147823, (1-t,1,t,1)), (134147823, ((t-1)*(1-t)/t,(t-1)/t,t,1)),
    (134151407, (1-t,1,t,(1-t)/t)), (134151407, (1-t,1,t,1-t)),
    (134151407, (t/(1-t),1,1-t,t)), (134151918, (1,t,-1,1-t)),
    (134151918, (1,t,-1,1)), (134151918, (1,t+1,t,1)),
    (134151918, (1-t,t,-1,1)),
)


class CornerCertificate:
    def __init__(self, mask, coefficients, field):
        self.field = field
        self.paths = CyclicPaths(mask, field, coefficients)
        self.exceptional = set()
        centers = self.paths.center(3)
        assert len(centers) == 1
        terms = centers[0]
        if field.is_FiniteField:
            self.center = terms
            return
        values = [field.to_sympy(v) for v in terms.values()]
        denominator = reduce(lcm, (fraction(v)[1] for v in values))
        numerators = [cancel(v * denominator) for v in values]
        common = reduce(gcd, numerators)
        self.center = {w: field.from_sympy(cancel(v / common)) for w,v in zip(terms,numerators)}

    @lru_cache(None)
    def unit(self, value):
        assert value
        if self.field.is_FractionField:
            for polynomial in (value.numer, value.denom):
                self.exceptional.update(p.as_expr() for p,_ in
                                        Poly(polynomial.as_expr(), t, domain=QQ).factor_list()[1])

    def regular(self, value):
        if self.field.is_FractionField and value:
            self.unit(self.field(value.denom))

    def eliminate(self, rows, width):
        """Pivot rows with recorded nonzero divisors and reduce modulo their span."""
        pivots = {}

        def subtract(row, col):
            multiplier = row[col]
            for k, coefficient in pivots[col].items():
                v = row.get(k, self.field.zero) - multiplier * coefficient
                if v:
                    row[k] = v
                else:
                    row.pop(k, None)

        for original in rows:
            row = original.copy()
            while row:
                col = min(row)
                if col in pivots:
                    subtract(row, col)
                else:
                    divisor = row[col]
                    self.unit(divisor)
                    pivots[col] = {k: v/divisor for k,v in row.items()}
                    break

        def reduce_row(original):
            row = original.copy()
            for col in sorted(pivots):
                if col in row:
                    subtract(row, col)
            return row

        return [k for k in range(width) if k not in pivots], reduce_row

    @lru_cache(None)
    def quotient(self, degree, rotation, start, end):
        paths = self.paths
        rows = []
        for tail in paths.space(degree-3, rotation, start, end)[0]:
            row = {}
            for (turn, loop), coefficient in self.center.items():
                if turn == rotation and loop[0] == start:
                    for k, v in paths.word(turn, loop + tail[1:]).items():
                        row[k] = row.get(k, self.field.zero) + coefficient*v
            rows.append({k: v for k,v in row.items() if v})
        return self.eliminate(rows, len(paths.space(degree, rotation, start, end)[0]))

    def loops(self, rotation, vertex):
        basis = self.paths.space(3, rotation, vertex, vertex)[0]
        free, _ = self.quotient(3, rotation, vertex, vertex)
        assert len(free) == 3
        return [basis[k] for k in free]

    def products(self, rotation, words):
        w = words[0]
        free, reduce_row = self.quotient(len(w)-1, rotation, w[0], w[-1])
        return [reduce_row(self.paths.word(rotation, word)) for word in words], free

    def check_generation(self):
        for r,i,j in TRIPLES:
            words = [tail + loop[1:] for tail in self.paths.space(2,r,i,j)[0]
                     for loop in self.loops((r+5) % 3, j)]
            rows, free = self.products(r, words)
            missing, _ = self.eliminate(rows, max(free)+1)
            assert set(missing).isdisjoint(free)

    def matrix(self, rows):
        return DomainMatrix.from_dok({(i,j): v for i,row in enumerate(rows)
                                      for j,v in enumerate(row) if v},
                                     (len(rows),len(rows[0])), self.field)

    def normal_form(self, base, normal):
        """Verify relations of a skew polynomial extension in the associated graded."""
        field = self.field
        first = next(i for i in I if normal[i])
        rows = [normal] + [[field(int(j == k)) for k in I] for j in I if j != first]
        change = self.matrix([[rows[i][a]*rows[j][b] for a,b in product(I,repeat=2)]
                              for i,j in product(I,repeat=2)])
        table = change * base
        cols = list(range(table.shape[1]))
        first_row = table.extract([0,1,2], cols)
        pivots = first_row.rref()[1]
        if len(pivots) != 3:
            return False
        minor = first_row.extract(I,pivots)
        left_products = table.extract([3,6],cols)
        coeff = left_products.extract([0,1],pivots) * minor.inv()
        if not (coeff * first_row - left_products).is_zero_matrix:
            return False
        action = coeff.extract([0,1],[1,2])
        if not action.det():
            return False
        kernel = first_row.vstack(table.extract([4,5,7,8],cols)).transpose().nullspace()
        if kernel.shape[0] != 1:
            return False
        entries = kernel.to_list()[0]
        relation = self.matrix([entries[3:5],entries[5:7]])
        if not relation.det():
            return False
        mapped = action.transpose() * relation * action
        nonzero = next((i,j) for i,j in product(range(2),repeat=2) if relation.to_list()[i][j])
        i,j = nonzero
        scalar = mapped.to_list()[i][j] / relation.to_list()[i][j]
        if not (mapped - relation.scalarmul(scalar)).is_zero_matrix:
            return False
        for value in list(normal) + sum(coeff.to_list(),[]) + entries:
            self.regular(value)
        for value in (normal[first], minor.det(), action.det(), relation.det(), scalar):
            self.unit(value)
        return True

    def check_corner(self, rotation, vertex):
        loops = self.loops(rotation,vertex)
        rows, free = self.products(rotation,[a + b[1:] for a,b in product(loops,repeat=2)])
        field = self.field
        table = self.matrix([[row.get(k,field.zero) for k in free] for row in rows])
        candidates = [[field(int(i == j)) for j in I] for i in I]
        if any(self.normal_form(table, normal) for normal in candidates):
            return
        if field.is_FiniteField:
            for i in I:
                for tail in product(range(int(field.mod)),repeat=2-i):
                    normal = [field.zero]*i + [field.one] + list(map(field,tail))
                    if self.normal_form(table,normal):
                        return
            raise AssertionError((rotation,vertex,'no normal linear form'))
        relations = table.transpose().nullspace().to_list()
        assert len(relations) == 3
        coordinates = symbols('X:3')
        forms = []
        for row in relations:
            poly = [field.to_sympy(v) for v in row]
            forms.append(Matrix([coordinates]) * Matrix(3,3,poly))
        determinant = factor(Matrix.vstack(*forms).det(method='domain-ge'))
        if determinant != 0:
            for polynomial,_ in factor_list(determinant,*coordinates)[1]:
                form = Poly(polynomial,*coordinates)
                if form.total_degree() == 1:
                    normal = [field.from_sympy(form.coeff_monomial(x)) for x in coordinates]
                    if self.normal_form(table,normal):
                        return
        raise AssertionError((rotation,vertex,'no normal linear form'))

    def validate_paths(self):
        for n in range(2,7):
            for r,i,j in TRIPLES:
                basis, maps = self.paths.space(n,r,i,j)
                assert len(basis) == expected(n,i,j)
                for row in maps.values():
                    for v in row.values():
                        self.regular(v)


def check(mask, values, field):
    coefficients = tensor(mask, dict(zip(params, values)))
    if field.is_FiniteField:
        c = {p: field(int(v.p)) / field(int(v.q)) for p,v in coefficients.items()}
        paths = CyclicPaths(mask, field, c)
        central = [(3,z) for z in paths.center(3)]
        if len(central) > 1:
            for degree in (5,6,7):
                if all(rank == dim for rank,dim in
                       (paths.product_rank(degree,r,i,j,central) for r,i,j in TRIPLES)):
                    return set()
    else:
        c = {p: field.from_sympy(v) for p,v in coefficients.items()}
    certificate = CornerCertificate(mask, c, field)
    certificate.check_generation()
    for r,i in product(I,repeat=2):
        certificate.check_corner(r,i)
    certificate.validate_paths()
    return certificate.exceptional


def admissible(mask, values, prime=None):
    """Check the coefficient torus and the v != 1 condition in the last support."""
    if any(not value.is_Rational or not value for value in values):
        return False
    if prime and any(int(value.p)*int(value.q) % prime == 0 for value in values):
        return False
    return mask != 134151918 or (values[1] != 1 if not prime
                                else int((values[1]-1).p) % prime != 0)


def check_closed_point(mask, expressions, divisor):
    polynomial = Poly(divisor,t,domain=QQ)
    required = list(expressions) + ([expressions[1]-1] if mask == 134151918 else [])
    for f in required:
        if any(Poly(z,t,domain=QQ).rem(polynomial).is_zero for z in fraction(cancel(f))):
            return 'outside open'
    if polynomial.degree() == 1:
        root = -polynomial.TC()/polynomial.LC()
        values = [f.subs(t,root) for f in expressions]
        try:
            check(mask,values,QQ)
            return 'Q'
        except AssertionError:
            pass
    for prime in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71):
        if int(polynomial.LC()) % prime == 0:
            continue
        for h in range(prime):
            if polynomial.eval(h) % prime or polynomial.diff().eval(h) % prime == 0:
                continue
            values = [f.subs(t,h) for f in expressions]
            if admissible(mask,values,prime):
                try:
                    check(mask,values,GF(prime))
                    return (prime,h)
                except AssertionError:
                    continue
    raise AssertionError((mask,expressions,divisor,'no verified fiber'))


def check_curve(index):
    mask, row = CURVES[index]
    values = list(map(S,row))
    divisors = check(mask,values,QQ.frac_field(t))
    for f in sorted(divisors,key=str):
        result = check_closed_point(mask,values,f)
        print("PASS exceptional factor", index, f, result, flush=True)
    print("PASS full pencil",index,mask,row,flush=True)


if __name__ == '__main__':
    for index in range(len(CURVES)):
        check_curve(index)
