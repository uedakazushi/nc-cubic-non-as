#!/usr/bin/env python3
"""Determine every divisorial failure of the four-zero complementary conic charts."""
from collections import Counter
from itertools import permutations, product

from sympy import Matrix, Poly, S, cancel, eye, factor, factor_list, fraction, gcd, prod

from verify_ou_corner_comparison import CURVES, t
from verify_ou_five_zero_cycle_pencils import verify_family

from verify_ou_complement_cycles import actions, face_dets, model, verify_factor
from verify_ou_complement_line_conic import boundary_factors, pencil_matrix, x, y
from verify_ou_four_zero_boundary import tensor, u, v, w, q
from verify_ou_zero_supports import I, TRIPLES, incidence, support


CASES = (
    (134135535, (0,2,2), u, {v-1,w-1,q-1}),
    (134135535, (1,2,1), v, {u-1,w-1,q-1}),
    (134135535, (1,1,2), w, {u-1,v-1,q-1}),
    (134135535, (0,1,1), q, {u-1,v-1,w-1}),
    (134147823, (1,2,1), v, {u+w-1,q-1}),
    (134147823, (0,1,1), q, {v-1,u+v*w-v}),
    (134151407, (1,2,1), v, {u+w-1,w+q-1}),
    (134151918, (0,2,2), u, {q-1,w+1}),
    (134151918, (0,1,1), q, {u-1,w+1}),
)


def support_symmetry():
    """Cyclic changes of the three blocks already act transitively on unordered pairs."""
    mask = 134135535
    supported = support(mask)
    weights = incidence(supported)
    rows, cols = weights.T.rref()[1], weights.rref()[1]
    free = [i for i in range(len(supported)) if i not in rows]
    normalize = weights[:,cols] * weights[rows,cols].inv()
    indices = {word:i for i,word in enumerate(supported)}
    pairs = set()
    for rot in I:
        for labels in product(permutations(I),repeat=3):
            permuted = [tuple(labels[j][word[j]] for j in I) for word in supported]
            permuted = [word[rot:]+word[:rot] for word in permuted]
            if set(permuted) != set(supported):
                continue
            exponents = eye(len(supported))[[indices[w] for w in permuted],free]
            action = (exponents - normalize*exponents[rows,:])[free,:]
            assert all(v in (0,1) for v in action)
            assert all(sum(action.row(i)) == sum(action.col(i)) == 1 for i in range(4))
            pairs.add(tuple(sorted(i for i in range(4) if action[i,0]+action[i,1])))
    assert len(pairs) == 6


def solve_factor(f, variables):
    for v in variables:
        p = Poly(f,v)
        if p.degree() == 1:
            # The two coefficients have no common zero in the remaining
            # one-variable case; for the first cut use a constant unit.
            if len(variables) > 2 and p.LC().free_symbols:
                continue
            if len(variables) == 2:
                assert gcd(p.TC(),p.LC()).is_number
            return {v: cancel(-p.TC()/p.LC())}
    raise AssertionError(f)


def check_residual_curve(mask, values):
    if not all(values) or (mask == 134151918 and values[1] == 1):
        return
    for candidate_mask,row in CURVES:
        if candidate_mask != mask:
            continue
        # Only one ordered pair of fixed parameters is used for the first mask.
        variants = permutations(row) if mask == 134135535 else [row]
        for pattern in variants:
            coordinate = next(i for i,v in enumerate(pattern) if v == t)
            at = values[coordinate]
            if all(cancel(S(v).subs(t,at)-b) == 0 for v,b in zip(pattern,values)):
                return
    raise AssertionError(('uncovered curve', mask, values))


def verify_uniform(mask, c, fixed, variables, f, count, chart):
    first = solve_factor(f,variables)
    remaining = [v for v in variables if v not in first]
    perm = next(p for r,i,j,p in actions() if (r,i,j) == chart)
    coefficients = {word: cancel(c[p].subs(first)) for word,p in perm.items()}
    selector = eye(5)
    square = pencil_matrix(coefficients).col_join(selector[1,:]).col_join(selector[3 if count==3 else 4,:])
    forms = (square.inv()*Matrix([0,0,0,x,y])).applyfunc(factor)
    ell,_,Q,R,Sform = forms
    lam, m = Q.diff(x), ell.diff(y)
    units = [lam,m]
    if count == 3:
        sigma,tau = Sform.diff(x),Sform.diff(y)
        units.extend([sigma,tau,tau-lam,ell.subs({x:tau-lam,y:-sigma},simultaneous=True)])
    else:
        mu = R.diff(x)
        units.extend([mu,ell.subs({x:1,y:lam*mu},simultaneous=True)])
    raw_factors = units + face_dets(coefficients) + [v for v in coefficients.values() if v != 0]
    invert = set(remaining)
    for value in raw_factors:
        assert value != 0
        for z in fraction(cancel(value)):
            invert.update(p for p,_ in factor_list(z)[1])
    for value in forms.jacobian((x,y)):
        invert.update(p for p,_ in factor_list(fraction(cancel(value))[1])[1])
    # These basis choices are identities in the localized parameter ring.
    assert all(cancel(v) == 0 for v in pencil_matrix(coefficients)*forms)
    matrices = model(forms,count)
    assert matrices is not None
    verify_family(coefficients,matrices,prod(invert))
    values = [S(1) if p==fixed else first.get(p,p) for p in (u,v,w,q)]
    for g in invert:
        second = solve_factor(g,remaining)
        specialized = [cancel(z.subs(second)) for z in values]
        # A pole lies outside this graph chart. Its leading coefficient
        # on the first cut is a constant, so no original point is lost.
        check_residual_curve(mask,specialized)


def main():
    support_symmetry()
    counts = Counter()
    for mask, zero, fixed, expected in CASES:
        coefficients = tensor(mask, {fixed: 1})
        labels = [[j] + [i for i in I if i != j] for j in zero]
        c = {t: coefficients[tuple(labels[i][t[i]] for i in I)] for t in TRIPLES}
        variables = [i for i in (u,v,w,q) if i != fixed]
        assert c[0,0,0] == 0 and all(c[t] == 1 for t in TRIPLES if all(t))
        factors = boundary_factors(c, variables)
        assert factors == expected
        for f in sorted(factors, key=str):
            count, chart = verify_factor(c, variables, factors, f)
            if mask != 134135535 or (fixed == u and f == v-1):
                verify_uniform(mask,c,fixed,variables,f,count,chart)
            counts[count] += 1
            print("PASS", mask, zero, f, count, chart)
    assert counts == {3: 6, 4: 16}
    print("PASS all 22 four-zero divisor incidences", dict(counts))


if __name__ == "__main__":
    main()
