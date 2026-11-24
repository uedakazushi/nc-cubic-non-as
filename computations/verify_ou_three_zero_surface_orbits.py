#!/usr/bin/env python3
"""Cover both remaining three-zero supports by twelve explicit rational surfaces."""
from functools import lru_cache
from itertools import product

from sympy import S, cancel, prod

from verify_ou_complement_line_conic import boundary_factors, normalize_binomial
from verify_ou_small_zero_corners import graph, residual_divisors
from verify_ou_zero_supports import I, TRIPLES, incidence, support


@lru_cache(None)
def chart(mask, zero):
    labels = [[i]+[j for j in I if j!=i] for i in zero]
    original = set(support(mask))
    sup = {w for w in TRIPLES if tuple(labels[i][w[i]] for i in I) in original}
    c,variables = normalize_binomial(sup)
    cube = sorted(w for w in sup if all(w))
    axes = sorted(w for w in sup if w.count(0)==2)
    rows = cube+axes+sorted(sup-set(cube+axes))
    matrix = incidence(rows)
    pivots = matrix.T.rref()[1]
    columns = matrix[list(pivots),:].rref()[1]
    weights = matrix[:,list(columns)]*matrix[list(pivots),list(columns)].inv()
    assert weights*matrix[list(pivots),:] == matrix and all(k.is_Integer for k in weights)
    coordinates = [next(w for w in sup if c[w]==v) for v in variables]
    monomials = [(w,[(rows[pivots[j]],-k) for j,k in enumerate(weights[rows.index(w),:]) if k])
                 for w in coordinates]
    return sup,c,variables,monomials


@lru_cache(None)
def reindexings(mask, source, target):
    src = chart(mask,source)[0]
    dest = chart(mask,target)[0]
    maps = []
    for rotation in I:
        for perms in product(((0,1,2),(0,2,1)),repeat=3):
            phi = {}
            for w in TRIPLES:
                v = tuple(p[i] for p,i in zip(perms,w))
                phi[w] = v[rotation:]+v[:rotation]
            if all((w in dest)==(phi[w] in src) for w in TRIPLES):
                maps.append(phi)
    return maps


def covered(mask, source, c, target, f, g):
    _,_,variables,monomials = chart(mask,target)
    for phi in reindexings(mask,source,target):
        values = {v:cancel(c[phi[w]]*prod(c[phi[h]]**k for h,k in powers))
                  for v,(w,powers) in zip(variables,monomials)}
        if all(cancel(q.subs(values,simultaneous=True))==0 for q in (f,g)):
            return True
    return False


def representatives(mask):
    result = []
    count = 0
    for zero in sorted(set(TRIPLES)-set(support(mask))):
        _,c,variables,_ = chart(mask,zero)
        for f in sorted(boundary_factors(c,variables),key=str):
            cc,params,_,cuts = residual_divisors(c,variables,f,require_single_chart=True)
            for g in cuts:
                count += 1
                # A rational graph at the generic point suffices for containment
                # in the closed target equations under the torus isomorphism.
                p,val = graph(g,params)
                generic = {w:cancel(S(a).subs(p,val)) for w,a in cc.items()}
                if any(covered(mask,zero,generic,target,rf,rg) for target,rf,rg in result):
                    continue
                graph(g,params,require_unit=True)
                result.append((zero,f,g))
                print('REP',mask,zero,f,g,flush=True)
    assert count == 42
    return result


def main():
    for mask,n in ((134151919,8),(134217454,4)):
        rows = representatives(mask)
        assert len(rows)==n
        print('PASS entire residual cover',mask,len(rows),flush=True)


if __name__ == '__main__':
    main()
