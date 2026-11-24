#!/usr/bin/env python3
"""Check simultaneous failure of the three determinant cuts over Q (Singular 4.4)."""
from itertools import combinations, product
import argparse
from collections import Counter
import os
import subprocess

from sympy import Matrix, factor_list, prod, symbols

from verify_ou_few_zero_cuts import SMOOTH_MASKS
from verify_ou_seven_zero_curves import determinant_matrix
from verify_ou_zero_supports import I, TRIPLES, normalize


def complement_factors(coeff):
    """Non-unit complementary cube binomials for the isolated zeros."""
    factors = set()
    for zero in (t for t in TRIPLES if coeff[t] == 0):
        axes = [[a for a in I if a != b] for b in zero]
        cube = [coeff[x, y, z] for x, y, z in product(*axes)]
        if not all(cube):
            continue
        delta = prod(cube[i] for i in (0, 3, 5, 6)) - prod(cube[i] for i in (1, 2, 4, 7))
        parts = [f for f, _ in factor_list(delta)[1] if not f.is_Symbol]
        assert len(parts) == 1 and delta != 0
        factors.update(parts)
    return factors


def failure_ideals(coeff):
    matrix = determinant_matrix(coeff)
    cols = range(6)
    for i, j in product((0, 3, 4), (1, 2, 5)):
        other = sorted(set(cols) - {i, j})
        yield [matrix[:, t].det() for t in combinations(other, 3)]
        yield [matrix[t, [i, j]].det() for t in combinations(I, 2)]
    for k in I:
        if all(matrix[k, j] != 0 for j in cols):
            reciprocal = Matrix([sign * prod(matrix[k, j] for j in cols if j != t)
                                 for t, sign in enumerate((1, -1, -1, 1, 1, -1))])
            yield list(matrix * reciprocal)


# Each list describes the closed union that can remain off the complementary
# binomials. On the one- and two-zero tori only the dimension bound is tested.
a, b, c, d, e = symbols('u0:5')
STRATA = {
    134217454: [],
    134184687: [(d-e-1, c+e, a*b-b+e+1)],
    134216943: [(c-e, b-e, a*e-d), (d+e-1, b-e, a*e+c-1)],
    134151919: [(c+d-1, b+d-1, a*e+d-e), (c-e, b-e, a*e+d-e),
                (d+e-1, b-e, a*e-c-e+1), (d+e-1, c-e, a*b-b-e+1)],
}


def singular_form(expr):
    return str(expr).replace('**', '^')


# Dimensions and counts of maximal components over Q after localization.
COMPONENTS = {
    134217471: {4: 39, 3: 12}, 134209279: {3: 12, 2: 3},
    134217455: {3: 8, 2: 1}, 134151919: {2: 4}, 134184687: {2: 1},
    134216943: {2: 2}, 134217454: {},
}
v = symbols("u0:7")
DIVISORS = {
    134217471: [v[0]-v[3]], 134209279: [], 134217455: [v[2]-1, v[5]-1],
    134151919: [v[1]-1, v[2]-1, v[4]-1], 134184687: [v[1]-1],
    134216943: [], 134217454: [v[1]-1, v[4]-1, v[0]*v[1]*v[4]-v[2]*v[3]],
}


def program(mask):
    original, field, _ = normalize(mask)
    coeff = {t: field.to_sympy(v) for t, v in original.items()}
    assert mask.bit_count() - len(field.symbols) == 19
    zeros = [t for t in TRIPLES if coeff[t] == 0]
    assert all(sum(x == y for x, y in zip(t, u)) <= 1 for t, u in combinations(zeros, 2))
    factors = complement_factors(coeff)
    assert {f for f in DIVISORS[mask] if f in factors or -f in factors} == set(DIVISORS[mask])
    assert len(factors) == len(DIVISORS[mask])
    units = factors | set(field.symbols)
    result = ['LIB "primdec.lib";', 'option(redSB);',
              'ring r=0,(' + ','.join(map(str, field.symbols)) + '),dp;',
              'ideal tor=' + singular_form(prod(units)) + ';',
              'list cuts0,cuts1,cuts2;']
    for k in I:
        rotation = {t: coeff[t[-k:] + t[:-k]] for t in TRIPLES}
        index = 0
        for gens in failure_ideals(rotation):
            reduced = []
            for g in gens:
                if g != 0:
                    reduced.append(prod(f**n for f, n in factor_list(g)[1] if f not in units))
            if 1 in reduced:
                continue
            assert reduced
            index += 1
            args = ','.join(singular_form(f) for f in reduced)
            result.append(f'cuts{k}[{index}]=ideal({args});')
    result.append('ideal target=1;')
    if mask in STRATA:
        for stratum in STRATA[mask]:
            args = ','.join(singular_form(f) for f in stratum)
            result.append(f'target=intersect(target,ideal({args}));')
    else:
        result.append('target=0;')
    # Localize and split each small failure ideal before taking sums. Minimal
    # primes have the same zero set as the input, also after extending scalars.
    result.extend([
        'proc pieces(list lefts,list rights,ideal unitprod){',
        'list L,A,keys; int i,j,k,p,seen; ideal Q,J; string key;',
        'for(i=1;i<=size(lefts);i++){for(j=1;j<=size(rights);j++){',
        'Q=std(sat(lefts[i]+rights[j],unitprod));',
        'if(size(reduce(1,Q))!=0){A=minAssGTZ(Q);',
        'for(k=1;k<=size(A);k++){J=std(A[k]);key=string(J);seen=0;',
        'for(p=1;p<=size(keys);p++){if(key==keys[p]){seen=1;break;}}',
        'if(seen==0){keys[size(keys)+1]=key;L[size(L)+1]=J;}}',
        '}}}return(L);}',
        'list zeroes=ideal(0);',
        'cuts0=pieces(cuts0,zeroes,tor); cuts1=pieces(cuts1,zeroes,tor);',
        'cuts2=pieces(cuts2,zeroes,tor);',
        'list pairs=pieces(cuts0,cuts1,tor); list triples=pieces(pairs,cuts2,tor);',
        'int k,j,dimen,contained; int worst=-1; ideal G;',
        'for(k=1;k<=size(triples);k++){G=std(triples[k]);',
        'if(size(reduce(target,G))!=0){ERROR("stratum containment failed");}',
        'dimen=dim(G);if(dimen>worst){worst=dimen;} contained=0;',
        'for(j=1;j<=size(triples);j++){if(dim(triples[j])>dimen){',
        'if(size(reduce(triples[j],G))==0){contained=1;break;}}}',
        'if(contained==0){print("COMPONENT "+string(dimen)+" "+string(G));}}',
        'print("RESULT "+string(worst)+" "+string(size(triples)));','quit;'])
    return '\n'.join(result)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--components', action='store_true', help='print defining ideals over Q')
    args = cli.parse_args()
    binary = os.environ.get('SINGULAR', 'Singular')
    for masks in SMOOTH_MASKS.values():
        for mask in masks:
            p = subprocess.run([binary, '-q'], input=program(mask), text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
            lines = p.stdout.splitlines()
            assert not any(line.lstrip().startswith('?') for line in lines), p.stdout
            dim, count = map(int, next(line for line in lines if line.startswith('RESULT ')).split()[1:])
            expected = -1 if mask == 134217454 else mask.bit_count() - 22
            assert dim == expected, (mask, dim, expected)
            components = [line for line in lines if line.startswith('COMPONENT ')]
            counts = Counter(int(line.split()[1]) for line in components)
            assert counts == COMPONENTS[mask], (mask, counts)
            if args.components:
                for line in components:
                    print(mask, line)
            print(f'PASS mask {mask}, maximum residual dimension {dim}, {count} distinct component ideals before inclusion pruning', flush=True)


if __name__ == '__main__':
    main()
