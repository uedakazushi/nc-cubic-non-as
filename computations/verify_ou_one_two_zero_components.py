#!/usr/bin/env python3
"""Exact path exclusions and symmetry orbits of the residual one/two-zero loci."""
from itertools import permutations, product
import os
import subprocess

from sympy import Matrix, QQ, S, cancel, prod, symbols

from verify_ou_few_zero_rotations import program, singular_form
from verify_ou_five_zero_supports import expected
from verify_ou_four_zero_supports import laurent_rank, rows_for
from verify_ou_zero_supports import I, incidence, normalize, path_dim, support


a, b, c, d = PARAMETERS = symbols("a b c d")
h = a*(1-c)/(a-c)
q = b*c+b*d-b-c*d
t = a*b+a*d-a-b*d
# Each tuple consists of the rational graph, an optional path witness (length,
# rotation, endpoints, dimension), and the number of different component ideals.
FAMILIES = {
    134217471: [
        ((a,b,c,1,1,c/a,1/b), (4,0,2,1,6), 12),
        ((a,b,c,d,h,d*(1-h)/(d-h),(a-h)/(a-b)), (5,0,2,1,8), 12),
        ((a,b,c,a*d,b*d,c*d,d), (7,1,1,1,13), 12),
        ((1,1,1,a,b,c,d), (7,1,0,0,13), 6),
        ((a,b,c,d,h,d*(1-b)/(a-b),a*(a-1)/((a-b)*(a-c))), None, 8),
        ((a,b,c,d*(a*d-a+b*c-c*d)/q,d,c*d*(b-1)/q,
          d*(a-1)/(a*d+b*c-b-c*d)), None, 1),
    ],
    134209279: [
        ((1,1,a,1,1,b), (4,1,0,0,6), 3),
        ((a,b,c,(b-1)*c/b,b+c-b*c,((b-1)*c/b-1)/(a-1)), (5,0,2,0,8), 3),
        ((1,a,1,1,b,c), (7,0,0,1,13), 6),
        ((a,b,b*d*(a-d)/t,d,b*d*(a-1)/t,-d/(a*b-a-b*d)), None, 2),
        ((a,b,c,c*(b-1)/b,c*(a-1)/a,-c/(a*b)), None, 1),
    ],
    134217455: [
        ((1,1,a,1,1,b), (4,1,0,0,6), 1),
        ((a,1-a,b,c,1-a*b,(1-c)/(1-a)), (5,0,1,0,8), 4),
        ((a,b,c,a*c,b*c,c), (7,1,1,1,13), 1),
        ((a,b,c,a*c-c+1,b*c-c+1,c), (5,2,0,0,8), 1),
        ((a,b,c,1-b*c,1-a*c,c), None, 2),
    ],
}


def run_singular(source):
    binary = os.environ.get('SINGULAR', 'Singular')
    p = subprocess.run([binary, '-q'], input=source, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, check=True)
    lines = p.stdout.splitlines()
    assert not any(line.lstrip().startswith('?') for line in lines), p.stdout
    return lines


def symmetry_maps(mask, coeff, free):
    sup = support(mask)
    matrix = incidence(sup)
    pivots = list(matrix.T.rref()[1])
    pivot_rows = matrix[pivots,:]
    cols = list(pivot_rows.rref()[1])
    exponents = matrix[[sup.index(t) for t in free],cols]*pivot_rows[:,cols].inv()
    assert all(e.is_integer for e in exponents)
    maps = set()
    perms = list(permutations(I))
    for axes in perms:
        for labels in product(perms, repeat=3):
            def change(t):
                return tuple(labels[i][t[axes[i]]] for i in I)
            if any(change(t) not in coeff or coeff[change(t)] == 0 for t in sup):
                continue
            piv = [coeff[change(sup[i])] for i in pivots]
            maps.add(tuple(coeff[change(t)]/prod(v**e for v,e in zip(piv,exponents.row(j)))
                           for j,t in enumerate(free)))
    return maps


def graph(values, xs):
    values = tuple(map(S, values))
    pars = sorted(set().union(*(v.free_symbols for v in values)),key=str)
    at = {p: xs[values.index(p)] for p in pars}
    equations, denominators = [], set()
    for x,v in zip(xs, values):
        num, den = cancel(x-v.xreplace(at)).as_numer_denom()
        if num:
            equations.append(num)
            denominators.add(den)
    return equations, prod(denominators), len(pars)


def path_check(coeff, free, values, witness):
    n, rot, i, j, dim = witness
    values = tuple(map(S, values))
    pars = sorted(set().union(*(v.free_symbols for v in values)),key=str)
    k = QQ.frac_field(*pars)
    at = dict(zip(free,values))
    tensor = {t: k.from_sympy(at.get(t,z)) for t,z in coeff.items()}
    if n < 7:
        value = path_dim(tensor,n,rot,i,j,k)
    else:
        try:
            from flint import fmpz_mpoly_ctx
            from verify_ou_collinear_zeros import flint_rank
        except ImportError:
            flint_rank = None
        ncols, rows = rows_for(tensor,n,rot,i,j,k)
        value = ncols-laurent_rank(rows,k,finish=flint_rank)
    assert value == dim > expected(n,i,j)


def main():
    for mask, families in FAMILIES.items():
        orig, field, free = normalize(mask)
        coeff = {t: field.to_sympy(v) for t,v in orig.items()}
        xs = field.symbols
        maps = symmetry_maps(mask, coeff, free)
        primes = [ln.split(' ',2)[1:] for ln in run_singular(program(mask))
                  if ln.startswith('COMPONENT ')]
        header = ['ring r=0,('+','.join(map(str,xs))+'),dp;', 'list G; ideal J; poly h;']
        for i, (_,gs) in enumerate(primes):
            header.append(f'G[{i+1}]=std(ideal({gs}));')
        seen = set()
        for values, witness, size in families:
            if witness:
                path_check(coeff,free,values,witness)
            gs, den, dim = graph(values,xs)
            code = header.copy()
            for image in maps:
                sub = dict(zip(xs,image))
                nums = [cancel(f.xreplace(sub)).as_numer_denom()[0] for f in gs]
                h = cancel(den.xreplace(sub)).as_numer_denom()[0]
                code.append('J=ideal('+','.join(map(singular_form,nums))+');h='+singular_form(h)+';')
                for i,(dimension,_) in enumerate(primes):
                    if int(dimension) == dim and i not in seen:
                        code.append(f'if(size(reduce(J,G[{i+1}]))==0 && reduce(h,G[{i+1}])!=0)'
                                    f'{{print("HIT {i}");}}')
            matches = {int(line[4:]) for line in run_singular('\n'.join(code+['quit;']))
                       if line.startswith('HIT ')}
            assert len(matches) == size, (mask,values,matches,size)
            seen.update(matches)
            print('PASS',mask,'path '+str(witness) if witness else 'remaining graph',size,flush=True)
        assert len(seen) == len(primes)


if __name__ == '__main__':
    main()
