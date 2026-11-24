"""Singular certificate for the specializations of the five three-line families.

The caller supplies the prime ideals whose entire loci already have excess
path dimension. All saturations take place in the original coefficient torus.
"""
from itertools import combinations, product

from sympy import Matrix, cancel, factor_list, prod

from verify_ou_few_zero_rotations import complement_factors, singular_form
from verify_ou_zero_supports import I, TRIPLES


# Special cell in the original tensor, number and dimension of maximal
# components where all three cyclic charts fail, in the graph order of the paper.
BOUNDARIES = {
    134217471: [((2,0,1),4,3), ((0,2,2),6,3)],
    134209279: [((0,2,2),3,2), ((2,0,0),1,2)],
    134217455: [((1,0,0),2,2)],
}


def numer(f):
    return cancel(f).as_numer_denom()[0]


def slice_equations(tensor, origin):
    labels = [[i] + [j for j in I if j != i] for i in origin]
    raw = {t: tensor[tuple(labels[i][j] for i,j in enumerate(t))] for t in TRIPLES}
    assert all(raw[i,j,0] for i,j in product(I,repeat=2) if (i,j) != (0,0))
    assert all(raw[i,j,k] for i,j in ((0,1),(0,2),(1,0),(2,0)) for k in I)
    z = {t: v/raw[t[0],t[1],0] if t[:2] != (0,0) else v for t,v in raw.items()}
    rows, factors = [], set()
    for k in (1,2):
        a,b,c,d,e,f,g,h,j = [z[i,m,k] for i,m in product(I,repeat=2)]
        K = a*(e*j-f*h)
        A = b*f*g-c*e*g+c*d*h-b*d*j
        B, C = b*(f*g-d*j), d*(c*h-b*j)
        rows.append([K,A,B,C])
        for expr in (K,B,C,B-A,C-A,K+z[0,0,0]*A):
            polynomial = numer(expr)
            assert polynomial != 0
            factors.update(v for v,_ in factor_list(polynomial)[1] if not v.is_Symbol)
    matrix = Matrix(rows)
    minors = [numer(matrix[:,cols].det()) for cols in combinations(range(4),2)]
    return minors, sorted(factors,key=str)


def boundary_program(coeff, xs, graph_eqs, denominator, origin, negative):
    alg = lambda fs: ','.join(map(singular_form,fs))
    units = prod(xs)*prod(complement_factors(coeff))
    lines = [
        'LIB "primdec.lib"; option(redSB);',
        'ring r=0,('+','.join(map(str,xs))+'),dp;',
        'ideal tor='+singular_form(units)+';',
        'ideal G=std(sat(ideal('+alg(graph_eqs)+'),ideal('+singular_form(denominator*units)+')));',
        'ideal chk,J,U; list parts,more,fac,pr,negative; int i,j,k,idx,hit;',
        'parts[1]=G;',
    ]
    for idx,prime in enumerate(negative,1):
        lines.append(f'negative[{idx}]=ideal({prime});')
    for rotation in I:
        rotated = {t: coeff[t[-rotation:]+t[:-rotation]] for t in TRIPLES}
        triple = origin[rotation:]+origin[:rotation]
        minors,factors = slice_equations(rotated,triple)
        lines.extend([
            'chk=ideal('+alg(minors)+');',
            'if(size(reduce(chk,G))!=0){ERROR("slice identities failed");}',
            'fac=list('+alg(factors)+'); more=list();',
            'for(i=1;i<=size(parts);i++){for(j=1;j<=size(fac);j++){',
            'U=std(sat(parts[i]+ideal(fac[j]),tor)); if(size(reduce(1,U))!=0){',
            'pr=minAssGTZ(U); for(k=1;k<=size(pr);k++){J=std(pr[k]);',
            'for(idx=1;idx<=size(more);idx++){if(string(more[idx])==string(J)){break;}}',
            'if(idx>size(more)){more[size(more)+1]=J;}}}}} parts=more;',
        ])
    lines.extend([
        'for(i=1;i<=size(parts);i++){hit=0;',
        'for(j=1;j<=size(negative);j++){if(size(reduce(negative[j],parts[i]))==0){hit=1;break;}}',
        'if(hit==0){ERROR("unexcluded component");} hit=0;',
        'for(j=1;j<=size(parts);j++){if(dim(parts[j])>dim(parts[i]) &&',
        'size(reduce(parts[j],parts[i]))==0){hit=1;break;}}',
        'if(hit==0){print("BOUNDARY "+string(dim(parts[i])));}}','quit;',
    ])
    return '\n'.join(lines)
