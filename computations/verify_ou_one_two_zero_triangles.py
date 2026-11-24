#!/usr/bin/env python3
"""Verify the relative three-line model on each of the five residual graph closures."""
from itertools import combinations
from sympy import Matrix, S, cancel, expand, factor, Poly, prod, symbols
from verify_ou_one_two_zero_components import FAMILIES, a, b, c, d, graph, run_singular
from verify_ou_zero_supports import I, normalize
from verify_ou_few_zero_rotations import singular_form
from verify_ou_seven_zero_curves import determinant_matrix

x,p,q,s,t=symbols('x p q s t')
def relabel(c,entry):
    rs=[entry[0]]+[i for i in I if i!=entry[0]]
    cs=[entry[1]]+[j for j in I if j!=entry[1]]
    z={t:c[rs[t[0]],cs[t[1]],t[2]] for t in c}
    return z


def model(c, entry):
    z = relabel(c, entry)
    L=determinant_matrix(z)
    # Normalize the proportional pair and the tangent hyperplane to the quadric.
    basis=L.nullspace()
    N=Matrix.hstack(*basis)
    assert N.shape==(6,3)
    rho=cancel(next(N[1,j]/N[0,j] for j in I if N[0,j]!=0))
    assert all(cancel(N[1,j]-rho*N[0,j])==0 for j in I)
    v=N[2:,:].T.nullspace();assert len(v)==1
    A,B,C,D=[factor(f) for f in v[0]]
    ss=[rho,1,A,-B,-C,D]
    assert all(ss) and cancel(rho*B*C-A*D)==0
    mult=Matrix([[cancel(A*D/rho),1,1],[1,1/A,-1/B],[1,-1/C,1/D]])
    ys=[(N[j,:]*ss[j]).applyfunc(cancel) for j in range(6)]
    h=(ys[0]*Matrix.vstack(ys[2],ys[3]-ys[2],ys[4]-ys[2]).inv()).applyfunc(factor)
    return [Matrix(3,3,lambda i,j:z[i,j,k]*mult[i,j]) for k in I],list(h)
def sample(c,entry):
    rows,H=model(c,entry)
    T=H[0]*x+H[1]*p+H[2]*q
    w=x+p+q
    assert H[0]!=0 and H[1]!=0 and H[2]!=0 # triangle independent and x nonzero at corners
    corners=[(0,-H[0],H[2]),(-H[0],0,H[1]),(0,0,1)] # (p,q,x)
    assert all(xx and xx+pp+qq for pp,qq,xx in corners)
    mats=[Matrix([[0,1,1],[1,x,x+p],[1,x+q,w]]).subs({p:s,q:t,x:-(H[1]*s+H[2]*t)/H[0]}),
          Matrix([[T,1,1],[w,1,1],[x,1,1]]).subs({p:0,x:s,q:t}),
          Matrix([[T,w,x],[1,1,1],[1,1,1]]).subs({q:0,x:s,p:t})]
    ks={};vs=((0,1,0),(1,0,0),(1,0,0));degs=[]
    for k in I:
        dd=[]
        for l,X in enumerate(mats):
            M=rows[k].multiply_elementwise(X)
            assert expand(M.det())==0
            cross=[M.row(i).cross(M.row(j)) for i,j in ((0,1),(1,2),(2,0))]
            vals=[expand(f) for v in cross for f in v if expand(f)!=0]
            g=Poly(vals[0],s,t)
            for f in vals[1:]:g=g.gcd(Poly(f,s,t))
            assert g.total_degree()==0
            v=next(v for v in cross if any(expand(f)!=0 for f in v))
            g=Poly(0,s,t)
            for f in v:g=g.gcd(Poly(f,s,t))
            v=[cancel(f/g.as_expr()) for f in v]
            dv = {sum(mon)+vs[j][l] for j,f in enumerate(v) if f
                  for mon in Poly(f,s,t).monoms()}
            assert len(dv)==1
            ks[k,l]=v;dd.append(dv.pop())
        degs.append(dd)
        for i in I:
            out=[]
            for l in I:
                polys=[Poly(ks[k,l][j]*mats[l][i,j],s,t) for j in I]
                exps=set().union(*(poly.monoms() for poly in polys))
                out.extend([[v.coeff_monomial(mon) for v in polys] for mon in exps])
            M=Matrix(out)
            assert M.rank()==2 and all(any(M[:,j]) for j in I)
    assert sorted(degs) == [[1,0,1],[1,1,0],[1,1,0]]
    for i,j in combinations(I,2):
        if degs[i]!=degs[j]:
            continue
        first, second = rows[i], rows[j]
        if any(bool(u)!=bool(v) for u,v in zip(first,second)):
            continue
        rectangles = [((i1,j1),(i2,j2),(i1,j2),(i2,j1))
                      for i1,i2 in combinations(I,2) for j1,j2 in combinations(I,2)]
        assert any(first[u]*first[v]*second[w]*second[z] != second[u]*second[v]*first[w]*first[z]
                   for u,v,w,z in rectangles)
    return H,degs


def generic_identities(original, entry, values, xs, witness):
    # Prove the plane projection and tangent-quadric identities in the
    # localized graph ideal, before computing bundle degrees at a single fibre.
    mat = determinant_matrix(relabel(original,entry))
    at = dict(zip(xs, [S(v).subs(witness) for v in values]))
    last = mat[:,2:]
    cross = next(last[:,i].cross(last[:,j]) for i,j in combinations(range(4),2)
                 if any(last[:,i].cross(last[:,j]).subs(at)))
    h = mat[:,0].cross(mat[:,1])
    A,B,C,D = [h.dot(last[:,j]) for j in range(4)]
    numer, denom = -cross.dot(mat[:,0]), cross.dot(mat[:,1])
    assert all(f.subs(at)!=0 for f in (numer,denom,A,B,C,D))
    identities = [last[:,list(cols)].det() for cols in combinations(range(4),3)]
    identities.append(numer*B*C-denom*A*D)
    equations, den, _ = graph(values,xs)
    code = ['LIB "primdec.lib";','ring r=0,('+','.join(map(str,xs))+'),dp;',
            'ideal graph='+','.join(map(singular_form,equations))+';',
            'ideal G=std(sat(graph,ideal('+singular_form(den*prod(xs))+')));',
            'ideal tests='+','.join(singular_form(expand(f)) for f in identities)+';',
            'if(size(reduce(tests,G))!=0){ERROR("three-line identity failed");}',
            'print("PASSED IDENTITIES");quit;']
    assert 'PASSED IDENTITIES' in run_singular('\n'.join(code))
# Original common arrow of the proportional even/odd pair. The families are
# ordered as the five unexcluded graphs in FAMILIES, in the paper's table.
ENTRIES = {134217471: [(2,0),(0,2)], 134209279: [(0,2),(2,0)], 134217455: [(1,0)]}


def main():
    witness = {a:3,b:5,c:11,d:19}
    for mask, families in FAMILIES.items():
        original, field, free = normalize(mask)
        coefficients = {u:field.to_sympy(v) for u,v in original.items()}
        remaining = [values for values,path,_ in families if path is None]
        assert len(remaining) == len(ENTRIES[mask])
        for entry,values in zip(ENTRIES[mask],remaining):
            generic_identities(coefficients,entry,values,field.symbols,witness)
            parameters = dict(zip(free, values))
            fibre = {u:S(parameters.get(u,v)).subs(witness) for u,v in coefficients.items()}
            plane,degrees = sample(fibre,entry)
            print('PASS',mask,'arrow',entry,'plane',plane,'multidegrees',degrees,flush=True)


if __name__ == '__main__':
    main()
