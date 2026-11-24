#!/usr/bin/env python3
"""Exact checks for four expanded line-conic presentations of five-zero tensors.

The paper proves gluing, the kernel bundle and cohomology statements.
All ideal checks are over Q, localized at the stated nonzero parameters.
"""
from itertools import combinations, product
from sympy import Matrix, symbols, expand, cancel, Poly, groebner, invert
from verify_ou_zero_supports import normalize, TRIPLES, I, PAIRS

a,b,s,t,h = symbols('a b s t h')


def unit_ideal(equations, variables, nonzero):
    eq = [cancel(f).as_numer_denom()[0] for f in equations]
    assert list(groebner(eq + [h*nonzero-1], h, *variables)) == [1], (equations, variables, nonzero)


def check(mask, subs, perm, C, S,T,U,V,R,D, kernels, extra):
    c,K,free=normalize(mask)
    u,v,w=symbols('u0 u1 u2')
    c={q:K.to_sympy(z).subs(subs,simultaneous=True) for q,z in c.items()}
    p,q,r=perm
    assert C==[Matrix(3,3,lambda i,j:c[p[i],q[j],r[k]]) for k in I]
    X=Matrix([[D,D,R],[S,T,1],[U,V,1]])
    nonzero=a*b*extra
    assert expand(D.subs({s:1,t:0}))!=0
    unit_ideal([D.subs({s:1,t:0})],(a,b),nonzero)
    assert expand(T*U-S*V) != 0
    # Every entry needed to glue is a unit on the length-two intersection.
    for l in (R,S,T,U,V):
        unit_ideal([D.subs(t,1),l.subs(t,1)],(s,a,b),nonzero)
    for k in I:
        M=C[k].multiply_elementwise(X)
        assert expand(M.det())==0
        assert all(expand(f)==0 for f in M*kernels[k])
        # The kernel vector has no common zero (correct component degrees),
        # and rank is two everywhere on N, including both affine charts.
        for coord,other in ((s,t),(t,s)):
            minors=[M.extract(list(i),list(j)).det().subs(coord,1)
                    for i,j in product(combinations(I,2),repeat=2)]
            unit_ideal(minors,(other,a,b),nonzero)
            unit_ideal([f.subs(coord,1) for f in kernels[k]],(other,a,b),nonzero)
    # Inverse linear forms in the quotient by D(s,1).  Invert S and U
    # using a 2x2 linear system; verify denominators are units on the locus.
    f=Poly(D.subs(t,1),s)
    from sympy import QQ
    K=QQ.frac_field(a,b)
    f=Poly(f,s,domain=K)
    L=[]
    for l in (S,U):
        g=Poly(l.subs(t,1),s,domain=K)
        v=invert(g,f).as_expr()
        unit_ideal([cancel(v).as_numer_denom()[1]],(a,b),nonzero)
        assert (Poly(v*l.subs(t,1)-1,s,domain=K)%f).is_zero
        L.append(expand(v).subs(s,s/t)*t)
    Xl=Matrix([[0,0,1],[1,1,L[0]],[1,1,L[1]]]).applyfunc(cancel)
    for k in I:
        M=C[k].multiply_elementwise(Xl)
        assert cancel(M.det())==0
        if k<2:
            assert M[0,2]==1 and M[1,0]==1
            assert M*Matrix([1,-1,0])==Matrix.zeros(3,1)
        else:
            assert M[0,:]==Matrix([[0,0,0]])
            unit_ideal([M[1:,:2].det()],(a,b),nonzero)
            W=M.row(1).cross(M.row(2)).T
            assert all(cancel(g)==0 for g in M*W)
    # All three products are nonzero for every parameter; at least two
    # generate, checked by coefficient-vector minors on N.
    for i,k in PAIRS:
        prods=[expand(X[i,j]*kernels[k][j]) for j in I]
        degree=max(Poly(f,s,t).total_degree() for f in prods)
        P=Matrix([[Poly(f,s,t).coeff_monomial(s**l*t**(degree-l))
                   for f in prods] for l in range(degree+1)])
        for j in I:
            unit_ideal(P[:,j],(a,b),nonzero)
        unit_ideal([P.extract(list(l),list(j)).det()
                    for l in combinations(range(degree+1),2)
                    for j in combinations(I,2)],(a,b),nonzero)
    # Distinguish equal-multidegree U and V using A10,A20 and A10,A11.
    for f,g in ((S,U),(S,T)):
        unit_ideal([Matrix([[f.coeff(s),f.coeff(t)],
                            [g.coeff(s),g.coeff(t)]]).det()],(a,b),nonzero)
    assert {p for p in PAIRS if C[0][p]==0}!={p for p in PAIRS if C[1][p]==0}
    print('PASS expanded curve',mask,'extra',extra,flush=True)


def main():
    u,v,w=symbols('u0 u1 u2')
    ident=(0,1,2);swap=(0,2,1)
    # 117374702: u=1, a=v, b=w, a != 1.
    S,T,U,V,R=t,s,s,((b+1)*s-b*t)/a,s
    D=s*s-t*V
    check(117374702,{u:1,v:a,w:b},(ident,ident,ident),
          [Matrix([[0,1,1],[1,1,1],[1,1,0]]),
           Matrix([[1,0,1],[1,1,0],[1,1,1]]),
           Matrix([[1,1,0],[1,1,1],[1,a,b]])],
          S,T,U,V,R,D,[Matrix([-V,s,-D]),Matrix([s,-t,-D]),Matrix([1,-1,s-t])],a-1)
    # 132054766: w=1, a=u, b=v, b != -1. Swap j,k labels 1,2.
    S,T,U,V,R=a*s-(b+1)*t,s,a*s-t,t,a*s-t
    D=T*U-S*V
    check(132054766,{u:a,v:b,w:1},(ident,swap,swap),
          [Matrix([[0,1,1],[1,1,1],[1,1,0]]),
           Matrix([[1,0,1],[1,1,1],[1,1,a]]),
           Matrix([[1,1,0],[1,0,1],[1,b,1]])],
          S,T,U,V,R,D,[Matrix([-t,U,-D]),Matrix([U,U-a*S,-D]),Matrix([1,-1,-S])],b+1)
    # 134135534: u=1, a=v, b=w, b != 1.
    S,T,U,V,R=t,s,a*s+t,a*s+(1-b)*t,s
    D=expand((T*U-S*V)/a)
    check(134135534,{u:1,v:a,w:b},(ident,ident,ident),
          [Matrix([[0,1,1],[1,1,1],[1,1,1]]),
           Matrix([[1,0,1],[1,1,0],[1,1,a]]),
           Matrix([[1,1,0],[1,0,1],[1,1,b]])],
          S,T,U,V,R,D,[Matrix([expand((D-s*s)/t).cancel(),s,-D]),Matrix([s,-t,-D]),Matrix([1,-1,-t])],b-1)
    # 134147822: w=1, a=u, b=v. Swap labels 1,2 on all axes.
    S,T,U,V,R=(b+1)*s-a*t,s,s-a*t,t,s-a*t
    D=T*U-S*V
    check(134147822,{u:a,v:b,w:1},(swap,swap,swap),
          [Matrix([[0,1,1],[1,1,1],[1,1,0]]),
           Matrix([[1,0,1],[1,1,a],[1,1,1]]),
           Matrix([[1,1,0],[1,b,1],[1,0,1]])],
          S,T,U,V,R,D,[Matrix([-t,U,-D]),Matrix([U,a*U-S,-D]),Matrix([1,-1,-U])],1)
    # Swap j,k labels 1,2 and rescale to transform (u,v,w) -> (w,1/v,u).
    c,K,free=normalize(134135534)
    c={q:K.to_sympy(z) for q,z in c.items()}
    rescale={2:1,4:-1,7:-1,12:1,13:1,16:-1,18:-1}
    for i,j,k in TRIPLES:
        e=sum(rescale.get(l,0) for l in (3*i+j,9+3*j+k,18+3*k+i))
        assert cancel(c[i,swap[j],swap[k]]*v**e-c[i,j,k].subs({u:w,v:1/v,w:u},simultaneous=True))==0
    U=Matrix([[-1,1],[0,0],[0,0]])
    V=Matrix([[1,0],[1,0],[0,1]])
    W=Matrix([[2,0],[2,0],[1,1]])
    A,B,C=[Matrix([[1,1,1]])*M for M in (U,V,W)]
    assert A+C-2*B==Matrix([[0,0]]) and 5*B-3*A-2*C==Matrix([[3,0]])
    for P,Q in ((U,V),(V,W),(U,W)):
        for i,j in PAIRS:
            assert min(Q.row(j)-P.row(i))>=-1 and sum(Q.row(j)-P.row(i))>0
    print('PASS gauge symmetry and common multidegrees',flush=True)

if __name__=='__main__':
    main()
