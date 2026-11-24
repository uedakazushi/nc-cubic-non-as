#!/usr/bin/env python3
"""Exact checks for complementary periodic objects and the last four-zero loci."""
from itertools import product, combinations
from sympy import Matrix, QQ, symbols, cancel, factor, Poly, groebner, expand
from verify_ou_zero_supports import normalize, I, TRIPLES, PAIRS
from verify_ou_five_zero_supports import recursive_dims, expected

u,v,w,q,h = symbols('u0 u1 u2 u3 h')
params=(u,v,w,q)

def tensor(mask, sub={}):
    c,K,free=normalize(mask)
    assert free==[(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    return {p:cancel(K.to_sympy(a).subs(sub,simultaneous=True)) for p,a in c.items()}

def check_complement():
    cases=[(134135535,(0,2,2),{u:1},(-1,w,v)),
           (134135535,(1,2,1),{v:1},(-1,q,u)),
           (134135535,(1,1,2),{w:1},(-1,u,q)),
           (134135535,(0,1,1),{q:1},(-1,v,w)),
           (134147823,(1,2,1),{v:1},(1,q,u)),
           (134147823,(0,1,1),{q:1},(-1,-1,w)),
           (134151407,(1,2,1),{v:1},(-1,q,u)),
           (134151918,(0,2,2),{u:1},(-1,w,v-1)),
           (134151918,(0,1,1),{q:1},(-1,v-1,w))]
    for mask,p,sub,dets in cases:
        c=tensor(mask,sub)
        rem=[(a,)+tuple(b for b in I if b!=a) for a in p]
        c={p:c[tuple(rem[j][p[j]] for j in I)] for p in TRIPLES}
        assert c[0,0,0]==0
        assert all(c[p]!=0 for p in TRIPLES if sum(a!=0 for a in p)==1)
        assert all(c[p]==1 for p in product((1,2),repeat=3))
        ds=tuple(factor(Matrix(2,2,[c[p] for p in product(*[(0,) if j==i else (1,2) for j in I])]).det()) for i in I)
        assert ds==dets,(mask,p,ds,dets)
        print('PASS complement',mask,sub,p,ds,flush=True)

def check_paths():
    for mask,sub,n,j,dim in ((134147823,{u:v-1,w:1-v,q:v},6,2,10),
                             (134151918,{u:1,v:1},8,2,19),
                             (134151918,{v:1,q:1},7,1,14)):
        a=tensor(mask,sub)
        variables=sorted(set().union(*(z.free_symbols for z in a.values())),key=str)
        K=QQ.frac_field(*variables)
        c={p:K.from_sympy(z) for p,z in a.items()}
        actual=dict(recursive_dims(c,K,0,0,n))[n][j]
        assert actual==dim and dim>expected(n,0,j)
        print('PASS fraction field',mask,sub,n,j,actual,expected(n,0,j),flush=True)

def check_line_conic():
    d,s,t=symbols('d s t')
    c={p:cancel(z.subs(v,1/(1-d))) for p,z in tensor(134147823,{w:-u,q:v}).items()}
    C=[Matrix(3,3,lambda i,j:c[i,j,k]) for k in I]
    target=[Matrix([[1,1,1],[1,0,1],[1,1,1]]),Matrix([[1,0,1],[1,1,0],[1,u,-u]]),Matrix([[1,1,0],[1,1,1],[1,1/(1-d),1/(1-d)]])]
    assert all((A-B).applyfunc(cancel)==Matrix.zeros(3) for A,B in zip(C,target))
    e=u*d*(1-d)*(u*(1-d)-d)
    H=(1-d)*s+d*t; J=s+u*t; X=-s*s-u*d*s*t+u*d*t*t
    Y=d*s*J; Z=d*t*J
    assert expand(Y*(d*X+Y)+u*d*Z*(X+Y-Z))==0
    a,b,c0=symbols('a b c0')
    F=b*(d*a+b)+u*d*c0*(a+b-c0)
    assert expand(Matrix([[F.diff(i).diff(j) for j in (a,b,c0)] for i in (a,b,c0)]).det()-2*d*d*u*(d*u+d-u))==0
    N=Matrix([[H,1,-u*d],[J/u,1,1],[X,-s,d*u*t]])
    L=Matrix([[1,t,s],[1,s,s-t],[0,1,1]])
    def unit(eq,variables):
        return list(groebner([cancel(f).as_numer_denom()[0] for f in eq]+[h*e-1],h,*variables))==[1]
    assert X.subs({s:1,t:0})==-1
    for f in (s,H,J):
        assert unit([X.subs(t,1),f.subs(t,1)],(s,u,d))
    left=Matrix.diag(1/H,u/J,-1/(s*H)); right=Matrix.diag(1,H,-s*H/(u*d))
    # On D in the t=1 chart, y/z=s. Verify gluing over the whole
    # intersection scheme by polynomial divisibility, including tangency.
    K=QQ.frac_field(u,d); modulus=Poly(X.subs(t,1),s,domain=K)
    for f in (left*N*right-L).subs(t,1):
        assert Poly(cancel(f).as_numer_denom()[0],s,domain=K).rem(modulus).is_zero
    for k in I:
        M=C[k].multiply_elementwise(N); ML=C[k].multiply_elementwise(L)
        assert cancel(M.det())==cancel(ML.det())==0
        ker=M.row(0).cross(M.row(1)).T.applyfunc(cancel)
        kerL=ML.row(0).cross(ML.row(2)).T.applyfunc(cancel)
        assert ker[0] in (1,u*d) and kerL[1] in (-1,u,1/(d-1))
        assert (M*ker).applyfunc(cancel)==Matrix.zeros(3,1)
        assert (ML*kerL).applyfunc(cancel)==Matrix.zeros(3,1)
        for i in I:
            prods=[cancel(N[i,j]*ker[j]) for j in I]
            degree=max(Poly(z,s,t).total_degree() for z in prods)
            P=Matrix([[Poly(z,s,t).coeff_monomial(s**l*t**(degree-l)) for z in prods] for l in range(degree+1)])
            for j in I: assert unit(P[:,j],(u,d))
            assert unit([P.extract(rr,cc).det() for rr in combinations(range(P.rows),2) for cc in combinations(I,2)],(u,d))
    print('PASS 134147823 line-conic, gluing, kernels and 27 products',flush=True)

if __name__=='__main__':
    check_complement()
    check_paths()
    check_line_conic()
