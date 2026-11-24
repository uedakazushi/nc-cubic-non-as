#!/usr/bin/env python3
"""Exact identities for the expanded line-conic presentation at 117358318, u=1.
The proof in cubic_zero_supports.tex gives the line bundles, their local
freeness, and cohomology; this script checks the algebraic identities.
"""
from itertools import product, combinations
from sympy import Matrix, symbols, expand, cancel, rem, Poly, QQ

I=range(3)
mask=117358318
v,s,t,r=symbols('v s t r')
D=s*s-s*t+v*t*t
C=[Matrix(3,3,lambda i,j: v if (i,j,k)==(2,2,2) else (mask>>(9*i+3*j+k))&1) for k in I]
assert C==[Matrix([[0,1,1],[1,1,1],[1,1,0]]),Matrix([[1,0,1],[1,1,0],[1,1,1]]),Matrix([[1,1,0],[1,0,1],[1,1,v]])]
# N is the conic, L is the line; both use homogeneous coordinates [s:t].
X=Matrix([[D,D,s],[t,s,1],[s,s-v*t,1]])
Y=Matrix([[0,0,1],[1,1,t],[1,1,(t-s)/v]])
# gluing over D=0 with t=1, r=s/t. This checks the length-two scheme,
# not merely the reduced points (including v=1/4).
K=QQ.frac_field(v)
def modD(f):
    n,d=cancel(f).as_numer_denom()
    p=Poly(r*r-r+v,r,domain=K)
    dn=Poly(d,r,domain=K)
    assert dn.gcd(p).degree()==0
    assert len(Poly(d,r,v).terms())==1
    return Poly(n,r,domain=K).rem(p).is_zero
R=Matrix.diag(1/r,1,1/r)
S=Matrix.diag(1,1/r,1)
assert all(modD(a) for a in R*X.subs({s:r,t:1})*S-Y.subs({s:r,t:1}))
# Six products give the normalization map of the determinant cubic.
ps=[(0,1,2),(0,2,1),(2,1,0),(1,2,0),(2,0,1),(1,0,2)]
from sympy import prod
assert [expand(prod(X[i,p[i]] for i in I)) for p in ps]==[expand(a) for a in [s*D,(s-v*t)*D,s**3,s*D,s*t*(s-v*t),t*D]]
# Kernel generators on each component, with their specified twists.
Kn=[Matrix([-(s-v*t),s,-D]),Matrix([s,-t,-D]),Matrix([1,-1,-t])]
Kl=[Matrix([1,-1,0]),Matrix([1,-1,0]),Matrix([-t,s,1])]
# Last component: second row b0+t*b2=0, third b0+b1+(t-s)*b2=0.
for k in I:
    M=C[k].multiply_elementwise(X);N=C[k].multiply_elementwise(Y)
    assert expand(M.det())==0 and cancel(N.det())==0
    assert all(expand(a)==0 for a in M*Kn[k])
    assert all(cancel(a)==0 for a in N*Kl[k])
# All compositions are nonzero and at least two independent on N.
for i,k in product(I,repeat=2):
    sections=[expand(X[i,j]*Kn[k][j]) for j in I]
    assert all(any(len(Poly(b,v).terms())==1 and b!=0 for b in Poly(a,s,t).coeffs()) for a in sections)
    assert expand(sum(C[k][i,j]*sections[j] for j in I))==0
    # Certify rank >=2 uniformly on v!=0 by a minor equal to a unit monomial.
    deg=Poly(sections[0],s,t).total_degree()
    P=Matrix([[Poly(a,s,t).coeff_monomial(s**(deg-d)*t**d) for a in sections] for d in range(deg+1)])
    minors=[expand(P.extract(rs,cs).det()) for rs in combinations(range(deg+1),2) for cs in combinations(I,2)]
    assert any(m!=0 and len(Poly(m,v).terms())==1 for m in minors),(i,k,minors)
# Multidegrees and the cohomological bounds.
U=Matrix([[-1,1],[0,0],[0,0]])
V=Matrix([[1,0],[1,0],[0,1]])
W=Matrix([[2,0],[2,0],[1,1]])
for block,d in [(U,0),(V,1),(W,2)]:
    assert all(sum(block.row(i))==d for i in I)
    assert all(abs(a)<=1 for i,j in product(I,repeat=2) for a in block.row(i)-block.row(j))
for X0,Y0 in [(U,V),(V,W),(U,W)]:
    for i,j in product(I,repeat=2):
        assert min(Y0.row(j)-X0.row(i))>=-1
A,B,C0=[Matrix([[1,1,1]])*T for T in [U,V,W]]
assert A+C0-2*B==Matrix([[0,0]])
assert 5*B-3*A-2*C0==Matrix([[3,0]])
print('PASS: six-zero residual, scheme gluing, kernel identities, all nine composition ranks, multidegrees')
