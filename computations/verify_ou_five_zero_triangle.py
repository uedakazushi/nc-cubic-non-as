#!/usr/bin/env python3
"""Exact checks for 133873391 at (u,v,w)=(1,1,1).

Component matrices and gluing over the common chart, kernels, rank witnesses,
composition generation, and multidegrees. The line-bundle proof is in the paper.
"""
from itertools import product
from sympy import Matrix, symbols, expand, cancel

MASK=133873391
I=range(3)
TRIPLES=list(product(I,repeat=3))


def main():
    c={a:(MASK>>n)&1 for n,a in enumerate(TRIPLES)}
    c={a:c[a[-1:]+a[:-1]] for a in TRIPLES}
    C=[Matrix(3,3,lambda i,j:c[i,j,k]) for k in I]
    assert C==[Matrix([[1,1,1],[1,0,1],[1,1,0]]),
               Matrix([[1,1,1],[1,1,0],[1,0,1]]),
               Matrix([[0,1,1],[1,1,1],[1,1,1]])]
    x,y,z=symbols('x y z')
    a,b,d=x+y-z,x,z
    T=x+y
    X=[Matrix([[0,1,1],[1,-z,x],[1,-x,z]]),
       Matrix([[x+z,1,1],[z,1,1],[x,1,1]]),
       Matrix([[y+z,z,y],[1,1,1],[1,1,1]])]
    K=[[Matrix([x,1,-1]),Matrix([z,1,-1]),Matrix([x+z,1,-1])],
       [Matrix([1,-x,-z]),Matrix([1,-z,-x]),Matrix([0,1,-1])],
       [Matrix([-1,1,1]),Matrix([-1,1,1]),Matrix([y-z,-y,z])]]
    # Six permutation products of component matrices, and cubic factorization.
    p=[T,T,a,x,y,z]
    assert expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]-(x+y)*(x-z)*(y-z))==0
    def products(A):
        return Matrix([A[0,0]*A[1,1]*A[2,2],A[0,0]*A[1,2]*A[2,1],
            A[0,2]*A[1,1]*A[2,0],A[0,1]*A[1,2]*A[2,0],
            A[0,2]*A[1,0]*A[2,1],A[0,1]*A[1,0]*A[2,2]])
    for l,sub in enumerate([{y:-x},{y:z},{x:z}]):
        assert products(X[l])==Matrix(p).subs(sub)
        for k in I:
            M=C[k].multiply_elementwise(X[l]);assert expand(M.det())==0
            assert M*K[l][k]==Matrix([0,0,0]) or all(expand(q)==0 for q in M*K[l][k])
            if l==0:assert M[:2,:2].det()==-1
            elif l==1 and k<2:assert M[1:,1:].det() in (1,-1)
            elif l==1:assert (M[:2,:2].det(),M.extract([0,2],[0,1]).det())==(-z,-x)
            elif k<2:assert M[1:,:2].det() in (1,-1)
            else:assert (M[:2,:2].det(),M.extract([0,1],[0,2]).det())==(-z,-y)
    # All three intersections lie on the chart a*d != 0.
    for pt in [{x:-1,y:1,z:1},{x:1,y:-1,z:1},{x:1,y:1,z:1}]:
        assert (a*d).subs(pt)!=0
    XG=Matrix([[T/d,1,1],[1,1,b/a],[1,y/a,d/a]])
    # Trivialize O(1) by a on each component and use diagonal units.
    X0=Matrix([[0,1,1],[1,1,x/(-z)],[1,(-x)/(-z),z/(-z)]])
    X1=Matrix([[(x+z)/x,1,1],[z/x,1,1],[1,1,1]])
    X2=Matrix([[(y+z)/y,z/y,1],[1,1,1],[1,1,1]])
    charts=[X0,Matrix.diag(1,1,z/x)*X1*Matrix.diag(x/z,1,1),
            Matrix.diag(y/z,1,1)*X2*Matrix.diag(1,1,z/y)]
    for Y,sub in zip(charts,[{y:-x},{y:z},{x:z}]):
        assert all(cancel(v)==0 for v in Y-XG.subs(sub))
    # On L0 the j=1,2 products are independent for i=1,2;
    # all j=0 products for i>0 are nonzero. For i=0 use another component.
    for k in I:
        for i in (1,2):
            f,g=[expand(X[0][i,j]*K[0][k][j]) for j in (1,2)]
            assert Matrix([[f.coeff(x),f.coeff(z)],[g.coeff(x),g.coeff(z)]]).det()!=0
            assert expand(X[0][i,0]*K[0][k][0])!=0
        assert all(X[0][0,j]*K[0][k][j]!=0 for j in (1,2))
        l=1 if k<2 else 2
        assert expand(X[l][0,0]*K[l][k][0])!=0
    U=Matrix([[1,0,-1],[0,0,0],[0,0,0]])
    V=Matrix([[0,1,0],[1,0,0],[1,0,0]])
    W=Matrix([[1,1,0],[1,1,0],[1,0,1]])
    A,B,C0=[Matrix([[1,1,1]])*Z for Z in (U,V,W)]
    assert A+C0-2*B==Matrix([[0,0,0]])
    assert 5*B-3*A-2*C0==Matrix([[1,1,1]])
    for P,Q in ((U,V),(V,W),(U,W)):
        for i,j in product(I,repeat=2):
            d=Q.row(j)-P.row(i)
            assert sum(d)>0 and min(d)>=-1
            # Any proper subcurve has dual degree <=1.
            assert all(sum(d[k] for k in I if mask>>k&1)>=-1 for mask in range(1,7))
    print('PASS 133873391 triangle: gluing, rank-two kernels, composition generation, multidegrees')

if __name__=='__main__':
    main()
