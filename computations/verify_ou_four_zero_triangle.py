#!/usr/bin/env python3
"""Exact checks for the residual 133889775 family on three lines."""
from itertools import combinations, product
from sympy import Matrix, symbols, expand, cancel, groebner, Poly
from verify_ou_zero_supports import normalize, TRIPLES, I, PAIRS
from verify_ou_seven_zero_curves import determinant_matrix


def main():
    u,v,w,q=symbols('u0 u1 u2 u3')
    x,y,z,h=symbols('x y z h')
    c,K,free=normalize(133889775)
    assert free==[(1,1,2),(2,1,2),(2,2,1),(2,2,2)]
    c={a:K.to_sympy(c[a[-1:]+a[:-1]]).subs({v:1,q:1,w:1-u},simultaneous=True) for a in TRIPLES}
    w=1-u
    C=[Matrix(3,3,lambda i,j:c[i,j,k]) for k in I]
    assert C==[Matrix([[1,1,1],[1,0,1],[1,1,0]]),Matrix([[1,1,1],[1,1,u],[1,0,1]]),Matrix([[0,1,1],[1,1,1],[1,w,1]])]
    a=x+w*y-z;T=x+y
    p=Matrix([w*T,T,a,x,y,z])
    assert all(expand(e)==0 for e in determinant_matrix(c)*p)
    assert expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]-T*(x-z)*(w*y-z))==0
    subs=[{y:-x},{z:w*y},{x:z}]
    coords=[(x,z),(x,y),(y,z)]
    X=[Matrix([[0,1,1],[1,u*x-z,x],[1,-x,z]]),
       Matrix([[w*(x+y),1,1],[w*y,1,1],[x,1/w,1]]),
       Matrix([[z+y,z,y],[1,w,1],[1,1,1]])]
    Y=[[Matrix([x,1,-1]),Matrix([z,1,-1]),Matrix([w*x+z,1,-1])],
       [Matrix([-1/w,x,y]),Matrix([-1,w*y-u*x,x]),Matrix([0,-1,1])],
       [Matrix([-1,1,1]),Matrix([-1,1,1]),Matrix([w*y-z,-y,z])]]
    def unit(eq,vs):
        assert list(groebner([cancel(e).as_numer_denom()[0] for e in eq]+[h*u*w-1],h,*vs,u))==[1]
    def perms(A):
        return Matrix([A[0,0]*A[1,1]*A[2,2],A[0,0]*A[1,2]*A[2,1],A[0,2]*A[1,1]*A[2,0],A[0,1]*A[1,2]*A[2,0],A[0,2]*A[1,0]*A[2,1],A[0,1]*A[1,0]*A[2,2]])
    for l,(A,B,sub,vs) in enumerate(zip(X,Y,subs,coords)):
        assert all(cancel(e)==0 for e in perms(A)-p.subs(sub))
        for k in I:
            M=C[k].multiply_elementwise(A)
            assert cancel(M.det())==0 and all(cancel(e)==0 for e in M*B[k])
            for j in (0,1):
                ms=[M.extract(rr,cc).det().subs(vs[j],1) for rr,cc in product(combinations(I,2),repeat=2)]
                unit(ms,(vs[1-j],))
                unit([e.subs(vs[j],1) for e in B[k]],(vs[1-j],))
    # All intersections are in D_+(a*z), including the concurrent case.
    for P in [{x:1,y:-1,z:-w},{x:1,y:-1,z:1},{x:w,y:1,z:w}]:
        unit([(a*z).subs(P,simultaneous=True)],())
    Xg=Matrix([[w*T/z,1,1],[1,1,x/a],[1,y/a,z/a]])
    X0=X[0].copy()
    for i,j in product((1,2),(1,2)): X0[i,j]/=u*x-z
    X1=X[1].copy()
    for i,j in [(0,0),(1,0),(2,0)]: X1[i,j]/=x
    X2=X[2].copy()
    for j in I: X2[0,j]/=w*y
    charts=[X0,Matrix.diag(1,1,w*y/x)*X1*Matrix.diag(x/(w*y),1,1),
            Matrix.diag(w*w*y/z,1,1)*X2*Matrix.diag(1,1/w,z/(w*y))]
    for A,sub in zip(charts,subs): assert all(cancel(e)==0 for e in A-Xg.subs(sub))
    for k in I:
        for i in (1,2):
            P=[expand(X[0][i,j]*Y[0][k][j]) for j in (1,2)]
            unit([Matrix([[f.coeff(x),f.coeff(z)] for f in P]).det()],())
            unit(Poly(X[0][i,0]*Y[0][k][0],x,z).coeffs(),())
        for j in (1,2): assert X[0][0,j]*Y[0][k][j]!=0
        l=1 if k<2 else 2
        unit(Poly(X[l][0,0]*Y[l][k][0],*coords[l]).coeffs(),())
    U=Matrix([[1,0,-1],[0,0,0],[0,0,0]])
    V=Matrix([[0,1,0],[1,0,0],[1,0,0]])
    W=Matrix([[1,1,0],[1,1,0],[1,0,1]])
    A,B,Csum=[Matrix([[1,1,1]])*M for M in (U,V,W)]
    assert A+Csum-2*B==Matrix([[0,0,0]]) and 5*B-3*A-2*Csum==Matrix([[1,1,1]])
    for M,N in ((U,V),(V,W),(U,W),(U,U),(V,V),(W,W)):
        for i,j in PAIRS:
            d=N.row(j)-M.row(i)
            assert all(sum(d[k] for k in I if mask>>k&1)>=-1 for mask in range(1,7))
    print('PASS 133889775: all parameters u(1-u) nonzero, including concurrent lines')

if __name__=='__main__':
    main()
