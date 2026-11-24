#!/usr/bin/env python3
"""Symbolic checks for the residual three-line nine-bundle construction.

Only polynomial identities and two examples are checked here. The general
line-bundle gluing, constant-rank and cohomology arguments are in the paper.
"""
from itertools import permutations, product
from sympy import Matrix, symbols, expand, factor, cancel, prod


def check_identities():
    a,b,c,d,l,m=symbols('a b c d l m')
    t,x,p,q=symbols('t x p q')
    y,z,w=x+p,x+q,x+p+q
    D=a*d-b*c
    S=a+d-b-c
    tau=1+l-m
    F1=(m-1)*D+l*S
    A=lambda v:1+v
    B=lambda v:1+l+m*v
    Aprod=A(a)*A(d)-A(b)*A(c)
    Bprod=B(a)*B(d)-B(b)*B(c)
    assert expand(Aprod-(D+S))==0
    assert expand(Bprod-m*(D+S)-m*F1)==0
    z0,z1=symbols('z0 z1')
    def ell(f):
        return f(a)*f(d)*z0-f(b)*f(c)*z1-f(a)*x+f(b)*y+f(c)*z-f(d)*w
    ell0,ell1,ell2=ell(lambda _:1),ell(A),ell(B)
    r=expand(ell2-m*ell1-tau*ell0)
    P=B(a)*B(d)-m*A(a)*A(d)-tau
    assert expand(r-P*(z0-z1)-m*F1*z1)==0
    H=(D+S)*t-S*x+(b-d)*p+(c-d)*q
    assert expand(ell1.subs({z0:t,z1:t})-H)==0
    assert expand((z0*y*z-z1*x*w).subs({z0:t,z1:t})-t*p*q)==0
    # Component matrices; projective linear entries are sections of O(1).
    X0=Matrix([[0,1,1],[1,x,y],[1,z,w]])
    X1=Matrix([[t,1,1],[w,1,1],[x,1,1]])
    X2=Matrix([[t,w,x],[1,1,1],[1,1,1]])
    e,f,g,h=symbols('e f g h')
    C=Matrix([[1,1,1],[1,e,f],[1,g,h]])
    E=(e*h-f*g)*t-e*x+f*y+g*z-h*w
    for X,var in [(X0,t),(X1,p),(X2,q)]:
        assert expand((C.multiply_elementwise(X).det()-E).subs(var,0))==0
    # At every singular point x*w!=0; trivialize O(1) by x.
    # The common local matrix on this chart is defined over O_E,
    # not merely on its normalization (important for concurrent lines).
    Xloc=Matrix([[t/w,1,1],[1,1,y/x],[1,z/x,w/x]])
    X0loc=Matrix([[0,1,1],[1,1,y/x],[1,z/x,w/x]])
    X1loc=Matrix([[t/x,1,1],[w/x,1,1],[1,1,1]])
    X2loc=Matrix([[t/x,w/x,1],[1,1,1],[1,1,1]])
    I=Matrix.eye(3)
    R1=Matrix.diag(1,1,w/x); S1=Matrix.diag(x/w,1,1)
    R2=Matrix.diag(x/w,1,1); S2=Matrix.diag(1,1,w/x)
    for X,R,Sg,var in [(X0loc,I,I,t),(X1loc,R1,S1,p),(X2loc,R2,S2,q)]:
        assert all(cancel(v.subs(var,0))==0 for v in R*X*Sg-Xloc)
    # Explicit kernel along L0 and rank-two witnesses on the other lines.
    Y0=C.multiply_elementwise(X0)
    k0=Matrix([f*y-e*x,1,-1])
    assert all(expand(v)==0 for v in Y0*k0-Matrix([0,0,E.subs(t,0)]))
    assert Y0[:2,:2].det()==-1
    assert C.multiply_elementwise(X1)[:2,1:].det()==f-e
    assert C.multiply_elementwise(X2)[1:,:2].det()==g-e
    # Multidegrees (columns are L0,L1,L2).
    U=Matrix([[1,0,-1],[0,0,0],[0,0,0]])
    V=Matrix([[0,1,0],[1,0,0],[1,0,0]])
    W=Matrix([[1,0,1],[1,1,0],[1,1,0]])
    assert [sum(U.row(i)) for i in range(3)]==[0]*3
    assert [sum(V.row(i)) for i in range(3)]==[1]*3
    assert [sum(W.row(i)) for i in range(3)]==[2]*3
    for X,Y in [(U,V),(V,W)]:
        for i,j in product(range(3),repeat=2):
            diff=list(Y.row(j)-X.row(i))
            assert sum(diff)==1 and all(-1<=v<=1 for v in diff)
    for i,j in product(range(3),repeat=2):
        assert min(W.row(j)-U.row(i))>=0
    A0,B0,C0=[Matrix([[1,1,1]])*Z for Z in (U,V,W)]
    assert A0+C0-2*B0==Matrix([[0,0,0]])
    assert 5*B0-3*A0-2*C0==Matrix([[1,1,1]])
    print('PASS: residual plane, three-line factorization, gluing, kernels, multidegrees',flush=True)


def check_examples():
    # Both geometries occur among coefficient tensors satisfying the indicated
    # low-degree conditions. No claim about all iterates of inverse Serre.
    I=range(3)
    P=sorted(permutations(I),key=lambda p:sum(p[i]>p[j] for i in I for j in range(i+1,3))%2)
    for a,b,c,d,l in [(2,3,5,8,-6),(2,3,5,6,2)]:
        from sympy import Rational
        D=a*d-b*c; S=a+d-b-c; m=1-Rational(l*S,D)
        coeff=dict.fromkeys(product(I,repeat=3),1)
        for ij,v in zip(product((1,2),repeat=2),(a,b,c,d)):
            coeff[(*ij,1)]=1+v; coeff[(*ij,2)]=1+l+m*v
        assert all(coeff.values()) and D*l*m*(1+l-m)*(D+S)!=0
        for r in I:
            C=[Matrix(3,3,lambda i,j:coeff[(i,j,k)[-r:]+(i,j,k)[:-r]]) for k in I]
            L=Matrix(3,6,lambda k,j:prod(C[k][i,P[j][i]] for i in I))
            assert L.rank()==3
            assert all(L[:,[e,o]].rank()==2 for e in range(3) for o in range(3,6))
            assert any(L[:,[j for j in range(6) if j not in (e,o)]].rank()<3 for e in range(3) for o in range(3,6))
            for k in I:
                assert any(Matrix(3,3,lambda i,j:C[h][i,j]/C[k][i,j]).det()!=0 for h in I if h!=k)
                for h in I:
                    if h!=k:
                        assert Matrix(3,3,lambda i,j:C[h][i,j]/C[k][i,j]).rank()>1
        print(f'PASS: example {a,b,c,d,l,m}, S={S}, all three cuts (ii),(iii), not (i)',flush=True)

if __name__=='__main__':
    check_identities()
    check_examples()
