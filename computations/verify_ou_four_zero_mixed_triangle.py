#!/usr/bin/env python3
"""Triangle presentations on equal and complementary loci of 134135535."""
from itertools import combinations, product
from sympy import Matrix, QQ, symbols, expand, cancel, groebner, Poly, ones, Integer
from verify_ou_zero_supports import normalize, TRIPLES, I, PAIRS, path_dim
from verify_ou_seven_zero_curves import determinant_matrix


def main():
    u,v,w,q=symbols('u0 u1 u2 u3')
    x,y,z,h=symbols('x y z h')
    c,K,free=normalize(134135535)
    assert free==[(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    c={a:K.to_sympy(b) for a,b in c.items()}
    ident=(0,1,2);cycle=(2,0,1);ks=(1,2,0)
    C0=Matrix([[1,1,1],[1,0,1],[1,1,0]])
    C1=Matrix([[1,1,1],[1,1,0],[1,0,1]])
    for kind in ('equal','diagonal','row','column'):
        if kind=='equal':
            r,P,Q,J=1,ident,ident,ident
            sub={v:u,w:u,q:u}
            C=[C0,C1,Matrix([[1,1,1],[1,u,u],[1,u,u]])]
        elif kind=='diagonal':
            r,P,Q,J=1,ident,ident,ident
            sub={v:1-u,w:1-u,q:u}
            C=[C0,C1,Matrix([[1,1,1],[1,u,1-u],[1,1-u,u]])]
        elif kind=='row':
            r,P,Q,J=2,ident,cycle,ks
            sub={v:u,w:1-u,q:1-u}
            C=[Matrix([[1,1,1],[u,0,1],[u,1,0]]),Matrix([[1,1,1],[1-u,1,0],[1-u,0,1]]),ones(3)]
        else:
            r,P,Q,J=0,cycle,ident,ks
            sub={v:1-u,w:u,q:1-u}
            C=[Matrix([[1,u,u],[1,0,1],[1,1,0]]),Matrix([[1,1-u,1-u],[1,1,0],[1,0,1]]),ones(3)]
        R={a:c[a[-r:]+a[:-r]].subs(sub,simultaneous=True) for a in TRIPLES}
        R={(i,j,k):R[P[i],Q[j],J[k]] for i,j,k in TRIPLES}
        assert C==[Matrix(3,3,lambda i,j:R[i,j,k]) for k in I]
        H=(1-u)/u if kind in ('row','column') else Integer(1)
        m=u if kind in ('row','column') else Integer(1)
        T=m*(x+y);a=H*(x+y)-z
        p=Matrix([H**2*T,T,a,x,y,z])
        assert all(cancel(f)==0 for f in determinant_matrix(R)*p)
        assert cancel(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]-T*(H*x-z)*(H*y-z))==0
        subs=[{y:-x},{z:H*y},{x:z/H}]
        coords=[(x,z),(x,y),(y,z)]
        X=[Matrix([[0,1,1],[1,-z,x],[1,-x,z]]),
           Matrix([[H*T,1,1],[H*y,H,1],[x,1/H,1]]),
           Matrix([[m*(z/H+y),z/H,y],[1,H,1],[1,1,H]])]
        Y=[[Matrix([x,1,-1]),Matrix([z,1,-1]),None],
           [Matrix([1,-x,-y]),Matrix([1,-y,-x]),None],
           [Matrix([-1,1,1]),Matrix([-1,1,1]),None]]
        if kind=='equal':
            Y[0][2]=Matrix([u*(x+z),1,-1]);Y[1][2]=Matrix([0,1,-1]);Y[2][2]=Matrix([u*(z-y),y-u*(y+z),u*(y+z)-z])
        elif kind=='diagonal':
            Y[0][2]=Matrix([u*z+(1-u)*x,1,-1]);Y[1][2]=Matrix([2*u-1,(1-u)*x-u*y,(1-u)*y-u*x]);Y[2][2]=Matrix([-1,1,1])
        else:
            Y[0][2]=Matrix([x+z,1,-1]);Y[1][2]=Matrix([H-1/H,x-H*y,y-H*x]);Y[2][2]=Matrix([-1/u,1,1])
            if kind=='row':
                Y[0][0]=Matrix([x/u,1,-1]);Y[0][1]=Matrix([z/(1-u),1,-1]);Y[1][0]=Matrix([-1/(1-u),x,y]);Y[1][1]=Matrix([-1/(1-u),y,x]);Y[2][0]=Y[2][1]=Matrix([-1/u,1,1])
            else:
                Y[1][0]=Matrix([-1/H,x,y]);Y[1][1]=Matrix([-1,y,x]);Y[2][1]=Matrix([-H,1,1])
        # In the first two cases z=y on L1; the kernels for k=0,1 have
        # negative sign on the degree-one entries.
        nonzero=u*(2*u-1)*(1 if kind=='equal' else 1-u)
        def unit(eq,vs=()):
            assert list(groebner([cancel(e).as_numer_denom()[0] for e in eq]+[h*nonzero-1],h,*vs,u))==[1],(kind,eq)
        def perms(A):
            return Matrix([A[0,0]*A[1,1]*A[2,2],A[0,0]*A[1,2]*A[2,1],A[0,2]*A[1,1]*A[2,0],A[0,1]*A[1,2]*A[2,0],A[0,2]*A[1,0]*A[2,1],A[0,1]*A[1,0]*A[2,2]])
        for l,(A,B,S,vs) in enumerate(zip(X,Y,subs,coords)):
            assert all(cancel(f)==0 for f in perms(A)-p.subs(S))
            for k in I:
                M=C[k].multiply_elementwise(A)
                assert cancel(M.det())==0 and all(cancel(f)==0 for f in M*B[k]),(kind,l,k,M*B[k])
                for j in (0,1):
                    ms=[M.extract(rr,cc).det().subs(vs[j],1) for rr,cc in product(combinations(I,2),repeat=2)]
                    unit(ms,(vs[1-j],));unit([f.subs(vs[j],1) for f in B[k]],(vs[1-j],))
        for P0 in [{x:1,y:-1,z:-H},{x:1,y:-1,z:H},{x:1,y:1,z:H}]:unit([(a*z).subs(P0,simultaneous=True)])
        Xg=Matrix([[H**2*T/z,1,1],[1,1,x/a],[1,y/a,z/a]])
        X0=X[0].copy()
        for i,j in product((1,2),(1,2)):X0[i,j]/=-z
        X1=X[1].copy()
        for i in I:X1[i,0]/=H*x
        X2=X[2].copy()
        for j in I:X2[0,j]/=H*y
        charts=[X0,Matrix.diag(H,1,H*y/x)*X1*Matrix.diag(x/y,1/H,1/H),
                Matrix.diag(H**3*y/z,1,1)*X2*Matrix.diag(1,1/H,z/(H**2*y))]
        for A,S in zip(charts,subs):assert all(cancel(f)==0 for f in A-Xg.subs(S)),(kind,A-Xg.subs(S))
        for k in I:
            for i in (1,2):
                P0=[expand(X[0][i,j]*Y[0][k][j]) for j in (1,2)]
                unit([Matrix([[f.coeff(x),f.coeff(z)] for f in P0]).det()])
                unit(Poly(X[0][i,0]*Y[0][k][0],x,z).coeffs())
            for j in (1,2):assert X[0][0,j]*Y[0][k][j]!=0
            l=2 if kind=='equal' and k==2 else 1
            unit(Poly(X[l][0,0]*Y[l][k][0],*coords[l]).coeffs())
        U=Matrix([[1,0,-1],[0,0,0],[0,0,0]])
        V=Matrix([[0,1,0],[1,0,0],[1,0,0]])
        W=Matrix([[1,1,0],[1,1,0],[1,0,1] if kind=='equal' else [1,1,0]])
        A,B,Csum=[Matrix([[1,1,1]])*M for M in (U,V,W)]
        assert A+Csum-2*B==Matrix([[0,0,0] if kind=='equal' else [0,1,-1]])
        for M,N in ((U,V),(V,W),(U,W),(U,U),(V,V),(W,W)):
            for i,j in PAIRS:
                d=N.row(j)-M.row(i)
                assert all(sum(d[k] for k in I if mask>>k&1)>=-1 for mask in range(1,7))
        print('PASS',kind,'triangle and multidegree',A+Csum-2*B,flush=True)
    c={a:QQ.convert(b.subs({u:QQ(1,2),v:QQ(1,2),w:QQ(1,2),q:QQ(1,2)},simultaneous=True)) for a,b in c.items()}
    assert path_dim(c,5,0,2,0,QQ)==8
    print('PASS all-half point: length-five quotient 8 > 7')

if __name__=='__main__':
    main()
