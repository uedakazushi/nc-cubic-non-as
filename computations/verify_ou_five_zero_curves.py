#!/usr/bin/env python3
"""Exact certificates for the five-zero smooth-grid and periodic-arc loci.
The paper proves the gluing and formal lifting arguments over characteristic zero.
"""
from itertools import combinations, product
from sympy import Matrix, symbols, Poly, expand, factor, cancel, groebner
from verify_ou_zero_supports import I, PAIRS, TRIPLES, normalize
from verify_ou_seven_zero_curves import determinant_matrix

u,v,w,x,y,z,h,t = symbols('u0 u1 u2 x y z h t')
XYZ=Matrix([x,y,z])


def tensor(mask):
    c,K,free=normalize(mask)
    assert len(free)==3
    return {a:K.to_sympy(b) for a,b in c.items()}


def rotate(c,r):
    return {a:c[a[-r:]+a[:-r]] for a in TRIPLES}


def relabel(c,P):
    return {(i,j,k):c[P[i][0],P[j][1],P[k][2]] for i,j,k in TRIPLES}


def unit(q):
    return q!=0 and len(Poly(expand(q),u,v,w).terms())==1


def check_grid():
    cases=[(134135534,0,(u-1)*(w-1)),(134147822,0,w-1),
           (117374702,0,u-1),(132054766,0,w-1),
           (133873391,1,w-1),(133873391,1,u-v),
           (133873391,0,v-1),(133873391,0,u-w),
           (133099247,1,1),(132054255,1,1)]
    for mask,r,extra in cases:
        c=rotate(tensor(mask),r);L=determinant_matrix(c)
        piv=next(js for js in combinations(range(6),3) if unit(L[:,list(js)].det()))
        oth=[j for j in range(6) if j not in piv];p=Matrix.zeros(6,1)
        for j,a in zip(oth,XYZ):p[j]=a
        for j,a in zip(piv,-L[:,list(piv)].inv()*L[:,oth]*XYZ):p[j]=factor(a)
        assert all(cancel(a)==0 for a in L*p)
        f=cancel(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]);num,den=f.as_numer_denom()
        assert unit(den)
        f=expand(num)
        for k in I:
            Z=[(i,j) for i,j in PAIRS if c[i,j,k]==0]
            assert 1<=len(Z)<=3 and len(set(a[0] for a in Z))==len(set(a[1] for a in Z))==len(Z)
        row=lambda l:Matrix([expand(l).coeff(a) for a in XYZ])
        for i,j in product((0,3,4),(1,2,5)):
            P=row(p[i]).cross(row(p[j]))
            G=[cancel(f.diff(a).subs(dict(zip(XYZ,P)),simultaneous=True)).as_numer_denom()[0] for a in XYZ]
            assert list(groebner(G+[h*u*v*w*extra-1],h,u,v,w))==[1]
        print('PASS grid',mask,r,extra,'pivot',piv,flush=True)


def check_arcs():
    p,a,b=symbols('p a b')
    for mask,P,LL,delta in [
      (125763310,[(0,0,0),(1,2,1),(2,1,2)],
       [[0,0,-1,1,1,-1],[0,-v,0,1,0,-u],[0,-w,-1,0,w,0]],-u*w-v),
      (132086511,[(0,2,2),(1,0,1),(2,1,0)],
       [[0,0,-w,u*w,1,-1],[0,-1,0,v,0,-1],[0,-1,-1,1,0,0]],-w-1)]:
        c=relabel(tensor(mask),P);L=determinant_matrix(c)
        assert L==Matrix(LL)
        assert all(c[i,i,i]==0 for i in I)
        assert all(c[i,j,k]!=0 for i,j,k in TRIPLES if len({i,j,k})==2)
        assert expand(L[:,[1,2,5]].det()-delta)==0
        # Generic arc: set the forward arrows to t and reverse arrows to t^2 b_i.
        r0,r1,r2=symbols('r0 r1 r2');R=[r0,r1,r2]
        A=Matrix.eye(3)
        for i in I:A[i,(i+1)%3]=t;A[(i+1)%3,i]=t*t*R[i]
        zz=Matrix([1,A[1,2]*A[2,1],A[0,2]*A[2,0],A[0,1]*A[1,2]*A[2,0],A[0,2]*A[1,0]*A[2,1],A[0,1]*A[1,0]])
        E=(L*zz/t**3).applyfunc(cancel)
        J=E.jacobian(R).subs(t,0)
        assert expand(J.det()**2-delta**2)==0
        assert all(cancel(Matrix(3,3,lambda i,j:c[i,j,k]*A[i,j]).det()-t**3*E[k])==0 for k in I)
        if mask==125763310:
            # v=-u*w. The linear equations imply (u-1)*a/u+(1-w)*b=0.
            zz=Matrix([1,p,w*(b-p),a,b,w*p+a/u])
            assert all(cancel(e)==0 for e in L.subs(v,-u*w)*zz-Matrix([(u-1)*a/u+(1-w)*b,0,0]))
            for sub,vec in [({u:1,v:-w},Matrix([1,0,0,t**3,0,t**3])),({w:1,v:-u},Matrix([1,0,t**3,0,t**3,0]))]:
                assert all(cancel(e)==0 for e in L.subs(sub,simultaneous=True)*vec)
                assert cancel(vec[3]*vec[4]-vec[1]*vec[2]*vec[5])==0
            lam=u*(w-1)/(u-1)
            zz=zz.subs({p:t**2,b:t**3*h,a:lam*t**3*h},simultaneous=True).applyfunc(cancel)
            assert all(cancel(e)==0 for e in L.subs(v,-u*w)*zz)
            H=lam*h**2-w*(t*h-1)*(w+lam*t*h/u)
        else:
            # w=-1. Set b=(u+v-1)*a.
            lam=u+v-1
            zz=Matrix([1,p,a-p,a,lam*a,v*a-p])
            assert all(cancel(e)==0 for e in L.subs(w,-1)*zz)
            vec=zz.subs({p:0,a:t**3,u:1-v},simultaneous=True)
            assert cancel(vec[3]*vec[4]-vec[1]*vec[2]*vec[5])==0
            zz=zz.subs({p:t**2,a:t**3*h},simultaneous=True).applyfunc(cancel)
            H=lam*h**2-(t*h-1)*(v*t*h-1)
        assert cancel((zz[3]*zz[4]-zz[1]*zz[2]*zz[5])/t**6-H)==0
        # Both branches lift by the scalar implicit-function theorem at a
        # nonzero root of H(h,0). The nonzero factors are stated in the paper.
        assert cancel(H.diff(h).subs(t,0)-2*lam*h)==0
        print('PASS periodic arcs',mask,flush=True)

if __name__=='__main__':
    check_grid()
    check_arcs()
