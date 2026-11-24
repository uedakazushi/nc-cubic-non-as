#!/usr/bin/env python3
"""Exact checks for the seven residual six-zero support orbits.

For the five geometric families, verify plane, irreducibility, smooth grid,
and discriminant identities; the paper explains the geometric consequences.
For the other two families, verify support conditions and formal arcs.
"""
from itertools import product
from sympy import (Matrix, QQ, symbols, expand, cancel, factor, Poly, groebner,
                   diff, discriminant)
from verify_ou_zero_supports import I, PAIRS, TRIPLES, normalize
from verify_ou_seven_zero_curves import determinant_matrix

u,v,w,x,y,z,t,h = symbols('u0 u1 u2 x y z t h')
XYZ=(x,y,z)

def tensor(mask):
    c,K,free=normalize(mask)
    return {a:K.to_sympy(b) for a,b in c.items()},free

def rotate(c,r):
    return {a:c[a[-r:]+a[:-r]] for a in TRIPLES}

def relabel(c,P):
    return {(i,j,k):c[P[i][0],P[j][1],P[k][2]] for i,j,k in TRIPLES}

def check_curves():
    # All maps are valid at every nonzero u,v; rotate 131005679 once.
    # Polynomial f is F times the stated nonzero factor.
    cases=[
        (117358318,0,[(2,1,1),(2,2,2)],
         [x+(1-u)*y,x-v*z,x+y,x,y,z],1,u-1,
         -256*u**2*v**4*(u-1)**2*(16*u*v-16*v**2-8*v-1)),
        (134145774,0,[(2,2,1),(2,2,2)],
         [(x-y)/u,(-u*y+u*z+v*x-v*y)/u,x,z-y,y,z],u,1,
         -256*u**3*v**3*((u+v-1)**3+27*u*v)),
        (117370606,0,[(2,1,2),(2,2,2)],
         [x-y,(v*x-v*y-v*z-y)/u,x,-y,y,z],u,1,
         -256*u**3*v**5*(16*u**2*v**2-u*v**3-20*u*v**2+8*u*v+(v+1)**3)),
        (133099246,0,[(2,2,1),(2,2,2)],
         [(x-y)/u,v*(-u*z+x-y)/u,x,z-y,y,z],u,1,
         -256*u**3*v**5*(16*u**2*v**2-u*v**3+20*u*v**2+8*u*v+(1-v)**3)),
        (131005679,1,[(2,2,1),(2,2,2)],
         [(u*x+u*y-x)/v,x+y,x-z,x,y,z],v,1,
         256*u**2*v**6*(16*u**2+8*u*v-32*u+v**2+8*v+16)),
    ]
    for mask,r,free,p,fac,extra,disc0 in cases:
        c,fs=tensor(mask);assert fs==free
        c=rotate(c,r);L=determinant_matrix(c)
        assert all(expand(e)==0 for e in L*Matrix(p))
        assert any(factor(L[:,list(js)].det()) in (1,-1,u,-u,v,-v) for js in __import__('itertools').combinations(range(6),3))
        f=expand(fac*(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]))
        A,B,C=Poly(f,x).all_coeffs()
        ay,az=A.coeff(y),A.coeff(z)
        assert groebner([ay,az,h*u*v*extra-1],h,u,v).polys[0].as_expr()==1
        B0=expand(B.subs({y:-az,z:ay},simultaneous=True));C0=expand(C.subs({y:-az,z:ay},simultaneous=True))
        assert list(groebner([B0,C0,h*u*v*extra-1],h,u,v))==[1]
        delta=expand(B**2-4*A*C)
        assert expand(discriminant(delta.subs(z,1),y)-disc0)==0
        assert disc0.subs({u:2,v:3})!=0
        co=[Poly(delta,y,z).coeff_monomial(y**(4-i)*z**i) for i in range(5)]
        assert list(groebner([co[0],h*u*v*extra-1],h,u,v))==[1]
        b=co[1]/(2*co[0]);d=co[2]/(2*co[0])-b*b/2
        eqs=[factor(cancel(co[3]-2*b*d*co[0]).as_numer_denom()[0]),factor(cancel(co[4]-d*d*co[0]).as_numer_denom()[0])]
        G=groebner(eqs+[h*u*v*extra*co[0]-1],h,u,v)
        if mask==134145774:
            assert list(G)==[h-1,u+1,v+1]
        else:
            assert list(G)==[1]
        # Every opposite-parity pair spans two independent lines with a
        # smooth grid point, throughout the indicated parameter region.
        row=lambda l:Matrix([expand(l).coeff(a) for a in XYZ])
        for i,j in product((0,3,4),(1,2,5)):
            P=row(p[i]).cross(row(p[j]));sub=dict(zip(XYZ,P))
            eqs=[cancel(diff(f,a).subs(sub,simultaneous=True)).as_numer_denom()[0] for a in XYZ]
            assert list(groebner(eqs+[h*u*v*extra-1],h,u,v))==[1]
        for k in I:
            M=Matrix(3,3,lambda i,j:c[i,j,k]);assert 0 in M
            assert all(any(M.row(i)) and any(M.col(i)) for i in I)
        if mask==134145774:
            p0=[cancel(a.subs({u:-1,v:-1})) for a in p]
            assert expand(p0[0]*p0[3]*p0[4]-p0[1]*p0[2]*p0[5]+(x*z-y*y)*(x-y+z))==0
            for k in I:
                M=Matrix(3,3,lambda i,j:c[i,j,k]); Z=[(i,j) for i,j in PAIRS if M[i,j]==0]
                assert len(Z)==2 and Z[0][0]!=Z[1][0] and Z[0][1]!=Z[1][1]
        print('PASS geometric family',mask,'rotation',r,'disc',factor(disc0),flush=True)

def order(f):
    a,b=cancel(f).as_numer_denom()
    return min(i[0] for i,c in Poly(a,t).terms() if c)-min(i[0] for i,c in Poly(b,t).terms() if c)

def check_periodic():
    c,free=tensor(131037935); assert free==[(1,1,2),(2,2,1),(2,2,2)]
    c=relabel(c,[(0,1,1),(1,2,0),(2,0,2)])
    assert all(c[i,i,i]==0 and c[i,(i+1)%3,(i+2)%3]==0 for i in I)
    assert all(c[a]!=0 for a in TRIPLES if len(set(a))==2)
    print('PASS periodic support 131037935',flush=True)
    c,free=tensor(125759214);assert free==[(2,2,1),(2,2,2)]
    c=relabel(c,[(0,0,0),(1,2,1),(2,1,2)])
    assert all(c[i,i,i]==0 for i in I)
    assert all(c[a]!=0 for a in TRIPLES if len(set(a))==2)
    L=determinant_matrix(c)
    assert L==Matrix([[0,0,-1,0,1,-1],[0,-u,0,1,0,-1],[0,-v,-1,0,v,0]])
    # Cases u+v!=0, u=-v with v!=1, and the point (-1,1).
    for sub,e,step in [({},1,1),({u:-v},3,1),({u:-1,v:1},1,-1)]:
        X=Matrix.eye(3)
        if step==1:
            R=-v+t**3
            D=1-(v-1)*R/v
            s=cancel(((u+D)*(R+D)/(R*D)).subs(sub))
            zz=Matrix([1,s,R*s,(u+D).subs(sub)*s,(R+D)*s,D*s]).applyfunc(cancel)
            assert all(cancel(a)==0 for a in L.subs(sub)*zz)
            assert cancel(zz[0]*zz[3]*zz[4]-zz[1]*zz[2]*zz[5])==0
            X[0,1]=X[1,2]=t**e
            X[2,0]=cancel(zz[3]/t**(2*e))
            X[1,0]=cancel(zz[5]/t**e)
            X[2,1]=cancel(zz[1]/t**e)
            X[0,2]=cancel(zz[2]/X[2,0])
        else:
            X[0,2]=X[1,0]=X[2,1]=t;X[2,0]=t**2
        assert all(order(X[i,j])>=e for i,j in PAIRS if i!=j and X[i,j]!=0)
        assert all(order(X[i,(i+step)%3])==e for i in I)
        for k in I:
            M=Matrix(3,3,lambda i,j:c[i,j,k].subs(sub)*X[i,j]);assert cancel(M.det())==0
            J=[j for j in I if j!=k];N=M.extract(J,J)
            assert order(N.det())==0
            b=-N.inv()*M.extract(J,[k]);Y=Matrix([1 if i==k else b[J.index(i)] for i in I])
            assert all(cancel(a)==0 for a in M*Y)
            assert all(cancel(Y[j]).subs(t,0)==0 for j in J)
            assert order(Y[(k-step)%3])==e
        print('PASS periodic arc 125759214',sub,'order',e,flush=True)

if __name__=='__main__':
    check_curves()
    check_periodic()
