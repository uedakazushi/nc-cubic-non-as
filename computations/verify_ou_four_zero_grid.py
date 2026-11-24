#!/usr/bin/env python3
"""Exact localized gradient ideals and rank-two plane cuts for four-zero tensors."""
from itertools import product
from sympy import Matrix, Poly, symbols, expand, cancel, groebner
from verify_ou_zero_supports import normalize, TRIPLES, PAIRS, I
from verify_ou_seven_zero_curves import determinant_matrix

u,v,w,q,h=symbols('u0 u1 u2 u3 h')
x,y,z=symbols('x y z')
XYZ=(x,y,z)
PARAMS=(u,v,w,q)


def tensor(mask,r):
    c,K,free=normalize(mask)
    assert len(free)==4
    c={a:K.to_sympy(b) for a,b in c.items()}
    return {a:c[a[-r:]+a[:-r]] for a in TRIPLES}


CASES=[
    (133660399,0,(0,1,3),u+w-1),
    (133660399,2,(0,1,3),u+w),
    (133889775,0,(0,1,3),u+w-1),
    (133889775,1,(0,1,2),q-1),
    (133889775,2,(0,1,2),v-1),
    (134147823,0,(0,1,2),(v-1)*(q-1)*(u+w)),
    (134147823,0,(0,1,2),(v-1)*(q-1)*(v-q)),
    (134151407,0,(0,1,5),v-1),
    (134151918,0,(0,1,2),(u-1)*(q-1)),
    (134183663,0,(0,1,3),1),
]


def check_grid():
    for mask,r,piv,e in CASES:
        c=tensor(mask,r); L=determinant_matrix(c)
        for k in I:
            zz=[(i,j) for i,j in PAIRS if c[i,j,k]==0]
            assert zz and len({i for i,j in zz})==len({j for i,j in zz})==len(zz)
        det=L[:,list(piv)].det()
        assert det!=0 and len(Poly(det,*PARAMS).terms())==1
        other=[i for i in range(6) if i not in piv]
        p=Matrix.zeros(6,1)
        for i,a in zip(other,XYZ):p[i]=a
        for i,a in zip(piv,-L[:,list(piv)].inv()*L[:,other]*Matrix(XYZ)):p[i]=cancel(a)
        assert all(cancel(a)==0 for a in L*p)
        num,den=cancel(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]).as_numer_denom()
        assert len(Poly(den,*PARAMS).terms())==1
        f=expand(num)
        row=lambda a:Matrix([expand(a).coeff(b) for b in XYZ])
        for i,j in product((0,3,4),(1,2,5)):
            P=row(p[i]).cross(row(p[j]))
            eq=[cancel(f.diff(a).subs(dict(zip(XYZ,P)),simultaneous=True)).as_numer_denom()[0] for a in XYZ]
            assert list(groebner(eq+[h*u*v*w*q*e-1],h,*PARAMS))==[1],(mask,r,e,i,j)
        print('PASS grid',mask,r,e,flush=True)


def check_rank_two():
    # 134135535: rank two exactly at the all-ones point in rotation 1.
    L=determinant_matrix(tensor(134135535,1))
    from itertools import combinations
    minors=[L[:,list(j)].det() for j in combinations(range(6),3)]
    assert list(groebner(minors+[h*u*v*w*q-1],h,*PARAMS))==[h-1,u-1,v-1,w-1,q-1]
    L=L.subs({u:1,v:1,w:1,q:1})
    p=Matrix([2*x+y+z,y+z,x,y,z,x+y+z])
    assert L.rank()==2 and all(cancel(a)==0 for a in L*p)
    assert all(Matrix.vstack(L[i,:],L[j,:]).rank()==2 for i,j in combinations(I,2))
    f=expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5])
    for a in XYZ:
        assert list(groebner([f.diff(b).subs(a,1) for b in XYZ],*[b for b in XYZ if b!=a]))==[1]
    print('PASS rank-two smooth plane 134135535',flush=True)
    # 133136111: rank two iff v=1,q=-1,w=1-u (u,w nonzero).
    for r,col in enumerate((3,2,5)):
        L=determinant_matrix(tensor(133136111,r))
        assert L[:,col]==Matrix.zeros(3,1)
        minors=[L[:,list(j)].det() for j in combinations(range(6),3)]
        G=groebner(minors+[h*u*v*w*q-1],h,*PARAMS)
        assert all(G.reduce(a)[1]==0 for a in (v-1,q+1,u+w-1))
        assert all(expand(a.subs({v:1,q:-1,w:1-u},simultaneous=True))==0 for a in minors)
    L=determinant_matrix(tensor(133136111,0))
    LL=L.subs({v:1,q:-1,w:1-u},simultaneous=True)
    A,B,C,D=symbols('A B C D')
    p=Matrix([A,C-u*A,A+B-C,D,B,C])
    assert all(expand(a)==0 for a in LL*p)
    assert all(Matrix.vstack(LL[i,:],LL[j,:]).rank()==2 for i,j in combinations(I,2))
    f=expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5])
    assert expand(f-(A*B*D-(C-u*A)*(A+B-C)*C))==0
    # On A=0, B!=0 the only possible singular point is C=D=0.
    G=groebner([f.diff(a).subs({A:0,B:1}) for a in (A,B,C,D)]+[h*u*(1-u)-1],C,D,h,u)
    assert all(G.reduce(a)[1]==0 for a in (C,D))
    # On B=0, A!=0 the cubic C(C-u)(1-C) has distinct roots.
    assert list(groebner([f.diff(a).subs({B:0,A:1}) for a in (A,B,C,D)]+[h*u*(1-u)-1],C,D,h,u))==[1]
    # On A=B=0, C!=0, f is nonzero; only [0:0:0:1] remains.
    assert f.subs({A:0,B:0,C:1})==1
    print('PASS rank-two finite singular locus 133136111',flush=True)


def check_full_slice():
    # 134135535 away from any parameter=1. Use signed maximal minors
    # of [L; e_i; e_j] for each grid point, so no chart denominator is needed.
    from itertools import combinations
    e=u*v*w*q*(u-1)*(v-1)*(w-1)*(q-1)
    pairings=[((u,w),(v,q)),((u,q),(v,w)),((u,v),(w,q))]
    full_indices=(0,2,0)
    for r,(pairs,k) in enumerate(zip(pairings,full_indices)):
        c=tensor(134135535,r); L=determinant_matrix(c)
        for l in I:
            zz=[(i,j) for i,j in PAIRS if c[i,j,l]==0]
            assert (not zz)==(l==k)
            assert len({i for i,j in zz})==len({j for i,j in zz})==len(zz)
        minors=[L[:,list(j)].det() for j in combinations(range(6),3)]
        assert list(groebner(minors+[h*e-1],h,*PARAMS))==[1]
        diffs=[a-b for a,b in pairs]
        for i,j in product((0,3,4),(1,2,5)):
            A=Matrix.vstack(L,Matrix.eye(6)[[i,j],:])
            p=Matrix([(-1)**a*A[:,[l for l in range(6) if l!=a]].det() for a in range(6)])
            grad=Matrix([[p[3]*p[4],-p[2]*p[5],-p[1]*p[5],p[0]*p[4],p[0]*p[3],-p[1]*p[2]]])
            B=Matrix.vstack(L,grad)
            eq=[expand(B[:,list(J)].det()) for J in combinations(range(6),4)]
            G=groebner(eq+[h*e-1],h,*PARAMS)
            assert all(G.reduce(expand(a**3))[1]==0 for a in diffs),(r,i,j)
        # A rank-one torus matrix forces permutation products proportional
        # to reciprocals of the six unsigned slice products.
        signs=(1,-1,-1,1,1,-1)
        p=Matrix([signs[j]/L[k,j] for j in range(6)])
        eq=[cancel(a).as_numer_denom()[0] for a in L*p]
        G=groebner(eq+[h*u*v*w*q-1],h,*PARAMS)
        assert all(G.reduce(a+b-1)[1]==0 for a,b in pairs)
        print('PASS mixed full-slice grid and torus rank',r,flush=True)


if __name__=='__main__':
    check_grid()
    check_rank_two()
    check_full_slice()
