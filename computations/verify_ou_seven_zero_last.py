#!/usr/bin/env python3
"""Exact rank and formal-arc identities for the remaining seven-zero supports.

The formal-arc argument, periodic simples, and bounded-head quotient
obstruction are proved in cubic_zero_supports.tex.  The finite checks alone
are not a proof of very strongness or its failure for the last special value.
"""
from itertools import product
from sympy import Matrix, QQ, symbols, cancel, Poly, eye, ones, zeros, BlockMatrix
from verify_ou_zero_supports import I, PAIRS, TRIPLES, normalize, path_dim
from verify_ou_seven_zero_curves import determinant_matrix


def order_at_zero(f, t):
    num, den = cancel(f).as_numer_denom()
    def order(p):
        return min(m[0] for m, c in Poly(p,t).terms() if c)
    return order(num)-order(den)


def main():
    u,t=symbols('u0 t')
    E,J,Z=eye(3),ones(3),zeros(3)
    D=BlockMatrix([[E,J,2*J],[Z,E,J],[Z,Z,E]]).as_explicit()
    S=D.inv()*D.T
    assert ((S**2).T*D)[:3,6:9]==15*J
    assert ((S**3).T*D)[:3,:3]==E+18*J
    c,K,free=normalize(117368558)
    assert free==[(2,2,2)]
    p,q,r=(0,1,2),(1,0,2),(2,0,1)
    c={(p[i],q[j],r[k]):K.to_sympy(v) for (i,j,k),v in c.items()}
    assert all(c[i,i,i]==0 for i in I)
    assert all(c[a]!=0 for a in TRIPLES if len(set(a))==2)
    assert determinant_matrix(c)==Matrix([[0,0,-1,0,1,-1],[0,-1,0,1,0,-u],[0,-1,-1,0,0,0]])
    # The arc works over Q(u), and separately at the exceptional value -1.
    # The only parameter-dependent leading factors are u+1 (and powers of it).
    R=1+t**3
    for value,e in [(u,1),(-1,3)]:
        v=cancel((R-1)*(R+value)/R**2)
        z=Matrix([1,R*v,-R*v,(R+value)*v,(1-R)*v,v])
        assert all(cancel(a)==0 for a in determinant_matrix(c).subs(u,value)*z)
        assert cancel(z[0]*z[3]*z[4]-z[1]*z[2]*z[5])==0
        X=Matrix.eye(3)
        X[0,1]=X[1,2]=t**e
        X[2,0]=cancel(z[3]/t**(2*e))
        X[1,0]=cancel(z[5]/X[0,1])
        X[2,1]=cancel(z[1]/X[1,2])
        X[0,2]=cancel(z[2]/X[2,0])
        assert all(order_at_zero(X[i,j],t)>0 for i,j in PAIRS if i!=j)
        assert all(order_at_zero(X[i,(i+1)%3],t)==e for i in I)
        for k in I:
            M=Matrix(3,3,lambda i,j:c[i,j,k].subs(u,value)*X[i,j])
            assert cancel(M.det())==0
            J=[j for j in I if j!=k]
            N=M.extract(J,J)
            assert order_at_zero(N.det(),t)==0
            B=-N.inv()*M.extract(J,[k])
            Y=Matrix(3,1,lambda i,j:1 if i==k else B[J.index(i)])
            assert all(cancel(a)==0 for a in M*Y)
            assert all(cancel(Y[j]).subs(t,0)==0 for j in J)
            assert order_at_zero(Y[(k-1)%3],t)==e
        print('PASS formal arc, parameter:',value,'head bound:',e,flush=True)
    c,K,free=normalize(125531767)
    assert free==[(2,2,2)]
    # Exact fraction-field computation: valid upper bound at every
    # specialization, not a rank estimate from sample parameters.
    assert path_dim(c,8,1,0,0,K)==17
    print('PASS 125531767: generic length-eight quotient 17 > 15',flush=True)
    c,K,free=normalize(116309742)
    assert free==[(2,1,1)]
    c={a:QQ.convert(K.to_sympy(v).subs(u,1)) for a,v in c.items()}
    # Further necessary-condition checks for the unresolved parameter only.
    for m,base in [(8,15),(9,18)]:
        for rot in I:
            for i,j in PAIRS:
                assert path_dim(c,m,rot,i,j,QQ)==base+int(m==9 and i==j)
    print('PASS 116309742 u=1: lengths eight and nine match Euler dimensions (not a sufficiency test)',flush=True)

if __name__=='__main__':
    main()
