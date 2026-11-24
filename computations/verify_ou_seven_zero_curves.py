#!/usr/bin/env python3
"""Exact plane-cut identities and singular-point checks for two seven-zero families."""
from sympy import Matrix, Rational, symbols, expand, factor, groebner, diff, hessian, prod
from verify_ou_zero_supports import I, PAIRS, normalize


def determinant_matrix(c):
    rows=[]
    for k in I:
        C=Matrix(3,3,lambda i,j:c[i,j,k])
        rows.append([C[0,0]*C[1,1]*C[2,2],-C[0,0]*C[1,2]*C[2,1],
                     -C[0,2]*C[1,1]*C[2,0],C[0,1]*C[1,2]*C[2,0],
                     C[0,2]*C[1,0]*C[2,1],-C[0,1]*C[1,0]*C[2,2]])
    return Matrix(rows)


def main():
    u,x,y,z,v=symbols('u0 x y z v')
    coords=(x,y,z)
    for mask,free_entry,expected_L,expected_map,expected_F,exception in (
        (116309742,(2,1,1),[[0,0,-1,1,1,0],[1,0,-1,0,u,0],[0,-1,0,0,0,-1]],
         [x+(1-u)*y,-z,x+y,x,y,z],x*x*y+(1-u)*x*y*y+(x+y)*z*z,1),
        (116834030,(2,2,2),[[0,0,-1,1,1,0],[1,0,0,0,1,0],[0,-1,0,1,0,-u]],
         [-y,x-u*z,x+y,x,y,z],u*(x+y)*z*z-x*x*z-x*y*y-x*y*z,Rational(1,16)),
    ):
        c,K,free=normalize(mask)
        assert free==[free_entry]
        c={t:K.to_sympy(q) for t,q in c.items()}
        L=determinant_matrix(c)
        assert L==Matrix(expected_L)
        assert L[:,[0,1,2]].det() in (1,-1)
        p=Matrix(expected_map)
        assert all(expand(a)==0 for a in L*p)
        F=expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5])
        assert expand(F-expected_F)==0
        # Invert u to restrict to the stated coefficient torus. Every singular
        # point is in one of the three affine charts. A univariate polynomial
        # in the Groebner basis vanishes at every singular parameter.
        for t in coords:
            equations=[diff(F,s).subs(t,1) for s in coords]+[v*u-1]
            G=groebner(equations,*(s for s in coords if s!=t),v,u,order='lex')
            p0=[g.as_expr() for g in G.polys if g.as_expr().free_symbols<={u}]
            assert p0 and (p0[0]==1 or factor(p0[0]/(u-exception)).is_constant() or
                           factor(p0[0]/(u-exception)**2).is_constant())
        if mask==116834030:
            f=F.subs(u,exception)
            for t,expected in [(x,[y+2,z-8]),(y,[x+Rational(1,2),z+4]),
                               (z,[x-Rational(1,8),y+Rational(1,4)])]:
                G=groebner([diff(f,s).subs(t,1) for s in coords],*(s for s in coords if s!=t))
                assert list(G)==expected
            point={x:1,y:-2,z:8}
            assert hessian(f.subs(x,1),(y,z)).subs(point).det()==Rational(1,4)
            assert p.subs({**point,u:exception})==Matrix([2,Rational(1,2),-1,1,-2,8])
            for k in I:
                C=Matrix(3,3,lambda i,j:c[i,j,k].subs(u,exception))
                assert 0 in C and all(any(C.row(i)) for i in I) and all(any(C.col(j)) for j in I)
        print('PASS: mask',mask,'smooth away from',exception,flush=True)

if __name__=='__main__':
    main()
