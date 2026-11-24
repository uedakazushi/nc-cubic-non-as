#!/usr/bin/env python3
"""Exact checks for the three seven-zero line-conic representations.

The general gluing, constant-rank, cohomology and multidegree obstruction
are proved in cubic_zero_supports.tex.  No finite-field inference is used.
"""
from itertools import product
from sympy import Matrix, QQ, symbols, expand, cancel, hessian
from verify_ou_zero_supports import I, TRIPLES, normalize, path_dim
from verify_ou_seven_zero_curves import determinant_matrix


def main():
    u, x, y, z, s, t = symbols('u0 x y z s t')
    U = Matrix([[1,-1],[0,0],[0,0]])
    V = Matrix([[0,1],[1,0],[1,0]])
    W = Matrix([[1,1],[1,1],[0,2]])
    A,B,C = [Matrix([[1,1,1]])*Z for Z in (U,V,W)]
    assert A+C-2*B == Matrix([[-1,1]])
    for X,Y in [(U,V),(V,W),(U,W)]:
        for i,j in product(I,repeat=2):
            assert min(Y.row(j)-X.row(i)) >= -1
    for X in (U,V,W):
        for i,j in product(I,repeat=2):
            assert min(X.row(j)-X.row(i)) >= -1
    for mask in (128957047,129972855,133118583):
        c,K,free=normalize(mask)
        assert free==[(2,2,2)]
        c={a:K.to_sympy(q) for a,q in c.items()}
        r={a:c[a[-2:]+a[:-2]] for a in TRIPLES}
        Cs=[Matrix(3,3,lambda i,j:r[i,j,k]) for k in I]
        eps=int(mask==128957047)
        assert Cs[2]==Matrix([[1,eps,1],[0,1,1],[0,1,u]])
        assert Cs[0]==Matrix([[1,1,1],[1,0,1],[1,1,int(mask==129972855)]])
        assert Cs[1]==Matrix([[0,1,1],[1,1,int(mask==133118583)],[1,0,1]])
        xi=(x+y-(z if mask==129972855 else 0))/u
        a=x-z if mask==133118583 else -z
        b,c0,d=x,y,z
        p=Matrix([xi,u*xi,a,b,c0,d])
        L=determinant_matrix(r)
        assert all(expand(v)==0 for v in L*p)
        assert L[:,[0,1,2]].det() in (u,-u)
        Q=expand(b*c0-u*a*d)
        assert hessian(Q,(x,y,z)).det()!=0
        assert expand(p[0]*p[3]*p[4]-p[1]*p[2]*p[5]-xi*Q)==0
        # On the line all pairs (a,b),(c,d) are independent linear forms.
        line_y=z-x if mask==129972855 else -x
        for f,g in [(a,b),(c0,d)]:
            v,w=f.subs(y,line_y),g.subs(y,line_y)
            assert Matrix([[v.coeff(x),v.coeff(z)],[w.coeff(x),w.coeff(z)]]).det()!=0
        # Conic parametrization a=l1*m1, b=l1*m2, c=l2*m1, d=l2*m2/u.
        l1,l2=s,t
        if mask==133118583:
            m1,m2=u*s-t,u*s
            xyz={x:u*s*s,y:t*(u*s-t),z:s*t}
        else:
            m1,m2=t,-u*s
            xyz={x:-u*s*s,y:t*t,z:-s*t}
        assert all(expand(f.subs(xyz)-g)==0 for f,g in
                   zip((a,b,c0,d),(l1*m1,l1*m2,l2*m1,l2*m2/u)))
        assert expand(Q.subs(xyz))==0
        assert Matrix([[m1.coeff(s),m1.coeff(t)],[m2.coeff(s),m2.coeff(t)]]).det()!=0
        h=expand(xi.subs(xyz))
        # The two gluing charts meet every intersection: h has no common
        # zero with any of l1,l2,m1,m2, for nonzero u.
        from sympy import resultant
        for f in (l1,l2,m1,m2):
            # Binary resultant computed by the homogeneous evaluation.
            res=expand(h.subs({s:-f.coeff(t),t:f.coeff(s)},simultaneous=True))
            assert res.as_coeff_exponent(u)[0] != 0 and res.as_coeff_exponent(u)[0].is_Rational
        Xl=Matrix([[0,1,1],[1,a,b],[1,c0,d]])
        Xc=Matrix([[u*h,m2,m1],[l2,1,1],[l1,1,1/u]])
        for k in I:
            assert expand(Cs[k].multiply_elementwise(Xl).det().subs(y,line_y))==0
            assert cancel(Cs[k].multiply_elementwise(Xc).det())==0
        for k in (0,1):
            assert cancel(Cs[k].multiply_elementwise(Xc)[1:,1:].det())!=0
            assert Cs[k].multiply_elementwise(Xl)[:2,:2].det()==-1
        h0=m1-eps*m2
        vec=Matrix([h0,u*h,-u*h])
        assert all(expand(v)==0 for v in Cs[2].multiply_elementwise(Xc)*vec)
        res=expand(h.subs({s:-h0.coeff(t),t:h0.coeff(s)},simultaneous=True))
        assert cancel(res/(u-1 if eps else 1)).as_coeff_exponent(u)[0].is_Rational
        # Equality on the conic chart. All linear forms here are in one
        # common local trivialization of O_{P1}(1).
        Xg=Matrix([[h/(l2*m2/u),1,1],[1,1,m2/m1],
                   [1,l2/l1,l2*m2/(u*l1*m1)]])
        R=Matrix.diag(1/m2,1,l2/l1)
        S=Matrix.diag(1/l2,1,m2/m1)
        assert all(cancel(v)==0 for v in R*Xc*S-Xg)
        print('PASS line-conic identities and multidegrees:',mask,flush=True)
    # In the one exceptional specialization the sixth-path dimension jumps.
    c,K,_=normalize(128957047)
    c={a:QQ.convert(K.to_sympy(q).subs(u,1)) for a,q in c.items()}
    assert path_dim(c,6,0,0,0,QQ)==11
    print('PASS 128957047 u=1: length 6, quotient 11 > 10',flush=True)

if __name__=='__main__':
    main()
