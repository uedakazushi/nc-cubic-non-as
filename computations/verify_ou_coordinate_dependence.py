#!/usr/bin/env python3
"""Exact polynomial identities for the residual dense coordinate-dependence locus.

The length-seven fraction-field rank proves an upper bound on rank at every
specialization of the polynomial matrix. The selected minors and all-minor
identities are checked over Q[a,b,c,d,lambda,mu].
"""
from itertools import combinations, permutations, product
from sympy import Matrix, QQ, prod, symbols, expand, cancel
from ou_path_rank import path_matrix, rotate, exact_rank, TRIPLES


def check_coordinate_identities():
    a,b,c,d,l,m=symbols('a b c d l m')
    t=1+l-m
    D=a*d-b*c
    S=a+d-b-c
    T=a*b*c-a*b*d-a*c*d+b*c*d
    F1=(m-1)*D+l*S
    F3=m*T-(l+1)*D
    F2=F3-F1
    assert expand(F2-(m*T-(l+m)*D-l*S)) == 0
    I=range(3)
    P=sorted(permutations(I),key=lambda p:sum(p[i]>p[j] for i in I for j in range(i+1,3))%2)
    coeff=dict.fromkeys(TRIPLES,1)
    for ij,x in zip(product((1,2),repeat=2),(a,b,c,d)):
        coeff[(*ij,1)]=1+x
        coeff[(*ij,2)]=1+l+m*x
    r=rotate(coeff,1)
    L=Matrix(3,6,lambda k,j:prod(r[i,P[j][i],k] for i in I))
    A=lambda x:1+x
    B=lambda x:1+l+m*x
    assert L==Matrix([[1]*6,[A(a)*B(b),B(a),A(b),B(a)*A(b),B(b),A(a)],
                            [A(c)*B(d),B(c),A(d),B(c)*A(d),B(d),A(c)]])
    # Delete columns e,o and take the three-column minor specified by the index.
    # The remaining columns are kept in increasing order; the four minors
    # correspond to triples 012,013,023,123.
    selected={(0,3):(1,F1),(0,4):(1,F1),(0,5):(1,-m*F1),
              (1,3):(3,-F1),(1,4):(1,-F2),(1,5):(0,t*F2),
              (2,3):(3,m*F1),(2,4):(0,t*F3),(2,5):(0,t*F3)}
    for (e,o),(idx,v) in selected.items():
        C=L[:,[j for j in range(6) if j not in (e,o)]]
        minors=[expand(C[:,list(js)].det()) for js in combinations(range(4),3)]
        assert expand(minors[idx]-v)==0, (e,o,idx,expand(minors[idx]-v))
        if (e,o) in [(0,3),(1,4),(2,5)]:
            f={(0,3):F1,(1,4):F2,(2,5):F3}[e,o]
            assert all(cancel(x/f).is_polynomial(a,b,c,d,l,m) for x in minors)
    # Three determinants sufficient for (iii) in the (j,k,i) rotation.
    ratio=lambda k,h:Matrix(3,3,lambda i,j:r[i,j,h]/r[i,j,k])
    assert cancel(ratio(0,1).det()-l*(a-b))==0
    assert cancel(ratio(1,2).det()+m*(a-b)*(a-c)*(b-d)*t/(A(a)*A(b)*B(a)*B(b)))==0
    assert cancel(ratio(2,1).det()+m*(a-c)*(b-d)*(c-d)*t/(A(c)*A(d)*B(c)*B(d)))==0
    # Transposition of M exchanges the other two rotations, after reordering columns.
    # All F_i are invariant under b<->c.
    for f in (F1,F2,F3):assert expand(f.xreplace({b:c,c:b})-f)==0
    print('PASS: selected minors, three all-minor divisibilities, determinants and transpose identities',flush=True)


def check_zero_slope_rank():
    vs=symbols('a b c d q')
    K=QQ.frac_field(*vs)
    a,b,c,d,q=map(K.from_sympy,vs)
    C=dict.fromkeys(TRIPLES,K.one)
    for ij,x in zip(product((1,2),repeat=2),(a,b,c,d)):
        C[(*ij,1)]=1+x
        C[(*ij,2)]=1+q
    n,rows=path_matrix(C,7,1,1,K)
    nr=len(rows)
    rank=exact_rank(rows,K)
    assert (nr,n,rank)==(1458,729,716),(nr,n,rank)
    print('PASS: constant second face, length 7, exact rank 716, quotient dimension 13',flush=True)

if __name__=='__main__':
    check_coordinate_identities()
    check_zero_slope_rank()
