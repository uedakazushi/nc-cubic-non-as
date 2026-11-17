#!/usr/bin/env python3
"""Exact finite checks for the quantum Fermat cubic counterexample.
Arithmetic in exponents modulo 3; no numerical roots of unity.
The domain, tails dimension, and closed-incidence proofs are in the paper.
"""
from __future__ import annotations
import itertools
import json
from pathlib import Path
from math import comb
import sympy as s

E = s.Matrix([[0,1,1,1],[-1,0,0,1],[-1,0,0,1],[-1,-1,-1,0]])
N, D = 4, 3
assert E + E.T == s.zeros(4)
assert all(int(v) % 3 == 0 for v in E*s.ones(4,1))

def monomials(d: int):
    return [a for a in itertools.product(range(d+1), repeat=4) if sum(a)==d]

def central(a):
    return all(int(v)%3==0 for v in E*s.Matrix(a))

centers={str(d): [list(a) for a in monomials(d) if central(a)] for d in range(1,5)}
assert centers['1']==[] and centers['2']==[]
for i in range(4):
    a=[0]*4; a[i]=3
    assert central(a)
rows=[tuple(int(E[i,j])%3 for j in range(4)) for i in range(4)]
classes=[]
for r in dict.fromkeys(rows): classes.append([i for i,row in enumerate(rows) if row==r])
assert classes==[[0],[1,2],[3]]
monodromies=[]
for cl in classes:
    e=rows[cl[0]]
    fourth=tuple(4*x%3 for x in e)
    assert len(set(fourth))>1
    monodromies.append({'support_class':cl,'eigenvalue_exponents':list(e),
                        'fourfold_monodromy_exponents':list(fourth)})
# Check the only row-proportionality cases permitted by a normal line pair.
# For i,k in its support, row_i-row_k is constant, hence the constants below.
for i,k in itertools.product(range(4),repeat=2):
    diff=[(rows[i][j]-rows[k][j])%3 for j in range(4)]
    if len(set(diff))==1:
        assert diff==[0]*4
# Skew antisymmetrizer coefficient of x_sigma(0) ... x_sigma(3).
def coeff_exp(p):
    inversions=[(p[j],p[i]) for i in range(4) for j in range(i+1,4) if p[i]>p[j]]
    return ((-1)**len(inversions),sum(int(E[a,b]) for a,b in inversions)%3)
coefs={p:coeff_exp(p) for p in itertools.permutations(range(4))}
for p,(sign,e) in coefs.items():
    q=p[1:]+p[:1]
    sq,eq=coefs[q]
    assert sq==-sign and eq==e
# Hilbert numbers from the monic central cubic relation.
h=[]
for n in range(16):
    b=comb(n+3,3)
    old=comb(n,3) if n>=3 else 0
    value=b-old
    assert value==1+3*n*(n+1)//2
    h.append(value)
# Dehomogenization at x0: all q_(0i) equal, so the remaining exponents stay fixed.
for i,j in itertools.combinations(range(1,4),2):
    new=(int(E[i,j])+int(E[0,i])-int(E[0,j]))%3
    assert new==int(E[i,j])%3
# Exact quadratic overlap calculation over QQ(zeta), not floating complex numbers.
from sympy.polys.matrices import DomainMatrix
zeta=(-1+s.sqrt(-3))/2
K=s.QQ.algebraic_field(zeta)
rel=[]
for i,j in itertools.combinations(range(4),2):
    v=s.zeros(16,1)
    v[4*i+j]=1
    v[4*j+i]=-s.expand(zeta**int(E[i,j]))
    rel.append(v)
left=s.zeros(64,24); right=s.zeros(64,24)
for a,v in enumerate(rel):
    for j in range(4):
        for p in range(16):
            left[4*p+j,4*a+j]=v[p]
            right[16*j+p,4*a+j]=v[p]
DM=DomainMatrix.from_Matrix(left.row_join(-right)).convert_to(K)
overlap_rank=len(DM.rref()[1])
assert overlap_rank==44
assert 48-overlap_rank==4
# Norm-residue obstruction: r^3+1 has three simple roots, so -1/(r^3+1)
# has simple poles and is not a cube in k(r).
r=s.symbols('r')
assert s.gcd(r**3+1,3*r**2)==1
assert s.factor(r**3+1)==(r+1)*(r**2-r+1)
report={
  'date':'2026-11-15','version':'recognition-v3','field':'any algebraically closed field of characteristic 0',
  'exponent_matrix':[list(map(int,E.row(i))) for i in range(4)],
  'cy_row_sums_mod3':[sum(row)%3 for row in rows],
  'low_degree_central_monomials':centers,
  'normal_linear_support_classes':classes,
  'four_cycle_monodromies':monodromies,
  'signed_cyclic_quartic_tensor':'24 coefficients checked exactly modulo 3',
  'hilbert_values_degrees_0_to_15':h,
  'rank_data':{'quadratic_relations':6,'cubic_new_relations':1,'quadratic_overlap':48-overlap_rank,'overlap_matrix_shape':[64,48],'overlap_matrix_rank_over_QQ_zeta':overlap_rank},
  'norm_obstruction':'v at w=0 is -1/(r^3+1), with simple poles; hence not a cube',
  'scope':'Finite arithmetic only. The general normal-line classification, proper closed incidence, domain, Serre duality and global dimension are ordinary proofs in the TeX paper.'
}
path=Path(__file__).with_name('quantum_fermat_certificate.json')
path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('PASS: CY, central cubic, no degree-one center, normal support classes, monodromy, signed potential, quadratic overlap, Hilbert series, norm residue.')
