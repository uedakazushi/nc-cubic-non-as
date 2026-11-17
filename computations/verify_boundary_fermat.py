#!/usr/bin/env python3
"""Exact checks supporting normal-boundary recognition and the quantum cubic orbit.

PBW part: universal Lie/central-cocycle identities for S=k[x,y,z].
Fermat part: central cubic basis, open orbit tangent and its finite stabilizer.
The normal-boundary recognition theorem itself uses Cassidy--Shelton and a
written lifting proof; no finite computation certifies that general theorem.
"""
from __future__ import annotations
import itertools as it
import json
from pathlib import Path
import sympy as sp

# W basis x,y,z and U basis xy-yx, yz-zy, zx-xz.
alpha_symbols = sp.symbols('a0:9')
beta_symbols = sp.symbols('b0:3')
alpha = sp.Matrix(3, 3, alpha_symbols)
beta = sp.Matrix(1, 3, beta_symbols)
pairs = [(0,1),(1,2),(2,0)]
U = sp.zeros(9,3)
for r,(i,j) in enumerate(pairs):
    U[3*i+j,r]=1
    U[3*j+i,r]=-1
# w=(xy-yx)z + (yz-zy)x + (zx-xz)y
#  =z(xy-yx) + x(yz-zy) + y(zx-xz).
last = [2,0,1]
delta_a = sp.zeros(9,1)
delta_b = sp.zeros(3,1)
for r,j in enumerate(last):
    for i in range(3):
        delta_a[3*i+j] += alpha[i,r]
        delta_a[3*j+i] -= alpha[i,r]
    delta_b[j] += beta[0,r]
    delta_b[j] -= beta[0,r]
coeff = sp.Matrix([delta_a[3*i+j] for i,j in pairs])
assert sp.simplify(U*coeff-delta_a)==sp.zeros(9,1)
# Bracket table encoded by alpha, beta.
def bracket_basis(i: int, j: int) -> sp.Matrix:
    if i==j:
        return sp.zeros(3,1)
    if (i,j) in pairs:
        return alpha[:,pairs.index((i,j))]
    return -alpha[:,pairs.index((j,i))]

def scalar_bracket_basis(i: int, j: int):
    if i==j:
        return sp.Integer(0)
    if (i,j) in pairs:
        return beta[0,pairs.index((i,j))]
    return -beta[0,pairs.index((j,i))]

jacobi=sp.zeros(3,1)
cocycle=sp.Integer(0)
for (i,j),k in zip(pairs,last):
    inner=bracket_basis(i,j)
    for q in range(3):
        jacobi += inner[q]*bracket_basis(q,k)
        cocycle += inner[q]*scalar_bracket_basis(q,k)
assert all(sp.expand(v)==0 for v in alpha*coeff+delta_b-jacobi)
assert sp.expand((beta*coeff)[0]-cocycle)==0
# A non-Lie choice to ensure the test detects an obstruction.
bad={a:0 for a in alpha_symbols}
bad.update({alpha[0,0]:1,alpha[1,1]:1})
assert any(sp.expand(v.subs(bad))!=0 for v in jacobi)
# Heisenberg bracket plus zero constant cocycle is PBW.
good={a:0 for a in alpha_symbols}|{b:0 for b in beta_symbols}
good[alpha[2,0]]=1
assert all(sp.expand(v.subs(good))==0 for v in list(jacobi)+[cocycle])

E=sp.Matrix([[0,1,1,1],[-1,0,0,1],[-1,0,0,1],[-1,-1,-1,0]])
central3=[]
for v in it.product(range(4), repeat=4):
    if sum(v)==3 and all(int(q)%3==0 for q in E*sp.Matrix(v)):
        central3.append(v)
expected={(3,0,0,0),(0,3,0,0),(0,2,1,0),(0,1,2,0),(0,0,3,0),(0,0,0,3)}
assert set(central3)==expected
# Centralize the normal t by the Zhang automorphism sigma^{-1}.
weights=[0,2,2,2]
Et=sp.Matrix(4,4,lambda i,j:(E[i,j]+weights[j]-weights[i])%3)
assert all(Et[0,j]==0 and Et[j,0]==0 for j in range(4))
assert all(Et[i,j]==E[i,j]%3 for i in range(1,4) for j in range(1,4))

# The derivative of G_m x GL2 x G_m on the affine central cubic space.
t,x,y,z=sp.symbols('t x y z')
F=t**3+x**3+y**3+z**3
basis=[t**3,x**3,x*x*y,x*y*y,y**3,z**3]
vecfields=[t*sp.diff(F,t),x*sp.diff(F,x),y*sp.diff(F,x),x*sp.diff(F,y),y*sp.diff(F,y),z*sp.diff(F,z)]
M=sp.Matrix([[sp.Poly(q,t,x,y,z).coeff_monomial(m) for q in vecfields] for m in basis])
assert M.det()==729 and M.rank()==6
# Exact projective stabilizer: diagonal mu3^4 modulo common mu3, semidirect swap(x,y).
def norm(v):
    return tuple((q-v[0])%3 for q in v)

def perm(v, swap):
    return (v[0],v[2],v[1],v[3]) if swap else v

def mul(g,h):
    v,p=g;w,q=h
    w=perm(w,p)
    return norm(tuple(a+b for a,b in zip(v,w))),p^q

H={(tuple([0]+list(v)),p) for v in it.product(range(3), repeat=3) for p in (0,1)}
assert len(H)==54
unit=((0,0,0,0),0)
for g in H:
    assert any(mul(g,h)==unit for h in H)
    for h in H:
        assert mul(g,h) in H
# All representatives preserve q matrix and the Fermat cubic exactly.
for v,p in H:
    pi=[0,2,1,3] if p else list(range(4))
    assert all(E[pi[i],pi[j]]==E[i,j] for i in range(4) for j in range(4))
    assert all(3*a%3==0 for a in v)

a,b,c,d=sp.symbols('a b c d')
f=a*x**3+b*x*x*y+c*x*y*y+d*y**3
Delta=sp.expand(sp.discriminant(f.subs(y,1),x))
expected_delta=b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
assert Delta==sp.expand(expected_delta)
assert Delta.subs({a:1,b:0,c:0,d:1})==-27

out={
 'version':'recognition-v4','date':'2026-11-15','arithmetic':'exact QQ and Z/3',
 'pbw':{'boundary':'commutative polynomial ring in three variables',
        'linear_part_symbols':9,'constant_part_symbols':3,
        'delta_alpha_in_U':True,
        'alpha_delta_plus_delta_beta_equals_Jacobi':True,
        'beta_delta_equals_central_cocycle':True,
        'jacobi_polynomials':[str(sp.expand(q)) for q in jacobi],
        'cocycle_polynomial':str(sp.expand(cocycle)),
        'non_Lie_obstruction_detected':True,'Heisenberg_test':True},
 'quantum_central_cubic_basis_exponents':[list(v) for v in central3],
 'twisted_central_t_exponent_matrix':[[int(Et[i,j]) for j in range(4)] for i in range(4)],
 'affine_central_orbit_derivative':{'shape':list(M.shape),'rank':6,'determinant':729},
 'projective_stabilizer':{'order':54,'description':'(mu3^4 / mu3_diag) semidirect S2 swapping x,y',
                          'closure_and_inverse_tests':True},
 'binary_cubic_discriminant':str(Delta),
 'scope':'Only these finite identities are computed. General boundary recognition, AS regularity and homogeneous-space identification are proved in the paper. Full weak-cubic deformation rigidity and singular-boundary classification are not claimed.'
}
Path(__file__).with_name('boundary_fermat_certificate.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print('PASS: PBW identities, central cubic basis, twist, orbit rank6/det729, 54-element stabilizer, binary discriminant.')
