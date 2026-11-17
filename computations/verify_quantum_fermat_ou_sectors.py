#!/usr/bin/env python3
"""Exact finite group and twisted-sector character checks for quantum Fermat.
H=(Z/3)^2 acts on the Fermat cubic X via weights
(t,x,y,z)=(a,b,b,-a+b). Discrete torsion has nondegenerate
commutator epsilon_g(h)=zeta**(h_a*g_b-g_a*h_b).
"""
from itertools import product
import json
from pathlib import Path
G=list(product(range(3),repeat=2))

def weights(g):
    a,b=g
    return [a%3,b%3,b%3,(-a+b)%3]

def torsion(g,h):
    a,b=g; c,d=h
    return (c*b-a*d)%3

def fixed(g):
    if g==(0,0):return ('surface',9)
    counts={k:weights(g).count(k) for k in set(weights(g))}
    # X intersect disjoint projective eigenspaces; each projective coordinate eigenspace of dimension r-1 meets Fermat cubic.
    # coordinate subspace with >=3 variables: Fermat plane cubic elliptic, Euler 0.
    # 2 variables: 3 points. single variable: no points on Fermat X.
    if 3 in counts.values():return ('elliptic',0)
    return ('points',3*sum(n==2 for n in counts.values()))

# Cocycle commutators on the (Z/3)^2 grading reproduce the quantum Fermat q matrix.
charweights=[(1,0),(0,1),(0,1),(2,1)]  # t,x,y,z
E=[[0,1,1,1],[2,0,0,1],[2,0,0,1],[2,2,2,0]]
for i in range(4):
    for j in range(4):
        a,b=charweights[i];c,d=charweights[j]
        assert (a*d-c*b)%3==E[i][j]

elliptic_t=[g for g in G if g!=(0,0) and weights(g)[1]==weights(g)[2]==weights(g)[3]]
elliptic_z=[g for g in G if g!=(0,0) and weights(g)[0]==weights(g)[1]==weights(g)[2]]
assert set(elliptic_t)=={(0,1),(0,2)}
assert set(elliptic_z)=={(1,1),(2,2)}
# H action on elliptic divisors: E_t x^3+y^3+z^3=0, z maps zeta^-a; on E_z t^3+x^3+y^3=0,
# t maps zeta^(a-b), hence H0(Omega) characters -a and +(a-b).
# HH_-1 receives H1(O_E): conjugate character; HH_+1 receives H0(Omega).
counts={-1:0,0:0,1:0}
sector=[]
# untwisted HH_0 invariants: H0(O)=1,H2(Omega^2)=1,H1(Omega^1)=1 from trace averaging on H^2(X).
traces={g:fixed(g)[1]-2 for g in G}
assert traces[(0,0)]==7
assert sum(traces.values())//9==1
counts[0]+=3
for g in G:
    if g==(0,0):continue
    kind,chi=fixed(g)
    record={'g':list(g),'weights':weights(g),'fixed':kind,'euler_characteristic':chi}
    if kind=='elliptic':
        assert g in elliptic_t+elliptic_z
        if g in elliptic_t:
            chOmega=lambda h:(-h[0])%3
        else:
            chOmega=lambda h:(h[0]-h[1])%3
        chO=lambda h:(-chOmega(h))%3
        for hhdeg,ch in ((1,chOmega),(-1,chO)):
            if all((ch(h)+torsion(g,h))%3==0 for h in G):
                counts[hhdeg]+=1
                record['HH_'+str(hhdeg)]=1
    else:
        # Fixed points on x,y eigenspace are individually H-fixed, so a nontrivial torsion character gives 0.
        # Fixed points on t,z eigenspace form H/K with |K|=3 and contribute 1 iff torsion char trivial on K.
        w=weights(g)
        tz=(w[0]==w[3])
        if tz:
            K=[h for h in G if weights(h)[0]==weights(h)[3]]
            assert len(K)==3
            assert all(torsion(g,h)==0 for h in K)
            counts[0]+=1
            record['HH_0']=1
    sector.append(record)
assert counts=={-1:2,0:5,1:2},counts
assert sum(v for v in counts.values())==9
print('PASS: elliptic sectors E_t,E_z each yield HH_-1=HH_1=1; HH_0=5.');print(counts)
Path(__file__).with_name('quantum_fermat_ou_sector_certificate.json').write_text(json.dumps({'group':'(Z/3)^2','weights':{'t':'a','x':'b','y':'b','z':'-a+b'},'torsion_exponent':'h_a*g_b-g_a*h_b','sectors':sector,'hh_dimensions':counts,'trace_on_invariant_picard':sum(traces.values())//9},indent=2,ensure_ascii=False)+'\n')