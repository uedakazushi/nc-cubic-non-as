#!/usr/bin/env python3
"""Exact path certificates for the 80 five-zero support orbits.

Each of the 71 witnesses is over the characteristic-zero fraction field.
For a fixed start, recursive_dims retains the quotient basis and right-arrow
maps. It gives the same quotient as the full path matrix. The nine residual
orbits are only necessary-condition survivors, not verified helices.
"""
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import os
import subprocess
from sympy import GF
from sympy.polys.matrices import DomainMatrix
from verify_ou_zero_supports import I, PAIRS, TRIPLES, normalize, support, path_dim


def recursive_dims(c,K,r,i,nmax):
    old=[int(g==i) for g in I]
    dim=[1]*3
    M=[[([[K.one]] if g==i else []) for h in I] for g in I]
    yield 1,dim
    for n in range(2,nmax+1):
        offset=[0,dim[0],dim[0]+dim[1]]; D=sum(dim)
        dims=[]; newmaps=[[None]*3 for _ in I]; rot=(r+n-2)%3
        for j in I:
            rows=[]
            for g in I:
                cs=[c[(g,h,j)[-rot:]+(g,h,j)[:-rot]] for h in I]
                for b in range(old[g]):
                    row={}
                    for h in I:
                        if cs[h]:
                            for k,val in enumerate(M[g][h][b]):
                                if val:row[offset[h]+k]=cs[h]*val
                    if row:rows.append(row)
            T=DomainMatrix(dict(enumerate(rows)),(len(rows),D),K)
            R,piv=T.rref(); RR=R.to_dok()
            free=[b for b in range(D) if b not in piv]
            where={p:a for a,p in enumerate(piv)}
            dims.append(len(free)); allimages=[]
            for b in range(D):
                if b in where:
                    allimages.append([-RR.get((where[b],k),K.zero) for k in free])
                else:
                    allimages.append([K.one if k==b else K.zero for k in free])
            for h in I:newmaps[h][j]=allimages[offset[h]:offset[h]+dim[h]]
        yield n,dims
        old,dim,M=dim,dims,newmaps


def expected(n,i,j):
    q,s=divmod(n,3)
    if s==0:return 3*q*(q+1)//2+int(i==j)
    if s==1:return (3*q*q+5*q+2)//2
    return (3*q*q+7*q+4)//2

# mask, length, rotation, start, end, dimension, expected, recursive flag
WITNESSES = (
    (129990239, 7, 1, 2, 2, 13, 12, 0),
    (130005631, 8, 1, 2, 0, 17, 15, 1),
    (130006647, 8, 1, 2, 0, 17, 15, 0),
    (130014815, 5, 2, 1, 1, 9, 7, 0),
    (132087391, 7, 1, 1, 2, 13, 12, 0),
    (132102783, 6, 2, 1, 1, 12, 10, 1),
    (132103295, 5, 1, 2, 0, 8, 7, 1),
    (132111967, 7, 1, 1, 2, 13, 12, 0),
    (133103199, 6, 2, 1, 2, 10, 9, 0),
    (133135967, 7, 1, 2, 2, 13, 12, 0),
    (133135991, 8, 1, 2, 0, 17, 15, 0),
    (133144191, 6, 2, 1, 1, 11, 10, 0),
    (133148287, 6, 0, 0, 1, 10, 9, 0),
    (133151871, 8, 1, 0, 0, 17, 15, 1),
    (133152375, 8, 1, 2, 0, 17, 15, 0),
    (133164663, 4, 2, 2, 0, 6, 5, 1),
    (133660255, 6, 1, 2, 2, 12, 10, 0),
    (133672575, 5, 1, 1, 1, 8, 7, 0),
    (133676159, 8, 1, 0, 0, 17, 15, 1),
    (133676663, 8, 1, 2, 0, 17, 15, 1),
    (133688927, 5, 2, 1, 1, 8, 7, 0),
    (133690999, 8, 1, 2, 0, 17, 15, 0),
    (133824095, 7, 0, 0, 0, 13, 12, 0),
    (133889631, 4, 0, 0, 2, 6, 5, 0),
    (133922399, 7, 0, 0, 2, 13, 12, 0),
    (133922423, 8, 1, 0, 0, 17, 15, 0),
    (133930623, 4, 0, 2, 0, 6, 5, 0),
    (133938807, 6, 2, 2, 2, 12, 10, 1),
    (133946975, 5, 2, 1, 1, 8, 7, 0),
    (133951071, 5, 1, 2, 0, 8, 7, 0),
    (133951095, 4, 2, 0, 0, 6, 5, 1),
    (133953143, 8, 0, 0, 2, 16, 15, 0),
    (134078047, 7, 0, 0, 2, 13, 12, 0),
    (134082143, 7, 0, 0, 2, 13, 12, 0),
    (134085727, 6, 0, 0, 0, 12, 10, 0),
    (134086223, 4, 0, 0, 1, 6, 5, 0),
    (134086231, 5, 0, 0, 2, 8, 7, 0),
    (134086238, 7, 0, 0, 2, 13, 12, 0),
    (134135391, 6, 0, 0, 0, 11, 10, 0),
    (134147679, 7, 0, 0, 1, 13, 12, 0),
    (134149727, 7, 0, 0, 1, 13, 12, 0),
    (134150751, 5, 0, 1, 2, 8, 7, 0),
    (134151263, 6, 0, 0, 1, 10, 9, 0),
    (134151767, 5, 0, 0, 2, 8, 7, 0),
    (134151774, 6, 0, 0, 0, 11, 10, 0),
    (134168183, 8, 1, 2, 0, 17, 15, 0),
    (134176351, 7, 1, 2, 2, 13, 12, 0),
    (134180447, 5, 1, 2, 0, 8, 7, 0),
    (134182495, 7, 1, 2, 2, 13, 12, 0),
    (134182519, 8, 1, 2, 0, 17, 15, 0),
    (134183519, 7, 1, 2, 1, 13, 12, 0),
    (134184031, 5, 1, 2, 0, 8, 7, 0),
    (134184055, 5, 1, 2, 0, 8, 7, 0),
    (134184527, 4, 0, 0, 1, 6, 5, 0),
    (134184535, 5, 0, 0, 2, 8, 7, 0),
    (134184542, 7, 1, 2, 2, 13, 12, 0),
    (134192759, 6, 2, 1, 1, 11, 10, 0),
    (134198903, 5, 0, 0, 2, 8, 7, 0),
    (134199927, 8, 1, 2, 0, 17, 15, 0),
    (134200439, 8, 1, 0, 0, 17, 15, 0),
    (134205023, 5, 1, 1, 0, 8, 7, 0),
    (134208607, 7, 0, 0, 2, 13, 12, 0),
    (134209118, 7, 1, 2, 2, 13, 12, 0),
    (134211191, 5, 1, 1, 0, 8, 7, 1),
    (134212191, 7, 1, 1, 1, 13, 12, 0),
    (134213214, 7, 1, 1, 2, 13, 12, 0),
    (134214263, 5, 0, 0, 1, 8, 7, 1),
    (134214775, 6, 0, 0, 0, 11, 10, 1),
    (134215255, 5, 0, 0, 2, 8, 7, 0),
    (134216279, 5, 0, 0, 1, 8, 7, 0),
    (134216798, 5, 0, 0, 0, 8, 7, 0),
)
RESIDUAL = [117374702, 125763310, 132054255, 132054766, 132086511, 133099247, 133873391, 134135534, 134147822]


def sample(mask,free,values=(5,11,17),prime=1009):
    F=GF(prime);T=set(support(mask))
    return {a:F.convert(values[free.index(a)] if a in free else int(a in T)) for a in TRIPLES},F


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--first',type=int,default=0)
    p.add_argument('--last',type=int,default=len(WITNESSES))
    args=p.parse_args()
    assert 0<=args.first<=args.last<=len(WITNESSES)
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'support_orbits';src=Path(__file__).with_name('ou_zero_supports.cpp')
        subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
        masks=[int(m) for z,m in map(str.split,subprocess.check_output([str(exe)],text=True).splitlines()) if z=='5']
    assert len(masks)==80 and len(WITNESSES)==71 and len(RESIDUAL)==9
    assert sorted(masks)==sorted([row[0] for row in WITNESSES]+RESIDUAL)
    assert Counter(row[1] for row in WITNESSES)=={4:6,5:20,6:12,7:19,8:14}
    # Cross-check both path computations in a finite field.
    for mask in (117374702,125763310,134199927):
        _,_,free=normalize(mask);c,F=sample(mask,free)
        for r,i in PAIRS:
            for n,d in recursive_dims(c,F,r,i,5):
                if n>=2:assert d==[path_dim(c,n,r,i,j,F) for j in I]
    for mask,n,r,i,j,d,e,rec in WITNESSES[args.first:args.last]:
        c,K,free=normalize(mask);assert len(free)==3
        actual=dict(recursive_dims(c,K,r,i,n))[n][j] if rec else path_dim(c,n,r,i,j,K)
        assert actual==d and d>e==expected(n,i,j)
        print('PASS witness',mask,n,r,i,j,d,e,flush=True)
    for mask in RESIDUAL:
        _,_,free=normalize(mask);assert len(free)==3
        c,F=sample(mask,free)
        for r,i in PAIRS:
            for n,d in recursive_dims(c,F,r,i,12):
                assert all(d[j]==expected(n,i,j) for j in I)
        print('PASS necessary dimensions to length 12 (finite-field sample)',mask,flush=True)
    print('PASS requested exact certificates',args.first,args.last,flush=True)

if __name__=='__main__':
    main()
