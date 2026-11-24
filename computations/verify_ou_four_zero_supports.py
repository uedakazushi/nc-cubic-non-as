#!/usr/bin/env python3
"""Four-zero support certificates and smooth-cut samples over characteristic zero.
The monomial and primitive polynomial elimination uses exact arithmetic.
"""
from collections import Counter, defaultdict
from heapq import heappop, heappush
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import os
import subprocess
from time import monotonic
from sympy import Matrix, QQ, symbols, groebner, cancel
from verify_ou_zero_supports import normalize, support, path_dim, I, PAIRS, TRIPLES
from verify_ou_seven_zero_curves import determinant_matrix
from verify_ou_five_zero_supports import recursive_dims, expected

def rows_for(c,n,r,i,j,K):
 words=list(product(I,repeat=n-1));idx={w:b for b,w in enumerate(words)};rows=[]
 for pos in range(n-1):
  rot=(r+pos)%3
  for w in words:
   if w[pos]:continue
   row={}
   for a in I:
    v=w[:pos]+(a,)+w[pos+1:];f=(i,)+v+(j,);t=f[pos:pos+3];q=c[t[-rot:]+t[:-rot]]
    if q:row[idx[v]]=q
   if row:rows.append(row)
 return len(words),rows


def laurent_rank(rows,K):
 from fractions import Fraction
 d=len(K.symbols)
 def from_field(p):
  nn,dd=p.numer.to_dict(),p.denom.to_dict();assert len(dd)==1
  e,b=next(iter(dd.items()));return {tuple(x-y for x,y in zip(a,e)):Fraction(int(c.numerator),int(c.denominator))/Fraction(int(b.numerator),int(b.denominator)) for a,c in nn.items()}
 rows=[{k:from_field(q) for k,q in r.items()} for r in rows]
 cols=defaultdict(set)
 for i,r in enumerate(rows):
  for j in r:cols[j].add(i)
 versions=[0]*len(rows);heap=[]
 def push(i):
  r=rows[i];mono=[j for j,p in r.items() if len(p)==1]
  if mono:
   j=min(mono,key=lambda j:(len(cols[j]),j))
   heappush(heap,((len(r)-1)*(len(cols[j])-1),len(r),i,j,versions[i]))
 for i in range(len(rows)):push(i)
 rank=0;tt=monotonic()
 while heap:
  score,n,i,j,v=heappop(heap)
  if versions[i]!=v or j not in rows[i] or len(rows[i][j])!=1:continue
  r=rows[i]
  # Reprioritize if col count changed significantly.
  mono=[k for k,p in r.items() if len(p)==1];jj=min(mono,key=lambda k:(len(cols[k]),k));ss=(len(r)-1)*(len(cols[jj])-1)
  if ss!=score or jj!=j:
   heappush(heap,(ss,len(r),i,jj,versions[i]));continue
  e,c=next(iter(r[j].items()));norm={k:{tuple(x-y for x,y in zip(f,e)):b/c for f,b in p.items()} for k,p in r.items() if k!=j}
  for k in r:cols[k].discard(i)
  rows[i]={};versions[i]+=1
  for ii in list(cols[j]):
   row=rows[ii];x=row.pop(j);cols[j].discard(ii)
   for k,w in norm.items():
    old=row.get(k,{});new=old.copy()
    for f,b in x.items():
     for g,c in w.items():
      z=tuple(a+b for a,b in zip(f,g));q=new.get(z,0)-b*c
      if q:new[z]=q
      else:new.pop(z,None)
    if new:
     row[k]=new
     if not old:cols[k].add(ii)
    elif old:del row[k];cols[k].discard(ii)
   versions[ii]+=1;push(ii)
  rank+=1
 # Finish exact over fraction field. Strip Laurent denominators rowwise.
 rem=[r for r in rows if r];cc=sorted(set(k for r in rem for k in r));idx={j:i for i,j in enumerate(cc)}
 print('monomial pivots',rank,'remainder',len(rem),len(cc),'time',round(monotonic()-tt,2),flush=True)
 # Clear Laurent monomials rowwise and remove duplicate rows.
 R=K.field.ring
 rr=[];seen=set()
 for r in rem:
  e=tuple(min(f[j] for p in r.values() for f in p) for j in range(d));f0,c=next(iter(sorted(r[min(r)].items())))
  row=tuple((idx[j],tuple(sorted((tuple(x-y for x,y in zip(f,e)),b/c) for f,b in p.items()))) for j,p in sorted(r.items()))
  if row in seen:continue
  seen.add(row);rr.append({j:R.from_dict({f:QQ.convert(b) for f,b in p}) for j,p in row})
 print('dedup',len(rr),flush=True)
 # Repeated row-wise primitive normalization and exact forward elimination.
 def primitive(row):
  if not row:return row
  ps=sorted(row.values(),key=len);g=ps[0]
  for p in ps[1:]:
   if g.is_ground:break
   g=g.gcd(p)
  if not g.is_ground:row={j:p.exquo(g) for j,p in row.items()}
  c=next(iter(row.values())).LC
  if c!=QQ.one:row={j:p.quo_ground(c) for j,p in row.items()}
  return row
 rr=[primitive(r) for r in rr]; rr=list({tuple(sorted((j,tuple(sorted(p.items()))) for j,p in r.items())):r for r in rr if r}.values());print('primitive',len(rr),'monomials',sum(any(len(p)==1 for p in r.values()) for r in rr),flush=True)
 pcount=0
 while rr:
  ccnt=Counter(j for r in rr for j in r)
  i,j=min(((i,j) for i,r in enumerate(rr) for j in r),key=lambda ij:(len(rr[ij[0]][ij[1]])!=1,(len(rr[ij[0]])-1)*(ccnt[ij[1]]-1),len(rr[ij[0]][ij[1]]),sum(len(p) for p in rr[ij[0]].values())))
  pivot=rr.pop(i);v=pivot[j];new=[]
  for row in rr:
   if j not in row:new.append(row);continue
   x=row[j];g=v.gcd(x);a=v.exquo(g);b=x.exquo(g)
   r={k:a*row.get(k,R.zero)-b*pivot.get(k,R.zero) for k in row.keys()|pivot.keys() if k!=j};r={k:p for k,p in r.items() if p}
   if r:new.append(primitive(r))
  rr=list({tuple(sorted((j,tuple(sorted(p.items()))) for j,p in r.items())):r for r in new}.values());pcount+=1
 return rank+pcount


WITNESSES = (
    (134217295, 4, 0, 0, 1, 6, 'full'),
    (134217303, 5, 0, 0, 2, 8, 'full'),
    (134217310, 7, 1, 2, 2, 13, 'poly'),
    (134086239, 7, 0, 0, 2, 13, 'full'),
    (134151775, 7, 0, 0, 1, 13, 'full'),
    (134184543, 7, 1, 2, 2, 13, 'poly'),
    (134209119, 7, 1, 2, 2, 13, 'poly'),
    (134213215, 7, 1, 1, 2, 13, 'full'),
    (134216799, 7, 0, 0, 2, 13, 'full'),
    (134184567, 8, 1, 2, 0, 17, 'poly'),
    (134200951, 8, 1, 2, 0, 17, 'poly'),
    (134213239, 5, 1, 1, 0, 8, 'full'),
    (134215287, 8, 1, 2, 0, 17, 'poly'),
    (134216823, 8, 1, 0, 0, 17, 'rec'),
    (134192767, 6, 2, 1, 1, 11, 'full'),
    (134196863, 6, 0, 0, 1, 10, 'full'),
    (134199935, 8, 1, 2, 0, 17, 'poly'),
    (134200447, 8, 1, 0, 0, 17, 'rec'),
    (133152383, 8, 1, 2, 0, 17, 'poly'),
    (133676671, 8, 1, 2, 0, 17, 'poly'),
    (133938815, 6, 2, 2, 2, 12, 'full'),
    (133951103, 5, 1, 2, 0, 8, 'full'),
)
RESIDUAL = [133136111, 133660399, 133889775, 134135535, 134147823, 134151407, 134151918, 134183663]
SMOOTH = [133660399, 133889775, 134135535, 134147823, 134151407, 134151918, 134183663]


def check_smooth():
    x,y,z=symbols('x y z')
    xyz=Matrix([x,y,z])
    values=(2,3,5,7)
    for mask in SMOOTH:
        c,K,free=normalize(mask)
        assert len(free)==4
        c={a:K.to_sympy(v).subs(dict(zip(K.symbols,values))) for a,v in c.items()}
        L=determinant_matrix(c)
        assert L.rank()==3
        p=Matrix.hstack(*L.nullspace())*xyz
        f=cancel(p[0]*p[3]*p[4]-p[1]*p[2]*p[5])
        for a in xyz:
            other=[b for b in xyz if b!=a]
            assert list(groebner([f.diff(b).subs(a,1) for b in xyz],*other))==[1]
        print('PASS smooth sample',mask,flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--first',type=int,default=0)
    parser.add_argument('--last',type=int,default=len(WITNESSES))
    args=parser.parse_args()
    assert 0<=args.first<=args.last<=len(WITNESSES)
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'orbits';src=Path(__file__).with_name('ou_zero_supports.cpp')
        subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
        masks=[int(m) for z,m in map(str.split,subprocess.check_output([str(exe)],text=True).splitlines()) if z=='4']
    assert len(masks)==30 and sorted(masks)==sorted([r[0] for r in WITNESSES]+RESIDUAL)
    assert Counter(r[1] for r in WITNESSES)=={4:1,5:3,6:3,7:7,8:8}
    for m in masks:assert len(normalize(m)[2])==4
    for m,n,r,i,j,d,method in WITNESSES[args.first:args.last]:
        c,K,free=normalize(m)
        if method=='poly':
            size,rows=rows_for(c,n,r,i,j,K)
            actual=size-laurent_rank(rows,K)
        elif method=='rec':actual=dict(recursive_dims(c,K,r,i,n))[n][j]
        else:actual=path_dim(c,n,r,i,j,K)
        assert actual==d and d>expected(n,i,j), (m,actual,d)
        print('PASS exact',m,n,r,i,j,d,expected(n,i,j),flush=True)
    check_smooth()

if __name__=='__main__':
    main()
