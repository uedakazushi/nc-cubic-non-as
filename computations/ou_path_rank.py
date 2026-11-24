from collections import Counter
from itertools import product
from sympy import QQ,symbols
from time import monotonic
I=range(3); TRIPLES=list(product(I,repeat=3)); CUBE=list(product((1,2),repeat=3))
def rotate(c,n):return {t:c[t[-n:]+t[:-n]] for t in TRIPLES}
def path_matrix(c,length,start,end,K):
    words=list(product(I,repeat=length-1));idx={w:i for i,w in enumerate(words)};rows=[]
    for pos in range(length-1):
        r=rotate(c,pos%3)
        for w in words:
            if w[pos]:continue
            row={}
            for a in I:
                v=w[:pos]+(a,)+w[pos+1:];f=(start,)+v+(end,);q=r[f[pos:pos+3]]
                if q:row[idx[v]]=q
            rows.append(row)
    return len(words),rows
def exact_rank(rows,K):
    cols={}
    for i,r in enumerate(rows):
        for j in r:cols.setdefault(j,set()).add(i)
    rank=0;t=monotonic()
    while cols:
        best=None
        for i,r in enumerate(rows):
            if not r:continue
            for j,v in r.items():
                score=(v not in (K.one,-K.one),(len(r)-1)*(len(cols[j])-1),i,j)
                if best is None or score<best:best=score
            if best[:2]==(False,0):break
        _,_,i,j=best;r=rows[i];v=r[j]
        norm={k:w/v for k,w in r.items() if k!=j}
        for k in r:
            cols[k].remove(i)
            if not cols[k]:del cols[k]
        rows[i]={}
        for ii in list(cols.get(j,())):
            row=rows[ii];x=row.pop(j)
            for k,w in norm.items():
                old=row.get(k,K.zero);new=old-x*w
                if new:
                    row[k]=new
                    if not old:cols.setdefault(k,set()).add(ii)
                elif old:
                    del row[k];cols[k].remove(ii)
                    if not cols[k]:del cols[k]
        cols.pop(j,None);rank+=1
        if rank%100==0:print('rank',rank,'time',round(monotonic()-t,1),flush=True)
    return rank
