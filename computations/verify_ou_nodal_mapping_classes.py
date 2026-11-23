#!/usr/bin/env python3
"""Exact free-group verification for normalized nodal-cycle line-bundle patterns.

Generators (a,b_1,...,b_c) are the horizontal and based vertical loops.
H and V_i are representatives of the induced outer automorphisms of the
standard horizontal and vertical Dehn twists (see the proof in the paper).
No assertion of realizability of the multidegree patterns is made.
"""
from collections import Counter
from itertools import product, permutations
from verify_ou_reduced_multidegrees import patterns, invariants, mul, normalize


def inv(word):
    return tuple(-x for x in reversed(word))


def red(*words):
    out = []
    for word in words:
        for x in word:
            if out and out[-1] == -x:
                out.pop()
            else:
                out.append(x)
    return tuple(out)


def apply(F, word):
    return red(*(F[x-1] if x > 0 else inv(F[-x-1]) for x in word))


def compose(F, G):
    return tuple(apply(F, word) for word in G)


def identity(c):
    return tuple((i,) for i in range(1, c+2))


def power_word(i, e):
    return ((i,) if e >= 0 else (-i,))*abs(e)


def H(c, e=1):
    return ((1,),) + tuple((j,)+power_word(1,-e) for j in range(2,c+2))


def V(c, i, e=1):
    j = i+2
    b = power_word(j,e)
    return (b+(1,),) + tuple((k,) if k <= j else b+(k,)+inv(b)
                             for k in range(2,c+2))


def product_of(maps, c):
    F = identity(c)
    for G in maps:
        F = compose(F,G)
    return F


def vertical(c, d):
    return product_of((V(c,i,e) for i,e in enumerate(d)),c)


def line_twist(c, d):
    return compose(compose(vertical(c,d),H(c)),vertical(c,mul(-1,d)))


def conjugating_prefix(word, target):
    """Return u in reduced word=u target u^{-1}, or None."""
    l, r = 0, len(word)
    while r-l > 1 and word[l] == -word[r-1]:
        l += 1
        r -= 1
    return word[:l] if word[l:r] == target else None


def is_inner(F):
    u = conjugating_prefix(F[0], (1,))
    if u is None:
        return False
    G = tuple(red(inv(u),word,u) for word in F)
    p = conjugating_prefix(G[1], (2,))
    # An inner automorphism fixing a conjugates by a power of a.
    if p is None or any(abs(x) != 1 for x in p):
        return False
    return all(word == red(p,(j,),inv(p)) for j,word in enumerate(G,1))


def cyclic(word):
    while len(word) > 1 and word[0] == -word[-1]:
        word = word[1:-1]
    return word


def cyclic_equiv(u,v):
    u,v=cyclic(u),cyclic(v)
    return len(u)==len(v) and any(u==v[i:]+v[:i] for i in range(len(v)))


def abelian(word,c):
    return tuple(word.count(i)-word.count(-i) for i in range(1,c+2))


def verify_generators(c):
    I=identity(c)
    boundary=[(-j,j+1) for j in range(2,c+1)]
    boundary.append((-(c+1),1,2,-1))
    assert compose(H(c),H(c,-1)) == I
    for i in range(c):
        v=V(c,i)
        assert compose(v,V(c,i,-1)) == I
        assert compose(compose(H(c),v),H(c)) == compose(compose(v,H(c)),v)
        for j in range(c):
            assert compose(v,V(c,j)) == compose(V(c,j),v)
    for F in (H(c),*(V(c,i) for i in range(c))):
        assert all(cyclic_equiv(apply(F,p),p) for p in boundary)


def line_word(c,d):
    """Free homotopy class of a line bundle of multidegree d."""
    return cyclic(vertical(c,d)[0])


def transform_lines(c,T,F):
    """Return the simultaneous new multidegrees, or None if not all are lines."""
    result=[]
    for block in T:
        B=[]
        for d in block:
            word=apply(F,line_word(c,d))
            rank,*x=abelian(word,c)
            x=tuple(x)
            if rank!=1 or sum(x)!=sum(d) or not cyclic_equiv(word,line_word(c,x)):
                return None
            B.append(x)
        result.append(tuple(B))
    return normalize(*result)


def half_twist():
    """Swap p_1 and p_2 in the three-punctured torus (DNB)."""
    s=((1,),(2,),(2,-3,4),(4,))
    si=((1,),(2,),(4,-3,2),(4,))
    p=((-2,3),(-3,4),(-4,1,2,-1))
    assert compose(s,si)==compose(si,s)==identity(3)
    assert all(cyclic_equiv(apply(s,p[i]),p[j]) for i,j in enumerate((1,0,2)))
    assert all(sum(abelian(s[j],3)[1:])==1 for j in (1,2,3))
    return s,si


def trial_maps(c):
    """A finite family of rank/degree preserving mapping classes.

    Enumerate first the pure maps G=V^v H^k V_i V_j^{-1} H^{-k} V^{-v}
    with v in {-1,0,1}^c, k=+/-1. For c=3 also compose with the two
    half-twist representatives, on either side. No classification of MCG is
    used by this search; its successful word equalities are certificates.
    """
    z=(0,)*c
    vs=[z]+[v for v in product(range(-1,2),repeat=c) if v!=z]
    pure=[]
    for v in vs:
        for i,j in permutations(range(c),2):
            for k in (1,-1):
                d=tuple(int(t==i)-int(t==j) for t in range(c))
                G=product_of((vertical(c,v),H(c,k),vertical(c,d),
                              H(c,-k),vertical(c,mul(-1,v))),c)
                pure.append(G)
    yield from pure
    if c==3:
        for s in half_twist():
            yield s
            for G in pure:
                yield compose(G,s)
                yield compose(s,G)


def base_case(c,T):
    delta,w=invariants(*T)
    assert not any(delta)
    D=product_of((line_twist(c,d) for block in T for d in block),c)
    assert abelian(D[0],c)==(1,)+mul(-1,w)
    assert all(abelian(D[j],c)==tuple(int(k==j) for k in range(c+1))
               for j in range(1,c+1))
    geometric=is_inner(compose(D,vertical(c,w)))
    zeros=sum(all(d[j]==0 for block in T for d in block) for j in range(c))
    if min(w)>0:
        assert geometric
        return 'positive'
    if geometric:
        assert min(w)==0
        return 'geometric_boundary'
    if zeros:
        assert zeros<c
        return f'contract_{zeros}'
    return None


def classification(c):
    verify_generators(c)
    admissible=[T for T in patterns(c) if not any(invariants(*T)[0])]
    cases={T:base_case(c,T) for T in admissible}
    targets={T:k for T,k in cases.items() if k is not None}
    maps=list(trial_maps(c))
    counts=Counter()
    for T in admissible:
        if T in targets:
            counts['direct_'+targets[T]]+=1
            continue
        for F in maps:
            new=transform_lines(c,T,F)
            if new in targets:
                counts['changed_'+targets[new]]+=1
                break
        else:
            raise AssertionError(('no simultaneous line-bundle change',T))
    expected={
      2:{'direct_positive':6,'direct_geometric_boundary':4,'direct_contract_1':2,
         'changed_positive':4},
      3:{'direct_positive':9,'direct_geometric_boundary':42,'direct_contract_1':42,
         'direct_contract_2':3,'changed_positive':18,'changed_geometric_boundary':36,
         'changed_contract_1':18},
    }
    assert dict(counts)==expected[c],counts
    return counts


if __name__=='__main__':
    for c in (2,3):
        print(c,dict(classification(c)))
