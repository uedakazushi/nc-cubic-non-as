#!/usr/bin/env python3
"""Exact free-group verification for normalized nodal-cycle line-bundle patterns.

Generators (a,b_1,...,b_c) are the horizontal and based vertical loops.
H and V_i are representatives of the induced outer automorphisms of the
standard horizontal and vertical Dehn twists (see the proof in the paper).
No assertion of realizability of the multidegree patterns is made.
"""
from collections import Counter
from verify_ou_reduced_multidegrees import patterns, invariants, mul


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


def classification(c):
    verify_generators(c)
    counts=Counter()
    for T in patterns(c):
        delta,w=invariants(*T)
        if any(delta):
            continue
        D=product_of((line_twist(c,d) for block in T for d in block),c)
        # h=D^{-1} sends a in H_1 to a+sum_i w_i b_i.
        # D sends a to a-sum_i w_i b_i, and fixes b_i in H_1.
        assert abelian(D[0],c)==(1,)+mul(-1,w)
        assert all(abelian(D[j],c)==tuple(int(k==j) for k in range(c+1))
                   for j in range(1,c+1))
        geometric=is_inner(compose(D,vertical(c,w)))
        zeros=sum(all(d[j]==0 for block in T for d in block) for j in range(c))
        residues=tuple(x%3 for x in w)
        if min(w)>0:
            assert geometric
            case='positive'
        elif geometric:
            assert min(w)==0
            case='geometric_boundary'
        elif zeros:
            case=f'contract_{zeros}'
        elif c==2:
            assert min(w)<0 and sorted(w)==[-1,4]
            assert any(len(set(block))>1 for block in T)
            case='two_axis_conjugacy'
        elif len(set(residues))==1 and residues[0]!=0:
            case='three_axis_conjugacy'
        else:
            case='remaining'
        counts[case]+=1
    expected={
        2: {'positive':6,'geometric_boundary':4,'contract_1':2,
            'two_axis_conjugacy':4},
        3: {'positive':9,'geometric_boundary':42,'contract_1':42,
            'contract_2':3,'three_axis_conjugacy':18,'remaining':54},
    }
    assert dict(counts)==expected[c],counts
    return counts


if __name__=='__main__':
    for c in (2,3):
        print(c,dict(classification(c)))
